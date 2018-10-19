
#include "graph.h"
#include "auxiliary.h"

using namespace std;
using namespace Rcpp;
using namespace arma;

using IM = IntegerMatrix;
using NV = NumericVector;
using NM = NumericMatrix;

double eps = 1e-9;

Graph::Graph(NM weight_matrix, IM fixed):
    generator_(initGenerator())
{
    checks(weight_matrix, fixed);
    m_ = weight_matrix.nrow();
    n_ = weight_matrix.ncol();
    vertices_ = vector<Vertex>(m_+n_);
    for (int i = 0; i != m_ + n_; ++i)
        vertices_[i].index = i;
    // allocate memory WITHOUT calling the constructor
    edges_ = (Edge**) malloc(m_ * sizeof(Edge*));
    for (int i = 0; i != m_; ++i)
        edges_[i] = (Edge*) malloc(n_ * sizeof(Edge));
    // initialise edges
    // applies constructor directly to final address (avoids copy constructor)
    for (int i = 0; i != m_; ++i)
        for(int j = 0; j != n_; ++j)
            new (&edges_[i][j]) Edge(&vertices_[m_+j], &vertices_[i], 
                                  weight_matrix(i,j), fixed(i,j));
    // initialise initial_vertices_
    for (int i = m_; i != m_ + n_; ++i)
        if (vertices_[i].edges.size())
            initial_vertices_.push_back(&vertices_[i]);
    // check we have something...     
    if (initial_vertices_.size() == 0)
        throw invalid_argument("Matrix fully determined by specification");
}


// performs multiple sampling steps, returning a weights matrix
List Graph::sample(int nsamples, int thin, int burnin, bool sparse) 
{
    List results(nsamples);
    for (int i = 0; i != nsamples; ++i) {
        for(int j = 0; j != (thin + 1); ++j)
            sampleStep();
        if(sparse) results(i) = sparse_weight_matrix();
        else results(i) = weight_matrix();
    }
    return results;
}

// performs a single sampling step
void Graph::sampleStep() 
{
    vector<Edge*> cycle;
    int discard = sampleKernel(cycle);
    reset(cycle);
    if (discard) return;
    updateWeights(cycle, sampleDelta(cycle));
}

// reset vertex positions for new sampling step
void Graph::reset(vector<Edge *> &vec)
{
  for(auto &e : vec)
  {
      e->head()->pos = e->head()->edges.size() - 1;
      e->tail()->pos = e->tail()->edges.size() - 1;
  }
}

// print out a representation of the internal state
void Graph::summary() const
{
    Rcout << "Initial Vertices: ";
    for(const auto &v: initial_vertices_)
        Rcout << v->index + 1 << " ";
    Rcout << endl;
    for(const auto &v: vertices_)
        printVertexData(v);
}

// reconstructs weight matrix from internal datastructure
NM Graph::weight_matrix() const
{
    NM wm(m_, n_);
    for (int i = 0; i != m_; ++i)
        for (int j = 0; j != n_; ++j)
            wm(i,j) = edges_[i][j].weight();
    return wm;
}

// reconstructs a sparse matrix representation from internal datastructure
sp_mat Graph::sparse_weight_matrix() const
{
    // get total number of edges
    int nedges = 0;
    for (int i = 0; i != m_; ++i)
        nedges += vertices_[i].edges.size();
    umat locations(2, nedges);
    vec values(nedges);

    int k = 0;
    for (int i = 0; i != m_; ++i)
        for (const auto e: vertices_[i].edges)
        {
            // store locations
            locations(0,k) = i;
            locations(1,k) = e->head()->index - m_;
            // store corresponding values
            values(k) = e->weight();
            k++;
        }
    return sp_mat(locations, values, m_, n_);
}

// reconstructs fixed matrix from internal datastructure
IM Graph::fixed() const
{
    IM fm(m_, n_);
    for (int i = 0; i != m_; ++i)
        for (int j = 0; j != n_; ++j)
            fm(i,j) = edges_[i][j].fixed();
    return fm;
}

// samples a new edge uniformly
// remove the edge from consideration here-on-in
int Graph::sampleEdge(Vertex* v, vector<Edge*>& vec)
{
    // check there are free edges
    if(v->pos < 0) return 1;
    // sample and add to vector
    uniform_int_distribution<int> dist(0, v->pos);
    Edge* e = v->edges[dist(generator_)];
    e->moveBack();
    vec.push_back(e);
    return 0;
}

// samples a kernel along which we perform an update
int Graph::sampleKernel(vector<Edge*>& vec)
{
    Vertex* u0 = sampleFromVector(initial_vertices_, generator_);
    Edge* e;
    if(sampleEdge(u0, vec)) return 1;
    if(sampleEdge(vec.back()->tail(), vec)) return 1;

    while(vec.back()->head() != u0)
    {
        if(sampleEdge(vec.back()->head(), vec)) return 1;
        // attempt to close cycle
        e = &edges_[vec.back()->tail()->index][u0->index - m_];
        if(!e->fixed()) vec.push_back(e);
        else if(sampleEdge(vec.back()->tail(), vec)) return 1;
    }
    return 0;
}

// given a cycle, computes required data from the edge cases
Boundary Graph::getBoundaryData(vector<Edge*> &vec)
{
    Boundary b;
    // start with computing the range of delta
    for (int i = 0; i < vec.size() - 1; i += 2)
    {
        b.dlow = max(b.dlow, -vec[i]->weight());
        b.dup = min(b.dup, vec[i+1]->weight());
    }
    // how many edges deleted at boundaries?
    for (int i = 0; i < vec.size() - 1; i += 2)
    {
        if (vec[i].weight() + b.dlow < eps) b.nlow++;
        if (vec[i+1].weight() + b.dup < eps) b.nup++;
    }
    // likelihood at the boundaries?
    b.llow = loglDelta(vec, b.dlow);
    b.lup = loglDelta(vec, b.dup);
    return b;
}

// given a vector of edge pointers, will sample a delta from its conditional distribution
double Graph::sampleDelta(vector<Edge *> &vec)
{
    double Delta;
    Boundary b = getBoundaryData(vec);
    if (b.nlow + b.nup < 3)
    {
        // compute lambda_marg
        double lambda_marg = 0.0;
        for (int i = 0; i < vec.size() - 1; i += 2)
            lambda_marg += vec[i].lambda() - vec[i+1].lambda();
        // correct for numerical errors
        if (fabs(lambda_marg) < eps) lambda_marg = 0.;
        // calculate case probabilities
        double pall= loglDelta(vec, (b.dlow + b.dup)/2.); 
        double len = b.dup - b.dlow;
        if (lambda_marg == 0.) pall += log(len);
        else pall += log(-1./lambda_marg*(exp(-lambda_marg*(len/2.))-exp(-lambdamarg*(-len/2.))));

        // normalise for better computational stability
        double maxval = max(pall, max(b.llow, b.lup));
        pall = exp(pall - maxval);
        b.llow = exp(b.llow - maxval);
        b.lup = exp(b.lup - maxval);

        uniform_real_distribution<double> dist(0.0, 1.0);
        double u = dist(generator_)*(pall + b.llow + b.lup);

        if (pall >= u) return extExp(b, lambda_marg); // sample from extendended exponential
        if (pall + b.llow >= u) return b.dlow; // return delta low
        else return b.dup; // return delta up
    }
    else if (b.nlow > b.nup) return b.dlow;
    else if (b.nlow < b.nup) return b.dup;
    else 
    {
        uniform_real_distribution<double> dist(0.0, 1.0);
        double u = dist(generator_)*(b.llow + b.lup);
        if(b.llow >= u) return b.dlow;
        else return b.dup;
    }
}

// applies delta along a cycle
void Graph::updateWeights(vector<Edge *> &vec, double delta)
{
  for (int i = 0; i < vec.size() - 1; i += 2)
  {
      vec[i]->weight(vec[i]->weight() + delta);
      vec[i+1]->weight(vec[i+1]->weight() - delta);
  }
}


// computes log unconditional density of delta along a vector
double Graph::loglDelta(vector<Edge*> &vec, double delta)
{
    double res = 0.;
    for (int i = 0; i != vec.size(); ++i)
    {
        double val = vec[i]->weight();
        if (i % 2) val -= delta;
        else val += delta;
        if (val < eps) res += log(1. - vec[i]->p());
        else    res += log(vec[i]->p()) + log(vec[i]->lambda()) - vec[i]->lambda()*val;
    }
    return res;
}

// generates a r.v. from an extended exponential distribution
// uses the inversion method
double Graph::extExp(Boundary b, double lambda_marg)
{
    uniform_real_distribution<double> dist(0.,1.);
    double u = dist(generator_);
    if (lambda_marg = 0.) return b.dlow + u*(b.dup - b.dlow);
    else return -log((1. - u)*exp(-lambda_marg*b.dlow) + u*exp(-lambda_marg*b.dup))/lambda_marg;
}





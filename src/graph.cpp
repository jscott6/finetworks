
#include "graph.h"
#include "auxiliary.h"

using namespace std;
using namespace Rcpp;
using namespace arma;

using IM = IntegerMatrix;
using NV = NumericVector;
using NM = NumericMatrix;

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
    DeltaRange dr = getDeltaRange(cycle);
    double delta = sampleDelta(dr);
    updateWeights(cycle, delta);
}

// reset vertex positions for new sampling step
void Graph::reset(vector<Edge *> &vec)
{
  for(auto &e : vec)
  {
      e->head()->pos = e->head()->edges.size()-1;
      e->tail()->pos = e->tail()->edges.size()-1;
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

DeltaRange Graph::getDeltaRange(vector<Edge *> &vec)
{
  // compute support for Delta
  DeltaRange dr;
  for (int i = 0; i != vec.size(); ++i)
      if (i % 2) dr.up = min(dr.up, vec[i]->weight());
      else dr.low = max(dr.low, -vec[i]->weight());
  return dr;
}

// samples delta from conditional distribution
double Graph::sampleDelta(const DeltaRange& dr)
{
    uniform_real_distribution<double> dist(dr.low, dr.up);
    return dist(generator_);
}

void Graph::updateWeights(vector<Edge *> &vec, double delta)
{
  for (int i = 0; i != vec.size(); ++i)
      if (i % 2) vec[i]->weight(vec[i]->weight() - delta);
      else vec[i]->weight(vec[i]->weight() + delta);
}













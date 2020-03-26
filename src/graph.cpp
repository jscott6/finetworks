
#include "graph.h"
#include "auxiliary.h"
#include <unistd.h>

using namespace std;
using namespace Rcpp;
using namespace arma;

using IM = IntegerMatrix;
using NV = NumericVector;
using NM = NumericMatrix;

int const FAIL = 1;
int const OK = 0;

void printEdge(Edge* const edge);

Graph::Graph(NM wm, NM p, NM lambda, IM fixed, double eps):
    m_(wm.nrow()),
    n_(wm.ncol()),
    eps_(eps),
    generator_(initGenerator())
{
    for (int i = 0; i != m_; ++i) rows_.push_back(Vertex(i, vrow));
    for (int j = 0; j != n_; ++j) cols_.push_back(Vertex(j, vcol));

    edges_ = vector<vector<Edge*> >(m_, vector<Edge*>(n_));
    for (int i = 0; i != m_; ++i)
        for(int j = 0; j != n_; ++j)
            edges_[i][j] = new Edge(&rows_[i], &cols_[j], &edge_list_,
                                  wm(i,j), fixed(i,j), p(i,j), lambda(i,j));
}

// // performs multiple sampling steps, returning a weights matrix
// List Graph::sample(int nsamples, int thin, int burnin, bool sparse) 
// {
//     List results(nsamples);
//     for (int i = 0; i != nsamples; ++i) {
//         for(int j = 0; j != (thin + 1); ++j)
//             sampleStep();
//         if(sparse) results(i) = sparse_weight_matrix();
//         else results(i) = weight_matrix();
//     }
//     return results;
// }

// // performs a single sampling step
void Graph::sampleStep() 
{
    int L = m_;
    vector<Edge*> cycle;
    cycle.reserve(2 * L);
    int discard = sampleKernel(cycle, L);
    Rcout << "Discarded: " << discard << endl;

    for(auto& e : cycle)
    {
        e->ends(vrow)->visited = false;
        e->ends(vcol)->visited = false;
    }

    for (const auto& e: cycle) printEdge(e);

    if (discard) return;
    
    //updateWeights(cycle, sampleDelta(cycle));
    return;
}

// reconstructs weight matrix from internal datastructure
NM Graph::weight_matrix() const
{
    NM wm(m_, n_);
    for (int i = 0; i != m_; ++i)
        for (int j = 0; j != n_; ++j)
            wm(i,j) = edges_[i][j]->weight();

    return wm;
}


// reconstructs a sparse matrix representation from internal datastructure
sp_mat Graph::sparse_weight_matrix() const
{
    int nedges = edge_list_.size();
    umat locations(2, nedges);
    vec values(nedges);

    int k = 0;
    for (int i = 0; i != m_; ++i)
        for (const auto e: rows_[i].edges)
        {
            // store locations
            locations(0,k) = i;
            locations(1,k) = e->ends(vcol)->index;
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
            fm(i,j) = edges_[i][j]->fixed();
    return fm;
}

// Sample an edge to an unvisited vertex.
// pos is the index of the final edge that can be sampled.  
int Graph::sampleEdge(Vertex* v, vector<Edge*>& vec, int pos)
{
    if (pos == 0) return FAIL;
    uniform_int_distribution<int> dist(0, pos - 1);
    int idx = dist(generator_);
    Edge* edge = v->edges[idx];
    if (edge == vec.back()) edge = v->edges[pos];

    VertexType a = static_cast<VertexType>(1 - v->vtype);
    if (edge->ends(a)->visited) {
        edge->setPos(idx, a);
        v->edges[pos]->setPos(pos, a);
        v->edges[idx] = v->edges[pos];
        v->edges[pos] = edge;
        return sampleEdge(v, vec, pos - 1);
    }
    vec.push_back(edge);
    return OK;
}

int Graph::sampleKernel(vector<Edge*>& vec, int L)
{
    vec.push_back(sampleFromVector(edge_list_, generator_));
    Vertex *u = vec.back()->ends(vrow), *v = vec.back()->ends(vcol);
    v->visited = true;
    u->visited = true;

    for (int k = 1; k != L; k++)
    {
        if (sampleEdge(u, vec, u->edges.size() - 1)) return FAIL;
        u = vec.back()->ends(vcol);
        u->visited = true;
        if (sampleEdge(u, vec, u->edges.size() - 1)) return FAIL;
        u = vec.back()->ends(vrow);
        u->visited = true;
    }
    Edge* e = edges_[u->index][v->index];
    if (e->fixed()) return FAIL; 
    vec.push_back(e);
    return OK;
}


void printEdge(Edge* const edge) 
{
    Rcout << "Weight: " << edge->weight() << " Loc: (" << edge->ends(vrow)->index + 1 << "," << edge->ends(vcol)->index + 1 << ")" << endl;
    usleep(100000);
}

// // given a cycle, computes required data from the edge cases
// Boundary Graph::getBoundaryData(vector<Edge*> &vec)
// {
//     Boundary b;
//     // start with computing the range of delta
//     for (int i = 0; i < vec.size() - 1; i += 2)
//     {
//         b.dlow = max(b.dlow, -vec[i]->weight());
//         b.dup = min(b.dup, vec[i+1]->weight());
//     }
//     // how many edges deleted at boundaries?
//     for (int i = 0; i < vec.size() - 1; i += 2)
//     {
//         if (vec[i]->weight() + b.dlow < eps_) b.nlow++;
//         if (vec[i+1]->weight() + b.dup < eps_) b.nup++;
//     }
//     // likelihood at the boundaries?
//     b.llow = loglDelta(vec, b.dlow);
//     b.lup = loglDelta(vec, b.dup);
//     return b;
// }

// // given a vector of edge pointers, will sample a delta from its conditional distribution
// double Graph::sampleDelta(vector<Edge *> &vec)
// {
//     Boundary b = getBoundaryData(vec);


//     // Rcout << "Delta Range: (" << b.dlow << "," << b.dup << ")" << endl;
//     // Rcout << "N: (" << b.nlow << "," << b.nup << ")" << endl;
//     // Rcout << "Likelihood: (" << b.llow << "," << b.lup << ")" << endl;


//     if (b.nlow + b.nup < 3)
//     {
//         // compute lambda_marg
//         double lambda_marg = 0.0;
//         for (int i = 0; i < vec.size() - 1; i += 2)
//             lambda_marg += vec[i]->lambda() - vec[i+1]->lambda();
//         // correct for numerical errors
//         if (fabs(lambda_marg) < eps_) lambda_marg = 0.;

//         //Rcout << "lambda_marg: " << lambda_marg << endl;
//         // calculate case probabilities
//         double pall= loglDelta(vec, (b.dlow + b.dup)/2.); 
//         double len = b.dup - b.dlow;
//         if (lambda_marg == 0.) pall += log(len);
//         else pall += log(-1./lambda_marg*(exp(-lambda_marg*(len/2.))-exp(-lambda_marg*(-len/2.))));

//         // normalise for better computational stability
//         double maxval = max(pall, max(b.llow, b.lup));
//         pall = exp(pall - maxval);
//         b.llow = exp(b.llow - maxval);
//         b.lup = exp(b.lup - maxval);

//         uniform_real_distribution<double> dist(0.0, 1.0);
//         double u = dist(generator_)*(pall + b.llow + b.lup);

//         if (pall >= u) return extExp(b, lambda_marg); // sample from extendended exponential
//         if (pall + b.llow >= u) return b.dlow; // return delta low
//         else return b.dup; // return delta up
//     }
//     else if (b.nlow > b.nup) return b.dlow;
//     else if (b.nlow < b.nup) return b.dup;
//     else 
//     {
//         uniform_real_distribution<double> dist(0.0, 1.0);
//         double u = dist(generator_)*(b.llow + b.lup);
//         if(b.llow >= u) return b.dlow;
//         else return b.dup;
//     }
// }

// // applies delta along a cycle
// void Graph::updateWeights(vector<Edge *> &vec, double delta)
// {
//     //Rcout << "Delta: " << delta << endl;
//   for (int i = 0; i < vec.size() - 1; i += 2)
//   {
//       vec[i]->weight(vec[i]->weight() + delta);
//       vec[i+1]->weight(vec[i+1]->weight() - delta);
//   }
// }


// // computes log unconditional density of delta along a vector
// double Graph::loglDelta(vector<Edge*> &vec, double delta)
// {
//     double res = 0.;
//     for (int i = 0; i != vec.size(); ++i)
//     {
//         double val = vec[i]->weight();
//         if (i % 2) val -= delta;
//         else val += delta;
//         if (val < eps_) res += log(1. - vec[i]->p());
//         else res += log(vec[i]->p()) + log(vec[i]->lambda()) - vec[i]->lambda()*val;
//     }
//     return res;
// }

// // generates a r.v. from an extended exponential distribution
// // uses the inversion method
// double Graph::extExp(Boundary b, double lambda_marg)
// {
//     uniform_real_distribution<double> dist(0.,1.);
//     double u = dist(generator_);
//     if (lambda_marg == 0.) return b.dlow + u*(b.dup - b.dlow);
//     else return -log((1. - u)*exp(-lambda_marg*b.dlow) + u*exp(-lambda_marg*b.dup))/lambda_marg;
// }





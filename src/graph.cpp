
#include "graph.h"
#include "auxiliary.h"
#include <unistd.h>

using namespace std;
using namespace Rcpp;
using namespace arma;

int const FAIL = 1;
int const OK = 0;

void printEdge(Edge* const edge);

Graph::Graph(NumericMatrix const &wm, NumericMatrix const &p, NumericMatrix const &lambda, IntegerMatrix const &fixed, double eps):
    eps_(eps),
    debug_(false)
{
    for (int i = 0; i != wm.nrow(); ++i) rows_.push_back(Vertex(i, vrow));
    for (int j = 0; j != wm.ncol(); ++j) cols_.push_back(Vertex(j, vcol));

    edges_ = vector<vector<Edge*> >(rows_.size(), vector<Edge*>(cols_.size()));
    for (int i = 0; i != rows_.size(); ++i)
        for (int j = 0; j != cols_.size(); ++j)
            edges_[i][j] = new Edge(&rows_[i], &cols_[j], &edge_list_,
                                  wm(i,j), fixed(i,j), p(i,j), lambda(i,j));

    IntegerVector cycle_lengths = seq(2, min(rows_.size(), cols_.size()));
    NumericVector cycle_length_prob = 1. / as<NumericVector>(cycle_lengths);
    cycle_length_prob = cycle_length_prob / sum(cycle_length_prob);
    cycle_length_cumprob_ = cumsum(cycle_length_prob).get();
}

Graph::~Graph()
{
    for (int i = 0; i != rows_.size(); ++i)
        for(int j = 0; j != cols_.size(); ++j)
            delete edges_[i][j];
}



// performs multiple sampling steps, returning a weights matrix
List Graph::sample(int nsamples, int thin, int burnin, bool sparse) 
{
    
    List results(nsamples);

    // burnin phase
    for (int i = 0; i != nsamples * (thin + 1); i++)
        sampleStep();

    // sampling phase
    for (int i = 0; i != nsamples; ++i) {
        for(int j = 0; j != (thin + 1); ++j)
            sampleStep();
        if(sparse) results(i) = sparse_weight_matrix();
        else results(i) = weight_matrix();
    }
    return results;
}

int Graph::sampleCycleLength()
{
    double u = R::runif(0.0, 1.0);
    int i = 0;
    for (; i != cycle_length_cumprob_.size(); i++)
        if (u <= cycle_length_cumprob_[i])
            break;

    return i + 2;
}

// // performs a single sampling step
void Graph::sampleStep() 
{
    int L = sampleCycleLength();
    vector<Edge*> cycle;
    cycle.reserve(2 * L);
    int discard = sampleKernel(cycle, L);
    if(debug_) Rcout << "Discarded: " << discard << endl;

    for(auto& e : cycle)
    {
        e->ends(vrow)->visited = false;
        e->ends(vcol)->visited = false;
    }

    if(debug_)
    {
        for (const auto& e: cycle) printEdge(e);
    }

    if (discard) return;

    double delta = sampleDelta(cycle);
    if(debug_) Rcout << "Delta: " << delta << endl;
    updateWeights(cycle, delta);

    return;
}

// reconstructs weight matrix from internal datastructure
NumericMatrix Graph::weight_matrix() const
{
    NumericMatrix wm(rows_.size(), cols_.size());
    for (int i = 0; i != rows_.size(); ++i)
        for (int j = 0; j != cols_.size(); ++j)
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
    for (int i = 0; i != rows_.size(); ++i)
        for (const auto e: rows_[i].edges)
        {
            // store locations
            locations(0,k) = i;
            locations(1,k) = e->ends(vcol)->index;
            // store corresponding values
            values(k) = e->weight();
            k++;
        }
    return sp_mat(locations, values, rows_.size(), cols_.size());
}

// reconstructs fixed matrix from internal datastructure
IntegerMatrix Graph::fixed() const
{
    IntegerMatrix fm(rows_.size(), cols_.size());
    for (int i = 0; i != rows_.size(); ++i)
        for (int j = 0; j != cols_.size(); ++j)
            fm(i,j) = edges_[i][j]->fixed();
    return fm;
}

// Sample an edge to an unvisited vertex.
// pos is the index of the final edge that can be sampled.  
int Graph::sampleEdge(Vertex* v, vector<Edge*>& vec, int pos)
{
    if (pos == 0) return FAIL;
    int idx = sampleInt(pos - 1);
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
    vec.push_back(sampleFromVector(edge_list_));
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

void Graph::printRows()
{
    for(const auto &v : rows_)
        printVertex(v);
}

void Graph::printCols()
{
    for(const auto &v : cols_)
        printVertex(v);
}

// given a cycle, computes required data from the edge cases
Boundary Graph::getBoundaryData(vector<Edge*> &vec)
{
    Boundary b(vec[0], vec[1]);
    double diff;

    for (int i = 2; i < vec.size() - 1; i += 2)
    {
        diff = b.dlow + vec[i]->weight();
        if (fabs(diff) < eps_) b.nlow++;
        else if (diff < 0)
        {
            b.dlow -= diff;
            b.nlow = 1;
            b.elow = vec[i];
        }
        diff = b.dup - vec[i+1]->weight();
        if (fabs(diff) < eps_) b.nup++;
        else if ( diff > 0)
        {
            b.dup -= diff;
            b.nup = 1;
            b.elow = vec[i+1];
        }    
    }
    return b;
}

// given a vector of edge pointers, will sample a delta from its conditional distribution
double Graph::sampleDelta(vector<Edge *> &vec)
{
    Boundary b = getBoundaryData(vec);
    if(debug_) printBoundary(b);

    if (b.nlow + b.nup < 3)
    {
        double lambda_marg = 0.0;
        double pint, plow, pup;
        double len = b.dup - b.dlow;

        // compute lambda_marg, correcting for numerical errors
        for (int i = 0; i < vec.size() - 1; i += 2)
            lambda_marg += vec[i]->lambda() - vec[i+1]->lambda();
        if (fabs(lambda_marg) < eps_) lambda_marg = 0.;

        // compute unnormalised log probability of intermediate and boundary cases
        if (lambda_marg == 0.) pint = log(len);
        else pint = log((exp(lambda_marg * len / 2.) - exp(-lambda_marg * len / 2.)) / lambda_marg);
        plow = log(1 - b.elow->p()) - log(b.elow->p() * b.elow->lambda()) + lambda_marg * len / 2.;
        pup = log(1 - b.eup->p()) - log(b.eup->p() * b.eup->lambda()) - lambda_marg * len / 2.;

        // normalise log probs for stability, and exponentiate
        double maxval = max(pint, max(plow, pup));
        pint = exp(pint - maxval);
        plow = exp(plow - maxval);
        pup = exp(pup - maxval);

        double u = R::runif(0.0, pint + plow + pup);
        if (u <= pint) return randExtExp(b, lambda_marg);
        if (u <= pint + plow) return b.dlow;
        else return b.dup;
    }
    else if (b.nlow >= b.nup) return b.dlow;
    else return b.dup;
}

// applies delta along a cycle
void Graph::updateWeights(vector<Edge *> &vec, double delta)
{
    //Rcout << "Delta: " << delta << endl;
  for (int i = 0; i < vec.size() - 1; i += 2)
  {
      vec[i]->weight(vec[i]->weight() + delta);
      vec[i+1]->weight(vec[i+1]->weight() - delta);
  }
}


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

// generates a r.v. from an extended exponential distribution
// uses the inversion method
double Graph::randExtExp(Boundary b, double lambda_marg)
{
    if (lambda_marg == 0.) return R::runif(b.dlow, b.dup);
    double u = R::runif(0., 1.);
    return -log((1. - u) * exp(-lambda_marg * b.dlow) + u * exp(-lambda_marg * b.dup)) / lambda_marg;
}





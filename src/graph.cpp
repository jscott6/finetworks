
#include "graph.h"
#include "auxiliary.h"

using namespace std;
using namespace Rcpp;
using namespace arma;

using IM = IntegerMatrix;
using NV = NumericVector;
using NM = NumericMatrix;

Graph::Graph(NM weight_matrix, NM p, NM lambda, IM fixed, double eps):
    generator_(initGenerator())
{
    m_ = weight_matrix.nrow();
    n_ = weight_matrix.ncol();
    eps_ = eps;
    vertices_ = vector<Vertex>(m_ + n_);
    for (int i = 0; i != m_ + n_; ++i)
        vertices_[i].index = i;

    edges_ = vector<vector<Edge*> >(m_, vector<Edge*>(n_));
    for (int i = 0; i != m_; ++i)
        for(int j = 0; j != n_; ++j)
            edges_[i][j] = new Edge(&vertices_[m_+j], &vertices_[i], &edge_list_,
                                  weight_matrix(i,j), fixed(i,j), p(i,j), lambda(i,j));
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
// void Graph::sampleStep() 
// {
//     vector<Edge*> cycle;
//     int discard = sampleKernel(cycle);
//     if (discard)
//     {
//         reset(cycle);
//         return;
//     }
//     //for(const auto & e: cycle)
//     //    printEdgeData(*e);
//     updateWeights(cycle, sampleDelta(cycle));
//     reset(cycle);
// }

// restore vertices for new sampling step
void Graph::reset(vector<Edge *> &vec)
{
  for(auto &e : vec)
  {
      e->head()->visited = false;
      e->tail()->visited = false;
  }
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
            fm(i,j) = edges_[i][j]->fixed();
    return fm;
}

// // // Sample an edge to an unvisited vertex.  
// Edge* Graph::sampleEdge(Vertex* v, Edge* e)
// {
//     uniform_int_distribution<int> dist(0, v->edges.size() - 1);
//     Edge* edge = v->edges[dist(generator_)];
//     if (edge == e ) edge = v->edges.back();
//     if (edge->head.visited) sampleEdge(v, e);
//     else return e;
// }

// int Graph::sampleKernel(vector<Edge*>& vec)
// {
//     int L = 3;
//     vec.push_back(sampleFromVector(edge_list_, generator_);
//     for (int k = 2; k != L + 1; k++)
//     {
//         if (vec.back()->tail()->edges.size() == 1) return 1;
//         sampleEdge(vec.back()->tail(), vec.back());
//     } 

// }


// // samples a kernel along which we perform an update
// int Graph::sampleKernel(vector<Edge*>& vec)
// {
//     Vertex* u0 = sampleFromVector(initial_vertices_, generator_);
//     Edge* e;
//     if(sampleEdge(u0, vec)) return 1;
//     if(sampleEdge(vec.back()->tail(), vec)) return 1;

//     while(vec.back()->head() != u0)
//     {
//         if(sampleEdge(vec.back()->head(), vec)) return 1;
//         // attempt to close cycle
//         e = &edges_[vec.back()->tail()->index][u0->index - m_];
//         if(!e->fixed()) vec.push_back(e);
//         else if(sampleEdge(vec.back()->tail(), vec)) return 1;
//     }
//     return 0;
// }

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






#ifndef GUARD_graph
#define GUARD_graph

#include <RcppArmadillo.h>
#include <random>
#include <float.h>
#include <algorithm>
#include "edge.h"

// [[Rcpp::depends(RcppArmadillo)]]

struct Boundary
{
    double dlow, dup, llow, lup;
    unsigned int nlow, nup;
    Boundary()
      : dlow(-DBL_MAX), 
        dup(DBL_MAX), 
        nlow(0), 
        nup(0),
        llow(0.),
        lup(0.)
  {
  }
}

class Graph 
{
public:
    Graph(Rcpp::NumericMatrix weight_matrix, Rcpp::IntegerMatrix fixed);
    Rcpp::List sample(int nsamples = 10000, int thin = 10, int burnin = 5000, bool sparse = FALSE);
    void summary() const;
    void sampleStep();
    Rcpp::NumericMatrix weight_matrix() const;
    arma::sp_mat sparse_weight_matrix() const;
    Rcpp::IntegerMatrix fixed() const;
private:
    Boundary getBoundaryData(std::vector<Edge*>& vec);
    int sampleKernel(std::vector<Edge*>& vec);
    int sampleEdge(Vertex* v, std::vector<Edge*>& vec);
    void reset(std::vector<Edge *> &vec);
    void updateWeights(std::vector<Edge *> &vec, double delta);
    double sampleDelta(const DeltaRange& dr);
    double loglDelta(vector<Edge*> &vec, double delta);
    double extExp(Boundary b, double lambda_marg);
    int m_, n_;
    std::default_random_engine generator_;
    std::vector<Vertex> vertices_;
    std::vector<Vertex*> initial_vertices_;
    Edge** edges_;
};

#endif

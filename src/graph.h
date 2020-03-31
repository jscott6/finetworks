
#ifndef GUARD_graph
#define GUARD_graph

#include <RcppArmadillo.h>
#include <float.h>
#include <algorithm>
#include "edge.h"

// [[Rcpp::depends(RcppArmadillo)]]

struct Boundary
{
    double dlow, dup;
    unsigned int nlow, nup;
    Edge *elow, *eup;
    Boundary()
      : dlow(-DBL_MAX), 
        dup(DBL_MAX),
        nlow(0), 
        nup(0),
        elow(NULL),
        eup(NULL)
  {
  }
};

class Graph 
{
public:
    Graph(Rcpp::NumericMatrix const &weight_matrix, 
          Rcpp::NumericMatrix const &p, 
          Rcpp::NumericMatrix const &lambda, 
          Rcpp::IntegerMatrix const &fixed,
          double eps = 1e-9);
    ~Graph();
    Rcpp::List sample(int nsamples = 10000, int thin = 10, int burnin = 5000, bool sparse = FALSE);
    //void summary() const;
    void sampleStep();
    Rcpp::NumericMatrix weight_matrix() const;
    arma::sp_mat sparse_weight_matrix() const;
    Rcpp::IntegerMatrix fixed() const;
private:
    //Boundary getBoundaryData(std::vector<Edge*>& vec);
    int sampleCycleLength();
    int sampleKernel(std::vector<Edge*>& vec, int L);
    int sampleEdge(Vertex* v, std::vector<Edge*>& vec, int pos);
    //void updateWeights(std::vector<Edge *> &vec, double delta);
    double sampleDelta(std::vector<Edge *> &vec);
    double randExtExp(Boundary b, double lambda_marg);
    double eps_;
    Rcpp::NumericVector cycle_length_cumprob_;
    std::vector<Vertex> rows_, cols_;
    std::vector<std::vector<Edge*> > edges_;
    std::vector<Edge*> edge_list_;
};

#endif

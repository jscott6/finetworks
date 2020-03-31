
#ifndef GUARD_graph
#define GUARD_graph

#include <RcppArmadillo.h>
#include <float.h>
#include <algorithm>
#include "boundary.h"

// [[Rcpp::depends(RcppArmadillo)]]

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
    void debug(bool b) { debug_ = true; }
    void printRows();
    void printCols();
    Rcpp::NumericMatrix weight_matrix() const;
    arma::sp_mat sparse_weight_matrix() const;
    Rcpp::IntegerMatrix fixed() const;
private:
    Boundary getBoundaryData(std::vector<Edge*>& vec);
    int sampleCycleLength();
    int sampleKernel(std::vector<Edge*>& vec, int L);
    int sampleEdge(Vertex* v, std::vector<Edge*>& vec, int pos);
    void updateWeights(std::vector<Edge *> &vec, double delta);
    double sampleDelta(std::vector<Edge *> &vec);
    double randExtExp(Boundary b, double lambda_marg);
    double eps_;
    bool debug_;
    Rcpp::NumericVector cycle_length_cumprob_;
    std::vector<Vertex> rows_, cols_;
    std::vector<std::vector<Edge*> > edges_;
    std::vector<Edge*> edge_list_;
};

#endif

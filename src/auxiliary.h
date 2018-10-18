
#ifndef GUARD_auxiliary
#define GUARD_auxiliary

#include <random>
#include <RcppArmadillo.h>
#include "graph.h"
#include "construct_network.h"


std::default_random_engine initGenerator();

void checks(Rcpp::NumericMatrix weight_matrix, Rcpp::IntegerMatrix fixed);
void printEdgeData(Edge const &e);
void printVertexData(Vertex const &v);
Rcpp::List constructNetwork(Rcpp::NumericVector in_strength, 
                            Rcpp::NumericVector out_strength, Rcpp::DataFrame df);

// generic function to sample from a vector
// given some vector, return a random (uniform) element from it
// Assumes vector is non-empty
template <class T>
T sampleFromVector(const std::vector<T> &vec, std::default_random_engine &gen)
{
  std::uniform_int_distribution<int> dist(0, vec.size() - 1);
  return vec[dist(gen)];
}

#endif
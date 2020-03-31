
#ifndef GUARD_auxiliary
#define GUARD_auxiliary

#include <RcppArmadillo.h>
#include "graph.h"
#include "construct_network.h"

//void checks(Rcpp::NumericMatrix weight_matrix, Rcpp::IntegerMatrix fixed);
//void printEdgeData(Edge const &e);
//void printVertexData(Vertex const &v);
Rcpp::List constructNetwork(Rcpp::NumericVector in_strength, 
                            Rcpp::NumericVector out_strength, Rcpp::DataFrame df);


// Uniformly samples an integer from 0 up to and including n
int sampleInt(int n);

// printing functions for debugging purposes
void printEdge(Edge* const edge);
void printVertex(Vertex const &v);
void printBoundary(Boundary const &b);

// generic function to sample from a vector
// given some vector, return a random (uniform) element from it
// Assumes vector is non-empty
template <class T>
T sampleFromVector(const std::vector<T> &vec)
{
  return vec[sampleInt(vec.size() - 1)];
}



#endif

#ifndef GUARD_auxiliary
#define GUARD_auxiliary

#include <random>
#include <Rcpp.h>
#include "graph.h"

using NM = Rcpp::NumericMatrix;
using IM = Rcpp::IntegerMatrix;

std::default_random_engine initGenerator();

void checks(NM weight_matrix, IM fixed);
void printEdgeData(Edge const &e);
void printVertexData(Vertex const &v);

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
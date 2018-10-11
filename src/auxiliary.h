
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

#endif
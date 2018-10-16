
#ifndef GUARD_construct_network
#define GUARD_construct_network

#include <Rcpp.h>
#include <vector>
#include <queue>

using NV = Rcpp::NumericVector;
using NM = Rcpp::NumericMatrix;
using DF = Rcpp::DataFrame;

namespace ConstructNetwork {

  enum VertexColor{white, gray, black};

  struct Vertex {
    unsigned int distance;
    VertexColor color;
    unsigned int predecessor;
    Vertex() :
      distance(0),
      color(white),
      predecessor(0) {};
  };

  class Graph {
    std::vector<Vertex> vertices_;
    std::vector<std::vector<double> > weight_matrix_, flow_, capacity_;
    unsigned int source_, sink_;
  public:
    Graph(NV, NV, DF);
    void addEdge(unsigned int, unsigned int, unsigned int);
    bool findPath();
    unsigned int calcPathFlow();
    void updateFlow(unsigned int);
    NM constructMatrix(NV, NV);
  };
}

#endif

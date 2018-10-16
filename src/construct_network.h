
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
public:
    Graph(NV, NV, DF);
    void addEdge(unsigned int, unsigned int, unsigned int);
    bool findPath();
    unsigned int calcPathFlow();
    void updateFlow(unsigned int);
    NM constructWeightMatrix();
    IM constructFixedMatrix();
private:
    std::vector<Vertex> vertices_;
    IV tail_, head_;
    NV weights;
    std::vector<std::vector<double> > adjacency_list_, flow_, capacity_;
    unsigned int source_, sink_;
    int m_, n_;
  };
}

#endif

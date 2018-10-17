
#ifndef GUARD_construct_network
#define GUARD_construct_network

#include <Rcpp.h>
#include <vector>
#include <queue>

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
    Graph(Rcpp::NumericVector, Rcpp::NumericVector, Rcpp::DataFrame);
    void addEdge(int, int, double);
    bool findPath();
    double calcPathFlow();
    void updateFlow(double);
    Rcpp::NumericMatrix constructWeightMatrix();
    Rcpp::IntegerMatrix constructFixedMatrix();
private:
    std::vector<Vertex> vertices_;
    Rcpp::IntegerVector tail_, head_;
    Rcpp::NumericVector weights_;
    std::vector<std::vector<double> > adjacency_list_, flow_, capacity_;
    unsigned int source_, sink_;
    int m_, n_;
  };
}

#endif


#ifndef GUARD_graph
#define GUARD_graph

#include <Rcpp.h>
#include <random>
#include "edge.h"

using IM = Rcpp::IntegerMatrix;
using NV = Rcpp::NumericVector;
using NM = Rcpp::NumericMatrix;

class Graph 
{
public:
    Graph(NM weight_matrix, IM fixed);
    //Graph(NV in_strength, NV out_strength, IM fixed);
    Rcpp::List sample(int nsamples = 10000, int thin = 10, int burnin = 5000);
    void summary() const;
    void sampleStep();
    NM weight_matrix() const;
    IM fixed() const;
private:
    int sampleKernel(std::vector<Edge*>& vec);
    int sampleEdge(Vertex& v, std::vector<Edge*>& vec);
    int m_, n_;
    std::default_random_engine generator_;
    std::vector<Vertex> vertices_;
    std::vector<Vertex*> initial_vertices_;
    Edge** edges_;
};

#endif

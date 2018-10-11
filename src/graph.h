
#ifndef GUARD_graph
#define GUARD_graph

#include <Rcpp.h>
#include <random>

using IM = Rcpp::IntegerMatrix;
using IV = Rcpp::IntegerVector;

class Graph 
{
public:
    Graph(IM weight_matrix, IM fixed);
    Graph(IV in_strength, IV out_strength, IM fixed);
    Rcpp::List sample(int nsamples = 10000, int thin = 10, int burnin = 5000);
    void summary() const;
    void sampleStep();
    IM weight_matrix() const;
    IM fixed() const;
private:
    std::default_random_engine generator_;
    std::vector<Vertex> vertices_;
    Edge** edges;
};

#endif
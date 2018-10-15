
#include "graph.h"

RCPP_MODULE(graph_module)
{
    using namespace Rcpp;

    class_<Graph>("FinNet")
    .constructor<NumericMatrix, IntegerMatrix>("Constructs network from a weight matrix and a Fixed Matrix")
    .const_method("summary", &Graph::summary, "Prints Internal Datastructure to the R console")
    .property("weight_matrix", &Graph::weight_matrix, "Returns matrix of edge weights")
    .property("fixed", &Graph::fixed, "Returns matrix showing which edges are fixed")
    .method("sampleStep", &Graph::sampleStep, "Performs a single sampling step")
    .method("sample", &Graph::sample, "Performs mutliple sampling steps")
    .method("printEdge", &Graph::printEdge, "Prints edge data");
}


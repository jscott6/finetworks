
#include "graph.h"

using namespace std;
using namespace Rcpp;


Graph::Graph(IM weight_matrix, IM fixed) 
{
    int nrow = weight_matrix.nrow();
    int ncol = weight_matrix.ncol();
    vertices_ = vector<Vertex>(nrow+ncol);
    for (int i = 0; i != nrow + ncol; ++i)
        vertices_[i].index = i;
    // allocate memory WITHOUT calling the constructor
    edges_ = (Edge**) malloc(nrow * sizeof(Edge*));
    for (int i = 0; i != nrow; ++i)
        edges_[i] = (Edge*) malloc(ncol * sizeof(Edge));
    // initialise edges
    // applies constructor directly to final address (avoids copy constructor)
    for (int i = 0! i != nrow; ++i)
        for(int j = 0; j != ncol; ++j)
            new (&edges_[i][j]) Edge(&vertices_[nrow+j], &vertices_[i], 
                                  &weight_matrix_(i,j), fixed_(i,j)) 
}

Graph::Graph(IV in_strength, IV out_strength, IM fixed) 
{

}

// performs multiple sampling steps, returning a weights matrix
List Graph::sample(int nsamples, int thin, int burnin) 
{
    List results(nsamples);
    for (int i = 0; i != nsamples; ++i) {
        for(int j = 0; j != (thin + 1); ++j)
            sampleStep();
        results(i) = weight_matrix();
    }
    return results;
}

// performs a single sampling step
void Graph::sampleStep() 
{
    Rcout << "Hello World!" << endl;
}

// print out a representation of the internal state
void Graph::summary() const
{
    Rcout << "Placeholder" << endl;
}

// reconstructs weight matrix from internal datastructure
IM Graph::weight_matrix() const
{
    IM wm(nrow, ncol);
    for (int i = 0; i != nrow; ++i)
        for (int j = 0; j != ncol; ++j)
            wm(i,j) = edges_[i][j].weight();
    return wm;
}

// reconstructs fixed matrix from internal datastructure
IM Graph::fixed() const
{
    IM fm(nrow, ncol);
    for (int i = 0; i != nrow; ++i)
        for (int j = 0; j != ncol; ++j)
            fm(i,j) = edges_[i][j].fixed();
    return fm;
}

#include "graph.h"
#include "auxiliary.h"

using namespace std;
using namespace Rcpp;


Graph::Graph(NM weight_matrix, IM fixed):
    generator_(initGenerator())
{
    checks(weight_matrix, fixed);
    m_ = weight_matrix.nrow();
    n_ = weight_matrix.ncol();
    vertices_ = vector<Vertex>(m_+n_);
    for (int i = 0; i != m_ + n_; ++i)
        vertices_[i].index = i;
    // allocate memory WITHOUT calling the constructor
    edges_ = (Edge**) malloc(m_ * sizeof(Edge*));
    for (int i = 0; i != m_; ++i)
        edges_[i] = (Edge*) malloc(n_ * sizeof(Edge));
    // initialise edges
    // applies constructor directly to final address (avoids copy constructor)
    for (int i = 0; i != m_; ++i)
        for(int j = 0; j != n_; ++j)
            new (&edges_[i][j]) Edge(&vertices_[m_+j], &vertices_[i], 
                                  weight_matrix(i,j), fixed(i,j));
}

/*
Graph::Graph(IV in_strength, IV out_strength, IM fixed) 
{

}
*/

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
    for(const auto &v: vertices_)
        printVertexData(v);
}

// reconstructs weight matrix from internal datastructure
NM Graph::weight_matrix() const
{
    NM wm(m_, n_);
    for (int i = 0; i != m_; ++i)
        for (int j = 0; j != n_; ++j)
            wm(i,j) = edges_[i][j].weight();
    return wm;
}

// reconstructs fixed matrix from internal datastructure
IM Graph::fixed() const
{
    IM fm(m_, n_);
    for (int i = 0; i != m_; ++i)
        for (int j = 0; j != n_; ++j)
            fm(i,j) = edges_[i][j].fixed();
    return fm;
}

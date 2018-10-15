
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
    // initialise initial_vertices_
    for (int i = m_; i != m_ + n_; ++i)
        if (vertices_[i].edges.size())
            initial_vertices_.push_back(&vertices_[i]);
    // check we have something...     
    if (initial_vertices_.size() == 0)
        throw invalid_argument("Matrix fully determined by specification");
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

// samples a new edge uniformly
// remove the edge from consideration here-on-in
int Graph::sampleEdge(Vertex& v, vector<Edge*>& vec)
{
    // check there are free edges
    if(v.pos < 0) return -1;
    // sample and add to vector
    uniform_int_distribution<int> dist(0, v.pos);
    Edge* e = v.edges[dist(generator_)];
    e->move_back();
    vec.push_back(e);
    return 0;
}

// samples a kernel along which we perform an update
void Graph::SampleKernel(std::vector<Edge*> vec)
{

}
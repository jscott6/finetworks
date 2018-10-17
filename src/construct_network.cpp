#include "construct_network.h"
#include <stdexcept>

using namespace std;
using namespace Rcpp;
using namespace ConstructNetwork;


using NV = Rcpp::NumericVector;
using IV = Rcpp::IntegerVector;
using NM = Rcpp::NumericMatrix;
using IM = Rcpp::IntegerMatrix;
using DF = Rcpp::DataFrame;

// adds an edge to adjacency matrix, while adjusting capacity_,
void Graph::addEdge(int i, int j, double capacity){
  adjacency_list_[i].push_back(j);
  capacity_[i][j] = capacity;
}

Graph::Graph(NV out_strength, NV in_strength, DF df)
{

  m_ = out_strength.size();
  n_ = in_strength.size();
  source_ = 0;
  sink_ = m_ + n_ + 1;
  tail_ = df["tail"];
  head_ = df["head"];
  weights_ = df["weight"];

  int nvertices = sink_ + 1;
  vertices_ = vector<Vertex>(nvertices);
  adjacency_list_ = vector<vector<double> >(nvertices);
  flow_ = vector<vector<double> >(nvertices, vector<double>(nvertices,0));
  capacity_ = vector<vector<double> >(nvertices, vector<double>(nvertices,0));

  // initialise edges between source_ and rows
  for (int i = 0; i != m_; ++i)
  {
    addEdge(source_, i + 1, out_strength(i));
    addEdge(i + 1, source_, 0);
  }
  // initialise edges between columns and sink_
  for (int j = 0; j != n_; ++j)
  {
    addEdge(sink_, m_ + j + 1, 0);
    addEdge(m_ + j + 1, sink_, in_strength(j));
  }
  // initialise all edges between rows and columns
  for (int i = 0; i != m_; ++i)
    for (int j = 0; j != n_; ++j)
    {
        // simple linear search of df to ensure not fixed
        bool fxd = false;
        for (int k = 0; k != tail_.size(); ++k)
            if (i + 1 == tail_(k) && j + 1 == head_(k))
                fxd = true;
    
        if (!fxd)
        {
            addEdge(i + 1, m_ + j + 1, DBL_MAX);
            addEdge(m_ + j + 1, i + 1, 0);   
        }
    }
}


// simple implementation of breadth first search
bool Graph::findPath(){

  // initialise all vertices_
  for (auto& v: vertices_)
  {
    v.color = white;
    v.distance = 0;
    v.predecessor = 0;
  }

  vertices_[source_].color = gray;
  vertices_[source_].distance = 0;

  queue<unsigned int> Q;
  Q.push(source_);

  while (Q.size() != 0 && vertices_[sink_].color == white) {
    unsigned int u = Q.front();
    Q.pop();
    //loop through the adjacency list of vertex u
    for (auto& i: adjacency_list_[u]) {
      if (vertices_[i].color == white) {
        if (capacity_[u][i] > (flow_[u][i]-flow_[i][u])) {
          vertices_[i].color = gray;
          vertices_[i].distance = vertices_[u].distance + 1;
          vertices_[i].predecessor = u;
          Q.push(i);
          if(i==sink_) return true;
        }
      }
    }
    vertices_[u].color = black;
  }
  // there is no augmenting path to the sink_...
  return false;
}

// given a path, calculates the flow_ across that path
double Graph::calcPathFlow()
{
    int v = sink_;
    double f, edge_flow;
    
    while (v != source_)
    {
        edge_flow = capacity_[vertices_[v].predecessor][v] - flow_[vertices_[v].predecessor][v]
        + flow_[v][vertices_[v].predecessor];
    if (v == sink_) f = edge_flow;
    else f = min(edge_flow, f);
    v = vertices_[v].predecessor;
  }
  return f;
}

// method updates the flow_ in the Graph according to the path found
void Graph::updateFlow(double f)
{
  unsigned int v = sink_;
  while(v != source_)
  {
    flow_[vertices_[v].predecessor][v] += f;
    v = vertices_[v].predecessor;
  }
}

NM Graph::constructWeightMatrix()
{
    NM x(m_,n_);
    for(int i = 0; i != m_; ++i)
        for(int j = 0; j != n_; ++j)
            x(i,j) = flow_[i + 1][m_ + j + 1] - flow_[m_ + j + 1][i + 1];

    for(int k = 0; k != tail_.size(); ++k)
        x(tail_(k) - 1,head_(k) - 1) = weights_(k);
    return x;
}

IM Graph::constructFixedMatrix()
{
    IM fixed(m_, n_);
    for(int k = 0; k != tail_.size(); ++k)
        fixed(tail_(k) - 1,head_(k) - 1) = 1;
    return fixed;
}



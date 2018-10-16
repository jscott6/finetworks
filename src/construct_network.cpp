#include "construct_network.h"
#include <stdexcept>

using namespace std;
using namespace Rcpp;
using namespace ConstructNetwork;

// adds an edge to adjacency matrix, while adjusting capacity_,
void Graph::addEdge(unsigned int i, unsigned int j, unsigned int capacity){
  adjacency_matrix_[i].push_back(j);
  capacity_[i][j] = capacity;
}

Graph::Graph(IV out_strength, IV in_strength, DF df){

  // useful constants
  source_=0;
  sink_ = in_strength.size() + out_strength.size() + 1;
  unsigned int nvertices = sink_ + 1;
  IV tail = df["tail"], head = df["head"];
  NV weights = df["weights"];

  // initialising data members
  vertices_ = vector<Vertex>(nvertices);
  weight_matrix_ = vector<vector<unsigned int> >(nvertices);
  flow_ = vector<vector<unsigned int> >(nvertices, vector<unsigned int>(nvertices,0));
  capacity_ = vector<vector<unsigned int> >(nvertices, vector<unsigned int>(nvertices,0));

  // initialise edges between source_ and rows
  for(unsigned int i=0;i!=in_degree.size();++i){
    addEdge(source_, i+1, in_degree(i));
    addEdge(i+1, source_, 0);
  }
  // initialise edges between columns and sink_
  for(unsigned int j=0;j!=out_degree.size();++j){
    addEdge(sink_, in_degree.size()+j+1, 0);
    addEdge(in_degree.size()+j+1, sink_, out_degree(j));
  }
  // initialise all edges between rows and columns
  for(unsigned int i=0;i!=in_degree.size();++i){
    for(unsigned int j=0;j!=out_degree.size();++j){
      addEdge(i+1, in_degree.size()+j+1, 1);
      addEdge(in_degree.size()+j+1, i+1, 1);
    }
  }
}

bool Graph::findPath(){

  // initialise all vertices_
  for (auto& v: vertices_) {
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
    for (auto& i: adjacency_matrix_[u]) {
      if (vertices_[i].color == white) {
        if ((int)capacity_[u][i] > (int)(flow_[u][i]-flow_[i][u])) {
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
unsigned int Graph::calcPathFlow(){
  unsigned int v, f, edge_flow;
  v = sink_;

  while(v!=source_){
    edge_flow = capacity_[vertices_[v].predecessor][v] - flow_[vertices_[v].predecessor][v]
      + flow_[v][vertices_[v].predecessor];
    if(v==sink_) f = edge_flow;
    else f = min(edge_flow, f);
    v = vertices_[v].predecessor;
  }
  return f;
}

// method updates the flow_ in the Graph according to the path found
void Graph::updateFlow(unsigned int f){
  unsigned int v=sink_;
  while(v!=source_){
    flow_[vertices_[v].predecessor][v] += f;
    v=vertices_[v].predecessor;
  }
}

IntegerMatrix Graph::constructMatrix(IntegerVector in_degree, IntegerVector out_degree){

  IntegerMatrix x(in_degree.size(),out_degree.size());

  for(unsigned int i=0; i!=in_degree.size();++i)
    for(unsigned int j=0; j!=in_degree.size();++j)
      x(i,j) = flow_[i+1][in_degree.size()+j+1] - flow_[in_degree.size()+j+1][i+1];

  // validate the matrix
  vector<int> rwrong, cwrong;
  for(unsigned int i=0; i!=in_degree.size();++i){
    int total = 0;
    for(unsigned int j=0; j!=out_degree.size();++j){
      total += x(i,j);
    }
    if(total!=in_degree(i))
      rwrong.push_back(i+1);
  }

  for(unsigned int j=0; j!=out_degree.size();++j){
    int total = 0;
    for(unsigned int i=0; i!=in_degree.size();++i){
      total += x(i,j);
    }
    if(total!=out_degree(j))
      cwrong.push_back(j+1);
  }

  if(rwrong.size()>0 || cwrong.size()>0){
    throw invalid_argument("Could not reconstruct graph from degrees. Please ensure at least one such graph exists");
  }

  return x;
}
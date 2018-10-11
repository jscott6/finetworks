
#include "auxiliary.h"

using namespace std;
using namespace Rcpp;

default_random_engine initGenerator() {
  auto seed = chrono::system_clock::now().time_since_epoch().count();
  return default_random_engine(seed);
}

void checks(NM weight_matrix, IM fixed) {
  if (weight_matrix.nrow() != fixed.nrow() ||
      weight_matrix.ncol() != fixed.ncol())
    throw invalid_argument("Dimension of weight matrix and fixed matrix do not match");
  bool valid_x = true;
  bool valid_f = true;
  for (int i = 0; i != weight_matrix.nrow(); ++i) {
    for (int j = 0; j != weight_matrix.ncol(); ++j) {
      if (weight_matrix(i,j) < 0) valid_x = false;
      if (fixed(i,j) != 0 && fixed(i,j) != 1) valid_f = false;
    }
  }
  if (!valid_x)
    throw invalid_argument("All entries of weight matrix must be non-negative");
  else if (!valid_f)
    throw invalid_argument("All entries of f must be binary");
}

// samples a new edge uniformly
// remove the edge from consideration here-on-in
Edge* sample(Vertex& v, std::default_random_engine &gen)
{
    std::uniform_int_distribution<int> dist(0, v.pos);
    Edge * e = v.edges[dist(gen)];
    e->move_back();
    return e;
}

void printEdgeData(Edge const &e)
{
  Rcout << "Edge: " << e.tail() + 1 << "->" << e.head() + 1 << endl;
  Rcout << "Weight: " << e.weight() << endl;
  Rcout << endl;
}

void printVertexData(Vertex const &v)
{
  Rcout << "------------------------" << endl;
  Rcout << "VERTEX " << v.index + 1 << endl;
  Rcout << "------------------------" << endl;
  Rcout << endl;
  Rcout << "Pos: " << v.pos << endl;
  Rcout << endl;
  Rcout << setw(8) << "Edges" << endl;
  Rcout << "------------" << endl;
  Rcout << setw(8) << "Tail: ";
  for (int i = 0; i != v.edges.size(); ++i)
    Rcout << setw(2) << v.edges[i]->tail() + 1 << " ";
  Rcout << endl;
  Rcout << setw(8) << "Head: ";
  for (int i = 0; i != v.edges.size(); ++i)
    Rcout << setw(2) << v.edges[i]->head() + 1 << " ";
  Rcout << endl;
  Rcout << setw(8) << "Weight: ";
  for (int i = 0; i != v.edges.size(); ++i)
    Rcout << setw(2) << v.edges[i]->weight() << " ";
  Rcout << endl;
  Rcout << endl;
}



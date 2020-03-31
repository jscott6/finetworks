
#include "auxiliary.h"

using namespace std;
using namespace Rcpp;

int sampleInt(int n)
{
    return floor(R::runif(0., n + 1));
}

void printEdge(Edge* const edge) 
{ 
    if (edge == NULL) Rcout << "NULL" << endl;
    else Rcout << "Weight: " << edge->weight() << " Loc: (" << edge->ends(vrow)->index + 1 << "," << edge->ends(vcol)->index + 1 << ")" << endl;
    usleep(100000);
}

void printVertex(Vertex const &v)
{
  Rcout << endl;
  if(v.vtype == vrow)
    Rcout << "Row " << v.index << endl;
  else 
    Rcout << "Col " << v.index << endl;

  Rcout << "Edges: " << endl;
  for (const auto &e : v.edges)
    printEdge(e);
  Rcout << endl;
  usleep(100000);
}

void printBoundary(Boundary const &b)
{
    Rcout << "BOUNDARY DATA" << endl;
    Rcout << "[nlow, nup] : [" << b.nlow << ", " << b.nup << "]" << endl;
    Rcout << "[dlow, dup] : [" << b.dlow << ", " << b.dup << "]" << endl;
    Rcout << "____Edges____" << endl;
    Rcout << "__Low__" << endl; 
    printEdge(b.elow);
    Rcout << "__Up__" << endl; 
    printEdge(b.eup);
    Rcout << endl;
    usleep(100000);
}



// [[Rcpp::export]]
Rcpp::List constructNetwork(NumericVector out_strength, NumericVector in_strength, DataFrame df)
{
  Rcpp::List res(2);
  bool sinkfound = true;
  ConstructNetwork::Graph net(out_strength, in_strength, df);
  // adjust flow until no path remains in the residual Graph
  while (sinkfound) {
    sinkfound = net.findPath();
    net.updateFlow(net.calcPathFlow());
  }
  res(0) = net.constructWeightMatrix();
  res(1) = net.constructFixedMatrix();
  return res;
}




#include "auxiliary.h"

using namespace std;
using namespace Rcpp;

int sampleInt(int n)
{
    return floor(R::runif(0., n + 1));
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


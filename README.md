
# finetworks

This package provides an efficient MCMC sampler for perfoming Bayesian network tomography on financial networks.

We consider networks with a finite node set (the institutions), and directed edges with real-valued weights (exposures/interbank liabilities). After specifying a prior distribution on the existence of, and weights of edges in the network, the algorithm samples from the conditional distribution given the observed vertex strengths.

One can fix an arbitrary set of known edges and non-edges in the network. The methods can easily be used to sample from both monopartite and bipartite network structures.


	
Please read the [introduction to cgsampr](./vignettes/introduction.md).
	

### Installation

The package is easily installed using devtools. For example, you could use the command

~~~~
devtools::install_github("jscott6/finetworks")
~~~~

## Example Usage

~~~~
library(finetworks)
n <- 10
m <- 10
# construct monopartite weight matrix
w <- matrix(rbinom(n*m,1, prob = 0.2)*rexp(n*m, rate =1/1e3), nrow=m, ncol=n)
diag(w) <- 0
fixed <- diag(m)
p <- matrix(0.2, m, n)
lambda <- matrix(1/1e3, m, n)
# initialise FinNet object from weight matrix
net <- FinNet_from_matrix(w, p, lambda, fixed)
# sample from posterior distribution
res <- net$sample(100, 5, 0, F)
~~~~




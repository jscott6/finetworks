

FinNet_from_matrix <- function(weight_matrix, p, lambda, fixed, eps = 1e-9)
{
  fixed <- get_fixed_matrix(weight_matrix, p, lambda, fixed)
  matrix_checks(weight_matrix, p, lambda, fixed)
  new(FinNet, weight_matrix, p, lambda, fixed, eps)
}

FinNet_from_strengths <- function(out_strength, in_strength, p, lambda, fixed, eps = 1e-9){
  
  if(!is.vector(out_strength, mode = "numeric")) stop('out strengths must be a numeric vector')
  if(!is.vector(in_strength, mode = "numeric")) stop('in strengths must be a numeric vector')
  if(!all.equal(sum(out_strength), sum(in_strength))) stop('sum of marginal vectors must equate')
  # test whether df
  if(!is.data.frame(fixed)) stop("fixed must be of type data frame")
  # test for presence of correct columns
  required_names = c("tail", "head", "weight")
  if(!identical(intersect(required_names, names(fixed)), required_names)) stop('fixed must contain columns "tail", "head" and "weight"')
  # test for correct vector types
  if(!is.vector(fixed$tail, mode = "integer")) stop('column "tail" in fixed must be an integer vector')
  if(!is.vector(fixed$head, mode = "integer")) stop('column "head" in fixed must be an integer vector')
  if(!is.vector(fixed$weight, mode = "numeric")) stop('column "weight" in fixed must be a numeric vector')
  
  # ensure specified entries exist...
  if(!identical(intersect(fixed$tail, 1:length(out_strength)), fixed$tail)) stop("Please ensure coordinates specified in fixed exist")
  if(!identical(intersect(fixed$tail, 1:length(in_strength)), fixed$tail)) stop("Please ensure coordinates specified in fixed exist")
  # 
  if(!all(fixed$weight >= 0)) stop('all "weights" in fixed must be non-negative')
  
  # alter marginals by removing fixed components
  in_strength_full <- in_strength
  out_strength_full <- out_strength
  tail_aggregate <- function(y) sum(fixed$weight[fixed$tail == y])
  head_aggregate <- function(y) sum(fixed$weight[fixed$head == y])
  in_strength <- in_strength - sapply(1:length(in_strength), tail_aggregate)
  out_strength <- out_strength - sapply(1:length(out_strength), head_aggregate)
  
  # check final margins are positive
  if(!all(out_strength >=0) && all(in_strength >=0)) stop("marginal vectors less fixed entries must be non-negative")
  
  res = constructNetwork(out_strength, in_strength, fixed)
  
  # check whether the constructed network is valid
  valid = all.equal(rowSums(res[[1]]), out_strength_full) && 
          all.equal(colSums(res[[1]]), in_strength_full) && 
          all(res[[1]]>=0)
  
  if(!valid) stop("No network exists with the specified data")
  
  print(res[[2]])
  res[[2]] <- get_fixed_matrix(res[[1]], res[[2]])
  matrix_checks(res[[1]], res[[2]], lambda, fixed)
  new(FinNet, res[[1]], res[[2]])
  
}

matrix_checks <- function(weight_matrix, p, lambda, fixed)
{
  if(!is.matrix(weight_matrix) || !is.matrix(fixed) || !is.matrix(p) || !is.matrix(lambda))
    stop('Arguments must be matrices')
  dim_ok = all(dim(weight_matrix) == dim(fixed)) &&
       all(dim(weight_matrix) == dim(p)) &&
       all(dim(weight_matrix) == dim(lambda))
  if(!dim_ok)
    stop('The dimensions of all matrix arguments must match')
  if(!all(fixed %in% 0:1))
    stop('All elements of the fixed matrix must be zero or one')
  if(!all(weight_matrix >=0))
    stop('All elements of the weights matrix must be positive')
  if(!(all(p >= 0.) && all(p <= 1.)))
    stop('Probability matrix must have all values in [0,1]')
  if(all(as.logical(fixed)))
    stop('Matrix is fully determined by specification!')
}


get_fixed_matrix <- function(weight_matrix, p, lambda, fixed)
  {
    matrix_checks(weight_matrix, p, lambda, fixed)
    m <- nrow(fixed)
    n <- ncol(fixed)
    # if p = 0, must be fixed by definition
    fixed <- 1*((p == 0.) | fixed)
    fixed_prev <- fixed + 1
    # continue until convergence 
    while(any(fixed_prev != fixed))
    {
      fixed_prev <- fixed
      # check zero conditions
      fixed <- 1*(replicate(n, rowSums(weight_matrix) - rowSums(fixed*weight_matrix)==0) | fixed)
      fixed <- 1*(t(replicate(m, colSums(weight_matrix) - colSums(fixed*weight_matrix) == 0)) | fixed)
      # check if only one undetermined
      fixed <- 1*(replicate(n, rowSums(fixed) >= n - 1) | fixed)
      fixed <- 1*(t(replicate(m, colSums(fixed) >= m - 1)) | fixed) 
    }
    fixed
  }

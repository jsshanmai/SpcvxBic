#' @title faster version of Sparse Convex Biclustering Algorithm
#'
#' @description Same algorithm as \code{sparse.biADMM}. Call python code to speed up the running time. Note that please get the py file in 'inst' fold.
#' @param X The data matrix to be clustered. The rows are the samples, and the columns are the features.
#' @param nu1 A turning parameter for the augmented term for fusing-row  shrinkage
#' @param nu2 A turning parameter for the augmented term for fusing-column  shrinkage
#' @param nu3 A turning parameter for the augmented term for feature selection shrinkage
#' @param gamma_1 A regularization parameter for fusing-row shrinkage
#' @param gamma_2 A regularization parameter for fusing-column shrinkage
#' @param gamma_3 A regularization parameter for feature selection shrinkage
#' @param feature.weight The adaptive weight
#' @param m m-nearest-neighbors in the weight function
#' @param phi The parameter phi in the weight function
#' @param niters Max iteration times
#' @param tol Stopping criterion
#' @param output When output = 1, print the results at each iteration. No print when output equals other value.

#' @return A list of results, containing matrix of A, v, z, lambda1, and lambda2
#' @export
#'
#' @examples
#' \donttest{
#' # generate dataset
#' set.seed(123)
#' X = data_gen(n = 100, true_p=40, p = 80)
#' # set parameters
#' nu1 = nu2 = nu3 = gamma_1 = gamma_2 = gamma_3 = 0.1
#' m = 5
#' phi = 0.5
#' # sparse.biADMM algorithm
#' res2 = sparse.biADMM.speed(X,
#' nu1, nu2, nu3, gamma_1, gamma_2, gamma_3,
#' feature_weight = rep(1,ncol(X)),
#' m, phi, niter = 10, tol = 0.0001, output = 0)
#' dim(res2$A)
#' }

sparse.biADMM.speed = function(X, nu1, nu2, nu3,
                               gamma_1, gamma_2, gamma_3,
                               feature_weight = rep(1,ncol(X)),
                               m = 5, phi=0.5,niter = 1000,tol = 0.1,output = 1){

  requireNamespace("cvxclustr")
  requireNamespace("reticulate")


  path <- paste("./R", "SBC_ADMM_python.py", sep="/")
  source_python(path)

  n <- dim(X)[1]
  p <- dim(X)[2]

  k_row <- m; k_col <- m
  w_row <- cvxclustr::kernel_weights(t(X), phi)
  w_col <- cvxclustr::kernel_weights(X, phi)
  w_row <- knn_weights(w_row, k_row, n)
  w_col <- knn_weights(w_col, k_col, p)
  w_row <- w_row/sum(w_row)
  w_col <- w_col/sum(w_col)
  w_row <- w_row/sqrt(p)
  w_col <- w_col/sqrt(n)

  w_l <- w_row
  u_k <- w_col
  w_l <- matrix(w_l, length(w_l),1)
  u_k <- matrix(u_k, length(u_k),1)

  feature_weight <- feature_weight/ sum(feature_weight)/sqrt(p)
  feature_weight <- matrix(feature_weight,length(feature_weight),1)


  res <- SCB_ADMM_python(X, nu1, nu2, nu3, gamma_1, gamma_2,gamma_3,
                       w_l, u_k, feature_weight, niter, tol, output = output)

  result <- list(A = res[[1]], v = res[[2]], z = res[[3]], g = res[[4]],
                 lambda_1 = res[[5]], lambda_2 = res[[6]], lambda_3 = res[[7]], iters = res[[8]])
  return(result)
}

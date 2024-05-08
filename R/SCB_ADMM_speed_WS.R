#' @title SCB-ADMM-WS：SCB Algorithm with Warm-Start
#'
#' @description Similar algorithm as \code{sparse.biADMM.speed}, add initial values for A,v,g,z to speed up the running time.
#'     Note that A,v,g,z needs to be filled in at the same time.
#'
#' @param X The data matrix to be clustered. The rows are the samples, and the columns are the features.
#' @param A Initial value of A. Using the result of previous grid point as Warm-Start to speed up. Note that A,v,g,z needs to be filled in at the same time.
#' @param nu1 A regularization parameter for row shrinkage
#' @param nu2 A regularization parameter for column shrinkage
#' @param nu3 A parameter for sparsity
#' @param v Initial value of v
#' @param z Initial value of z
#' @param g Initial value of g
#' @param gamma_1 A regularization parameter for row shrinkage
#' @param gamma_2 A regularization parameter for column shrinkage
#' @param gamma_3 A parameter for sparsity
#' @param m m-nearest-neighbors in the weight function
#' @param phi The parameter phi in the weight function
#' @param niters Max iteraion times
#' @param tol Stopping criterion
#' @param output When output = 1, print the results at each iteration. No print when output equals other value.
#'
#'
#' @return A list of results, containing matrix of A, v, z, lambda1, lambda2, lambda_3 and iters
#' @export
#'
#' @examples
#' \donttest{
#' #Note that please get the py file in 'inst' fold.
#' # generate dataset
#' set.seed(123)
#' X = data_gen(n = 100, p = 80 , true_p = 40)
#' # set parameters
#' nu1 = nu2 = nu3 = gamma_1 = gamma_2 = gamma_3 = 0.1
#' m = 5
#' phi = 0.5
#' # Without fill A,v,g,z
#' res1 = sparse.biADMM.speed.WS(X, A= NULL,
#' nu1, nu2, nu3,v=NULL ,z=NULL ,g=NULL ,
#' gamma_1, gamma_2, gamma_3,
#' feature_weight = rep(1,ncol(X)),
#' m , phi,niter = 2000,tol = 1e-6,output = 0)
#' # change some parameters, use previous results as initial values
#' gamma_1 = gamma_2 = 0.5
#' res2 = sparse.biADMM.speed.WS(X,
#' A= res1$A, nu1, nu2, nu3,v=res1$v ,
#' z=res1$z ,g=res1$g ,
#' gamma_1, gamma_2, gamma_3,
#' feature_weight = rep(1,ncol(X)),
#' m , phi,niter = 2000,tol = 1e-6,output = 0)
#' #A rough comparison of computational efficiency
#' print (res1$iters)
#' print (res2$iters)
#' dim(res2$A)
#' }

sparse.biADMM.speed.WS = function(X, A= NULL, nu1, nu2, nu3,v=NULL ,z=NULL ,g=NULL ,
                               gamma_1, gamma_2, gamma_3,
                               feature_weight = rep(1,ncol(X)),
                               m = 5, phi=0.5,niter = 1000,tol = 0.1,output = 1){

  requireNamespace("cvxclustr")
  requireNamespace("reticulate")
  path <- paste("./R", "SBC_ADMM_WS.py", sep="/")
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

  res <- SCB_ADMM_python_WS(X,A,v, z, g, nu1, nu2, nu3, gamma_1, gamma_2, gamma_3,
                            w_l, u_k, feature_weight,  tol, niter,output = output)

  result <- list(A = res[[1]], v = res[[2]], z = res[[3]], g = res[[4]],
                 lambda_1 = res[[5]], lambda_2 = res[[6]], lambda_3 = res[[7]],
                 gamma_1 =res[[8]],gamma_2 =res[[9]],gamma_3=res[[10]],iters = res[[11]])
  return(result)
}

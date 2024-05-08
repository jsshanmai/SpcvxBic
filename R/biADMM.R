
############################################################################################################
# please change the dic for the py file
############################################################################################################
# path <- "/Users/yaolanqiu/Desktop/NYU/rotation/admmpackage0902/biclusterADMM/biclusterADMM/inst/biADMM.python.py"

path = "./R/biADMM.python_04.py"


######################################################
# biADMM
######################################################

#' bi-ADMM: a Biclustering Algorithm for the General Model
#'
#' @param X The data matrix to be clustered. The rows are the samples, and the columns are the features.
#' @param nu1 A regularization parameter for row shrinkage
#' @param nu2 A regularization parameter for column shrinkage
#' @param gamma_1 A regularization parameter for row shrinkage
#' @param gamma_2 A regularization parameter for column shrinkage
#' @param m m-nearest-neighbors in the weight function
#' @param phi The parameter phi in the weight function
#' @param prox The proximal maps. Could calculate L1 norm, L2 norm, or L-infinity, use "l1", "l2", or "l-inf", respectively.
#' @param niters Iteraion times
#' @param tol Stopping criterion
#' @param output When output = 1, print the results at each iteration. No print when output equals other value.
#'
#' @return A list of results, containing matrix of A, v, z, lambda1, and lambda2
#' @export
#'
#' @examples
#' \donttest{
#' #Note that please get the py file in 'inst' fold.
#' # generate dataset
#' set.seed(123)
#' X = data_gen(n = 100, p = 80, true_p=0)
#' # set parameters
#' nu1 = nu2 = gamma_1 = gamma_2 = 0.1
#' m = 5
#' phi = 0.5
#' # biADMM algorithm
#' res1 = biADMM(X, nu1, nu2, gamma_1, gamma_2,
#'  m, phi, niter = 10, tol = 0.0001, output = 0)
#' dim(res1$A)
#'}



biADMM = function(X, nu1, nu2,
                  gamma_1, gamma_2, m = 5, phi=0.5,
                  prox = 'l2',
                  niters = 1000,tol = 0.1,output = 1){

  requireNamespace("reticulate")
  #requireNamespace(cvxbiclustr)
  requireNamespace("cvxclustr")
  requireNamespace("Matrix")
  requireNamespace("MASS")

  n <- dim(X)[1]; p <- dim(X)[2]

  n2 <- n*(n-1)/2
  p2 <- p*(p-1)/2

  elks <- elk(n,p)
  el1 <- elks$el1
  el2 <- elks$el2
  ek1 <- elks$ek1
  ek2 <- elks$ek2

  k_row <- m
  k_col <- m

  w_row <- kernel_weights(t(X), phi)
  w_col <- kernel_weights(X, phi)
  w_row <- knn_weights(w_row, k_row, n)
  w_col <- knn_weights(w_col, k_col, p)
  w_row <- w_row/sum(w_row)
  w_col <- w_col/sum(w_col)
  w_row <- w_row/sqrt(p)
  w_col <- w_col/sqrt(n)

  w_l <- w_row; u_k <- w_col

  A <- matrix(0,n,p)
  v <- matrix(0,p,n2)
  z <- matrix(0,n,p2)
  lambda_1 <- matrix(0,p,n2)
  lambda_2 <- matrix(0,n,p2)

  for(iter in 1:niters){

    A_old <- A; v_old <- v; z_old <- z; lambda_1_old <- lambda_1; lambda_2_old <- lambda_2

    # update A

    En <- diag(0:(n - 1)) + diag((n - 1):0) - matrix(1, n, n) + diag(1, n, n)
    Ep <- diag(0:(p - 1)) + diag((p - 1):0) - matrix(1, p, p) + diag(1, p, p)

    M <- diag(1,n,n) + nu1 * En

    N <- nu2 * Ep

    lv <- lambda_1+ nu1 * v
    lz <- lambda_2 + nu2 * z

    C2 <- (el1-el2) %*% t(lv)
    C3 <- lz %*% t(ek1-ek2)
    C <- X +  C2 + C3

    A <- sylvester(M,t(N),C)

    al1 <- t(A) %*% el1
    al2 <- t(A) %*% el2
    ak1 <- A %*% ek1
    ak2 <- A %*% ek2

    # update vz

    if(prox == 'l1'){

      # update v
      sigma_1 <- gamma_1 * w_l/nu1
      vtemp <- al1 - al2 - 1/nu1 * lambda_1
      v <- prox_L1_m(vtemp, sigma_1)

      # update z
      ztemp <- ak1 - ak2 - 1/nu2 * lambda_2
      sigma_2 <- gamma_2 * u_k/nu2
      z <- prox_L1_m(ztemp, sigma_2)

    }else if(prox == 'l2'){

      # update v
      sigma_1 <- gamma_1 * w_l/nu1
      vtemp <- al1 - al2 - 1/nu1 * lambda_1
      v <- prox_L2_m(vtemp, sigma_1)

      # update z
      ztemp <- ak1 - ak2 - 1/nu2 * lambda_2
      sigma_2 <- gamma_2 * u_k/nu2
      z <- prox_L2_m(ztemp, sigma_2)

    }else if(prox == 'l-inf'){

      # update v
      sigma_1 <- gamma_1 * w_l/nu1
      vtemp <- al1 - al2 - 1/nu1 * lambda_1
      v <- prox_Linf_m(vtemp, sigma_1)

      # update z
      ztemp <- ak1 - ak2 - 1/nu2 * lambda_2
      sigma_2 <- gamma_2 * u_k/nu2
      z <- prox_Linf_m(ztemp, sigma_2)

    }else{
      print('Error: please specify the norms of the proximal mapping')
      break
    }

    # update lambda
    lambda_1 <- lambda_1 + nu1 * (v - al1 + al2)

    # update lambda 2
    lambda_2 <- lambda_2 + nu2 * (z - ak1 + ak2)

    if(output == 1){
      print('iter')
      print(iter)

      print(paste('A',mean(abs(A - A_old))))
      print(paste('v',mean(abs(v - v_old))))
      print(paste('z',mean(abs(z -z_old))))
      print(paste('l',mean(abs(lambda_1 - lambda_1_old))))
      print(paste('2',mean(abs(lambda_2 - lambda_2_old))))
    }


    # whether coverage
    if(mean(abs(A - A_old)) < tol &
       mean(abs(v - v_old)) < tol&
       mean(abs(z - z_old)) < tol &
       mean(abs(lambda_1 - lambda_1_old)) < tol &
       mean(abs(lambda_2 - lambda_2_old)) <tol){
      return(list(A = A,
                  v = v,
                  z = z,
                  lambda_1 = lambda_1,
                  lambda_2 = lambda_2,
                  niters = iter))
      break
    }
  }

  if(iter == niters){

    return(list(A = A,
                v = v,
                z = z,
                lambda_1 = lambda_1,
                lambda_2 = lambda_2,
                niters = iter))
  }
}





######################################################
# biADMM speed
######################################################

#' bi-ADMM: a Biclustering Algorithm for the General Model (faster version)
#'
#' Same algorithm as \code{biADMM}. Call python code to speed up the running time.
#' @param X The data matrix to be clustered. The rows are the samples, and the columns are the features.
#' @param nu1 A regularization parameter for row shrinkage
#' @param nu2 A regularization parameter for column shrinkage
#' @param gamma_1 A regularization parameter for row shrinkage
#' @param gamma_2 A regularization parameter for column shrinkage
#' @param m m-nearest-neighbors in the weight function
#' @param phi The parameter phi in the weight function
#' @param prox The proximal maps. Could calculate L1 norm, L2 norm, or L-infinity, use "l1", "l2", or "l-inf", respectively.
#' @param niters Iteraion times
#' @param tol Stopping criterion
#' @param output When output = 1, print the results at each iteration. No print when output equals other value.
#'
#'
#' @return A list of results, containing matrix of A, v, z, lambda1, and lambda2
#' @export
#'
#' @examples
#' # generate dataset
#' set.seed(123)
#' X = data_gen(n = 100, p = 80,true_p=0)
#' # set parameters
#' nu1 = nu2 = gamma_1 = gamma_2 = 0.1
#' m = 5
#' phi = 0.5
#' # biADMM algorithm
#' res2 = biADMM.speed(X, nu1, nu2, gamma_1, gamma_2,
#'  m, phi, niter = 10, tol = 0.0001, output = 0)
#' dim(res2$A)


biADMM.speed = function(X, nu1,nu2, gamma_1, gamma_2, m,
                        phi,  prox = 'l2', niters = 10, tol = 0.1, output = 1){

  requireNamespace(reticulate)
  #require(cvxbiclustr)
  requireNamespace(cvxclustr)
  requireNamespace(Matrix)
  requireNamespace(MASS)
  requireNamespace(reticulate)


  # path <- paste(system.file(package="biclusADMM"), "biADMM_python.py", sep="/")
  source_python(path)

  n <- dim(X)[1]
  p <- dim(X)[2]

  k_row <- m; k_col <- m
  w_row <- kernel_weights(t(X), phi)
  w_col <- kernel_weights(X, phi)
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

  res <- biADMM_python(X, nu1, nu2, gamma_1, gamma_2,
                       w_l, u_k,
                       prox,
                       niters, tol, output = output)

  result <- list(A = res[[1]], v = res[[2]], z = res[[3]],
                 lambda_1 = res[[4]], lambda_2 = res[[5]], iters = res[[6]])
  return(result)
}








######################################################
# projection functions
######################################################
#' Projection of L2 norm
#'
#' @param x data matrix
#' @param tau tau
#'
#' @return px
#' @export
#'


proj_L2_m = function(x, tau){
  lv = sqrt(apply(x^2,2,sum))
  px = x
  temp = matrix(c(tau/lv), dim(x)[1], dim(x)[2], byrow = TRUE) * x
  px[,(lv > tau)] = temp[, (lv > tau)]
  return(px)
}

#' Proximal of L2
#'
#' @param x data matrix
#' @param tau tau
#'
#' @return px
#' @export

prox_L2_m = function(x, tau){
  temp1 <- 1 - tau/sqrt(apply(x^2,2,sum))
  temp1 <- ifelse(temp1 < 0, 0, temp1)
  px <- matrix(temp1,dim(x)[1],dim(x)[2],byrow = TRUE) * x
  return(px)
}

#' Projection of infinity norm
#'
#' @param x data matrix
#' @param tau tau
#'
#' @return px
#' @export

proj_Linf_m = function(x, tau){
  ax = abs(x)
  temp = matrix(tau, dim(ax)[1], dim(ax)[2], byrow = TRUE)
  ax[which(ax >= temp)] = temp[which(ax >= temp)]
  px = sign.c(ax, x)
  return(px)
}

prox_L1_m = function(x, tau){
  tau = matrix(tau, dim(x)[1], dim(x)[2], byrow = TRUE)
  px = abs(x) - tau
  px[which(px < 0)] = 0
  px = sign.c(px, x)
  return(px)
}

project_to_simplex = function(x, z){
  n = length(x)
  mu = sort(x,  decreasing = TRUE)
  print(mu)
  cumsum = mu[1]
  for(j in 2:n){
    cumsum = cumsum + mu[j]
    if((j * mu[j] - cumsum + z) <= 0){
      rho = j-1
      # print(rho)
      break
    }
    rho = j
  }
  theta = (sum(mu[1:rho]) - z) / rho
  x = ifelse((x - theta) < 0, 0, x-theta)
  return(x)
}

proj_L1_m = function(x, tau){
  px = abs(x)
  for(i in 1:dim(px)[2]){
    px[,i] = project_to_simplex(px[,i], tau[i])
  }
  px = sign.c(px, x)
  return(px)
}

prox_Linf_m = function(x, tau){
  # This function performs the proximal mapping of tau * L-infinity norm.
  # It is computed via Moreau decomposition and Projection onto
  # an L1-ball of radius tau.
  # px = x - project_L1(x,tau)
  px = matrix(1/tau, dim(x)[1], dim(x)[2], byrow = TRUE) * x
  px = x - proj_L1_m(x, tau)
  return(px)
}


#' Adjusts the sign of a vector according to specified conditions.
#'
#' @param px The condition for adjusting the sign.
#' @param x The vector whose sign needs adjustment.
#'
#' @return a
#'
#' @export

sign.c = function(px, x){
  a = ifelse(px == abs(x), x, ifelse(x > 0, abs(px), -abs(px)))
  return(a)
}

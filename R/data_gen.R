#' Simulated data set generation
#'
#' \code{data_gen} generates dataset contains uninformative features in the simulation study.
#'
#' @param seed.cluster seed.cluster is used to control the cluster assignment
#' @param seed.data seed.data is used to control the data generation given clustering structure
#' @param n number of subjects
#' @param true_p number of informative features
#' @param p number of total features
#' @param theta the standard deviation of the x, which is sampled from a normal distribution
#' @param theta_noise the standard deviation of noisy features
#' @param mu.lower lower bound on the mean of the normal distribution
#' @param mu.upper upper bound on the mean of the normal distribution
#' @param row_group number of row clusters
#' @param col_group number of informative column clusters
#'
#' @return Output is the simulated dataset
#' @export
#'
#' @examples
#' set.seed(123)
#' X = data_gen(n = 100, true_p=40, p = 80)

data_gen <- function(seed.cluster = 123, seed.data = 654, n, true_p, p, mu.lower=-10, mu.upper=10, theta = 2, theta_noise=0.1, row_group = 4, col_group = 4){

  mu <- seq(mu.lower,mu.upper,0.2)
  x <- matrix(0, n, p)

  set.seed(seed.cluster)
  mu_rc <- sample(mu, size = row_group*col_group, replace = T)
  dim(mu_rc) <- c(row_group,col_group )

  row_assign <- c()
  col_assign <- c()
  for(i in 1:n){
    row_assign <- c(row_assign, sample(1:row_group,1))
  }
  for(i in 1:true_p){
    col_assign <- c(col_assign, sample(1:col_group,1))
  }

  ################
  set.seed(seed.data)
  for(i in 1:n){
    for(j in 1:true_p){
      r <- row_assign[i]
      c <- col_assign[j]
      mu_temp <- mu_rc[r,c]
      x[i,j] <- rnorm(1,mu_temp,theta)
    }
  }
  set.seed(seed.data*2)
  x[,(true_p+1):p] = matrix(rnorm(n*(p-true_p),0,theta_noise),n,p-true_p)

  rownames(x) <- paste('row', row_assign, sep = '')
  colnames(x) <- c(paste('col', col_assign, sep = ''),rep("noisy",p-true_p))

  return(x)
}

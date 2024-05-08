#' Cluster_assignment: get the cluster labels after ADMM convex-biclustering algorithm
#'
#' The function: tri2vec, vec2tri, kernel_weights, knn_weights, get_subgroup_means_full, biclust_smooth,
#' create_adjacency, find_clusters are from the cvxbiclustr package.


tri2vec <- function(i,j,n) {
  return(n*(i-1) - i*(i-1)/2 + j -i)
}

vec2tri <- function(k,n) {
  i <- ceiling(0.5*(2*n-1 - sqrt((2*n-1)^2 - 8*k)))
  j <- k - n*(i-1) + i*(i-1)/2 + i
  return(as.matrix(cbind(i,j)))
}



#' Compute Gaussian Kernel Weights
#'
#' \code{kernel_weights} computes Gaussian kernel weights given a data matrix \code{X} and a scale parameter \code{phi}. Namely,
#' the lth weight \code{w[l]} is given by
#' \deqn{
#' w[l] = exp(-phi ||X[,i]-X[,j]||^2)
#' }, where the lth pair of nodes is (\code{i},\code{j}).
#' @param X The data matrix to be clustered. The rows are the features, and the columns are the samples.
#' @param phi The nonnegative parameter that controls the scale of kernel weights
#' @return A vector \cite{w} of weights for convex clustering.
#'
kernel_weights <- function(X,phi=1) {
  requireNamespace("biADMM")
  requireNamespace("cvxclustr")
  storage.mode(X) <- "double"
  p <- as.integer(nrow(X))
  n <- as.integer(ncol(X))
  phi <- as.double(phi)
  w <- double(n*(n-1)/2)
  sol <- .C('kernel_weights',X=X,p=p,n=n,phi=phi,w=w)
  return(weights=as.vector(sol$w))
}


#' "Thin" a weight vector to be positive only for its k-nearest neighbors
#'
#' \code{knn_weights} takes a weight vector \code{w} and sets the ith
#' component \code{w[i]} to zero if either of the two corresponding nodes
#' is not among the other's \code{k} nearest neighbors.
#'
#' @param w A vector of nonnegative weights. The ith entry \code{w[i]} denotes the weight used between the ith pair of centroids. The weights are in dictionary order.
#' @param k The number of nearest neighbors
#' @param n The number of data points.
#' @return A vector \cite{w} of weights for convex clustering.

knn_weights <- function(w,k,n) {
  i <- 1
  neighbors <- tri2vec(i,(i+1):n,n)
  keep <- neighbors[sort(w[neighbors],decreasing=TRUE,index.return=TRUE)$ix[1:k]] #find the k nearst
  for (i in 2:(n-1)) {
    group_A <- tri2vec(i,(i+1):n,n)
    group_B <- tri2vec(1:(i-1),i,n)
    neighbors <- c(group_A,group_B)
    knn <- neighbors[sort(w[neighbors],decreasing=TRUE,index.return=TRUE)$ix[1:k]]
    keep <- union(knn,keep)
  }
  i <- n
  neighbors <- tri2vec(1:(i-1),i,n)
  knn <- neighbors[sort(w[neighbors],decreasing=TRUE,index.return=TRUE)$ix[1:k]]
  keep <- union(knn,keep)
  if (length(keep) > 0)
    w[-keep] <- 0
  return(Matrix(data=w,ncol=1,sparse=TRUE))
}

#' Get Subgroup Means
#'
#' \code{get_subgroup_means_full} computes the subgroup means.
#'
#' @param X Data matrix
#' @param clusters_row Row cluster object
#' @param clusters_col Column cluster object
#' @export
#'
get_subgroup_means_full <- function (X, clusters_row, clusters_col){
  Y <- X
  num_clusters_row <- length(clusters_row$size)
  num_clusters_col <- length(clusters_col$size)
  M <- matrix(NA, num_clusters_row, num_clusters_col)
  for (i in 1:num_clusters_row) {
    ix_row <- which(clusters_row$cluster == i) #verify the rows belong to this subgroup
    for (j in 1:num_clusters_col) {
      ix_col <- which(clusters_col$cluster == j)
      M[i, j] <- mean(Y[ix_row, ix_col], na.rm = TRUE) # M[i, j] is the mean of one block
    }
  }
  return(M)
}

#' Bicluster Smooth
#'
#' \code{biclust_smooth} computes the bicluster estimates given row and column partitions
#'
#' @param X Data matrix
#' @param clusters_row Row cluster object
#' @param clusters_col Column cluster object
#' @export
#'
biclust_smooth <- function(X,clusters_row,clusters_col) {
  p <- nrow(X); n <- ncol(X)
  Y <- matrix(NA,p,n)
  M <- get_subgroup_means_full(X,clusters_row,clusters_col)
  num_clusters_row <- length(clusters_row$size)
  num_clusters_col <- length(clusters_col$size)
  for (i in 1:num_clusters_row) {
    ixi <- which(clusters_row$cluster == i)
    for (j in 1:num_clusters_col) {
      ixj <- which(clusters_col$cluster == j)
      Y[ixi,ixj] <- M[i,j]
    }
  }
  return(Y)
}

#' Create adjacency matrix from V
#'
#' \code{create_adjacency} creates an n-by-n sparse adjacency matrix from the matrix of centroid differences.
#'
#' @param V Matrix of centroid differences
#' @param Phi Edge-incidence matrix
#' @import Matrix
#' @export
#'
create_adjacency <- function(V,Phi) {
  differences <- apply(V,2,FUN=function(x) {norm(as.matrix(x),'f')})
  connected_ix <- which(differences == 0)
  n <- ncol(Phi)
  m <- length(connected_ix)
  A <- Matrix(0, nrow = n, ncol = n, sparse = TRUE)

  if (m > 0) {
    ix <- integer(m)
    jx <- integer(m)
    for (i in 1:m) {
      ix[i] <- which(Phi[connected_ix[i],]==1)
      jx[i] <- which(Phi[connected_ix[i],]==-1)
    }
    A[(jx-1)*n + ix] <- 1
  }
  return(as(A,"generalMatrix"))
  #return(as(A,"dgCMatrix")) ## convert the matrix class to dgCMatrix, and A has an original form ddiMatrix which
                            ## could be used by graph.adjacency
}

#' Find clusters
#'
#' \code{find_clusters} uses breadth-first search to identify the connected components of the corresponding
#' adjacency graph of the centroid differences vectors.
#'
#' @param A adjacency matrix
#' @export
#' @import igraph
find_clusters <- function(A) {
  G <- graph.adjacency(A, mode = 'upper')
  n <- nrow(A)
  node_seen <- logical(n) # true false
  cluster <- integer(n)
  k <- 1
  for (i in 1:n) {
    if (!node_seen[i]) {
      connected_set <- graph.bfs(G, root=i, unreachable = FALSE)$order
      node_seen[connected_set] <- TRUE
      cluster[connected_set] <- k
      k <- k + 1
    }
  }
  nClusters <- k - 1
  size <- integer(nClusters)
  for (j in 1:nClusters) {
    size[j] <- length(which(cluster == j))
  }
  return(list(cluster=cluster, size=size))
}
#'
#' \code{cluster_assign}
#' @param X Data matrix
#' @param k the number of nearest-neighbors
#' @param result the list, output of bi_ADMM function.
#' @export
#'
cluster_assign <- function(X, k, result,phi=0.5){
  phi <- phi
  k_row = k_col = k
  n <- nrow(X)
  p <- ncol(X)

  # set weight
  w_row <- kernel_weights(t(X), phi)
  w_col <- kernel_weights(X, phi)
  w_row <- knn_weights(w_row, k_row, n)
  w_col <- knn_weights(w_col, k_col, p)
  w_row <- w_row/sum(w_row)
  w_col <- w_col/sum(w_col)
  w_row <- w_row/sqrt(p)
  w_col <- w_col/sqrt(n)

  nrow <- dim(X)[1]; ncol <- dim(X)[2]
  A <- result$A
  V1 <- result$v[,which(w_row!=0)]
  Z1 <- result$z[,which(w_col!=0)]

  w <- as.vector(w_row)
  P <- which(w!=0)
  P <- vec2tri(P,nrow)
  E <- matrix(0,nrow,dim(P)[1])
  for(i in 1:dim(P)[1]){
    E[P[i,1],i] <- 1
    E[P[i,2],i] <- -1
  }
  E1 <- E

  ncol <- dim(X)[2]
  w <- as.vector(w_col)
  P <- which(w!=0)
  P <- vec2tri(P,ncol)
  E <- matrix(0,ncol,dim(P)[1])
  for(i in 1:dim(P)[1]){
    E[P[i,1],i] <- 1
    E[P[i,2],i] <- -1
  }
  E2 <- E

  clusters_row <- find_clusters(create_adjacency(V1,t(E1)))
  clusters_col <- find_clusters(create_adjacency(Z1,t(E2)))

  M <- biclust_smooth(X,clusters_row,clusters_col)
  return(list(M = M, clusters_row = clusters_row, clusters_col = clusters_col))
}


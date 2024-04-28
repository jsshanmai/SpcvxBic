# function to generate cluster label for biclustering problem using row and col names

#' @title generate cluster label for bicluster
#'
#' @param X The data matrix to be clustered, containing rownames and colnames.
#'
#' @return bicluster labels
#' @export
#'
#' @examples
#'
#' x <-matrix(c(1,2,3,4),2,2)
#' rownames(x) <-c('row1','row2')
#' colnames(x) <-c('col1','col2')
#' bicluster.label(x)

bicluster.label = function(X){
  rname = rownames(X)
  cname = colnames(X)
  col1 = rep(rname, dim(X)[2])
  col2 = rep(cname, each = dim(X)[1])
  col3 = paste(col1, col2)
  groups = data.frame(col1 = col1, col2 = col2, col3 = col3)
  groups$num = NULL
  col3.name = names(table(col3))
  for(i in 1:dim(groups)[1]){
    groups$num[i] =  which(col3.name == col3[i])
  }
  return(groups)
}

#' @title generate cluster label for cluster using row or col names
#'
#' @param X The data matrix to be clustered
#' @param labels using row or col names
#'
#' @return cluster labels
#' @export
#'
#' @examples
#' x.rowlab<-cluster.label(x,labels = 'row')$num
#' x.collab<-cluster.label(x)$num
#' # adjustedRandIndex(x.rowlab,x.collab)
#'
cluster.label = function(X ,labels ='col'){
  if (labels == 'row'){
    rname = rownames(X)
    col1 = rep(rname, dim(X)[2])
    groups = data.frame(col1 = col1)
    groups$num = NULL
    col3 =col1
    col3.name = names(table(col3))
    for(i in 1:dim(groups)[1]){
      groups$num[i] =  which(col3.name == col3[i])
    }
  }

  if (labels == 'col'){
    cname = colnames(X)
    col2 = rep(cname, each = dim(X)[1])
    groups = data.frame(col2 = col2)
    groups$num = NULL
    col3 =col2
    col3.name = names(table(col3))
    for(i in 1:dim(groups)[1]){
      groups$num[i] =  which(col3.name == col3[i])
    }
  }

  return(groups)
}


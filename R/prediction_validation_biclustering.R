#' Predict new data for biclustering problem
#'
#' report Rand index and variable stability in validate dataset
#' @param X (n by p)
#' @param clusters_row cluster assignment for each subject
#' @param clusters_col cluster assignment for each feature
#' @param data_validate validate dataset
#'
predict_bi_validate <- function(X,clusters_row,clusters_col, X_validate,label_valdiate=Null,label.known=1){

  size.row <- table(clusters_row)
  size.col <- table(clusters_col)
  n.cluster.row <- length(size.row) #类个数
  n.cluster.col <- length(size.col)

  p1  <- ncol(X_validate)
  n1  <- nrow(X_validate)

  # assign row clusters
  center.row <- matrix(0, n.cluster.row, p1)
  dist.row <- matrix(0, n.cluster.row, n1)
  for(i.cluster.row in 1:n.cluster.row){
    center0.row <- apply( as.matrix( X[clusters_row == i.cluster.row,]), 2, mean) # 同一类的中心算出来
    dist0.row <- apply(X_validate, 1, function(x) sqrt( sum( (x - center0.row)^2 ) ) ) # sd
    center.row[i.cluster.row, ] <- center0.row
    dist.row[i.cluster.row, ] <- dist0.row
  }
  cluster.validate.row <- apply(dist.row, 2, which.min)

  # assign col clusters
  center.col <- matrix(0, n1,n.cluster.col)
  dist.col <- matrix(0, p1,n.cluster.col)
  for(i.cluster.col in 1:n.cluster.col){
    center0.col <- apply( as.matrix( X[,clusters_col == i.cluster.col]), 1, mean)
    dist0.col <- apply(X_validate, 2, function(x) sqrt( sum( (x - center0.col)^2 ) ) )
    center.col[,i.cluster.col] <- center0.col
    dist.col[,i.cluster.col ] <- dist0.col
  }
  cluster.validate.col <- apply(dist.col, 1, which.min)

  # make combined row and col names
  dat3 = X_validate
  rownames(dat3) = cluster.validate.row
  colnames(dat3) = cluster.validate.col

  validate.X.groups = bicluster.label(dat3)
  # with known cluster label, calculate rand index; otherwise return predicted clusters
  if (label.known ==1){
    adj.rand<- adjustedRandIndex(validate.X.groups$num, label_valdiate)
    return(data.frame(adj.rand=adj.rand))
  } else {
    return(validate.X.groups)
  }
}


# prepare for python version: calculate w_l and u_k

kernal.knn = function(X){
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
  
  result <- list(w_l =w_l,u_k=u_k)
  return(result)
}
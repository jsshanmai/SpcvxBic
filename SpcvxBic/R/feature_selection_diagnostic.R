#' Evaluate diagnostic statistics for estimated clusters
#'
#' @param X Data matrix
#' @param A Biclustering result
#' @param tol threshold of noisy variable
#'
#' @return diagnostic statistics
#' @export
#'
#' @examples
#' set.seed(123)
#' X = data_gen(n = 100, true_p=40, p = 80)
#' #Pretend A is the Biclustering result
#' set.seed(111)
#' A = data_gen(n = 100, true_p=40, p = 80)
#' feature.selection.diagnostic(X, res1$A)

feature.selection.diagnostic <- function(X, A,tol = 1e-6)
{
  require(caret)
  true.feature <- colnames(X) != "noisy" # select not noisy feature
  # use L2 norm for column estimates to decide whether it is a noisy variable
  estimated.feature <- apply(A, 2, sd) > tol

  temp <- caret::confusionMatrix(factor(estimated.feature,levels=c("FALSE","TRUE")),
                          factor(true.feature,levels=c("FALSE","TRUE")), positive = "TRUE",mode = "sens_spec")
  return(temp$byClass)
}


# define a new function diag sensetivity and specificity, when true_p/p is known
feature.selection.diagnostic.quatile<-function(X, A,p,true_p,method =sd){
  require(caret)
  true.feature <- colnames(X) != "noisy" # select not noisy feature
  tol = quantile(apply(A, 2, method),1-true_p/p) # use the label ratio
  estimated.feature <- apply(A, 2, method) > tol
  temp <- caret::confusionMatrix(factor(estimated.feature,levels=c("FALSE","TRUE")),
                                 factor(true.feature,levels=c("FALSE","TRUE")), positive = "TRUE",mode = "sens_spec")
  return(temp$byClass[c(1,2)])
}

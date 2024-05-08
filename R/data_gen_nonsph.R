# 产生非球面数据的code
data_gen_nonsph <- function(n=200,n_threshold=100, p =400, true_p=200,seed.data = 654) {
    set.seed(seed.data)
  x <- matrix(0, n, p)
  epu <-rnorm(n,0,0.2)
  theta <- runif (500, min = 0, max = pi)
  for (j in 1:true_p) {

    if (j %% 2 == 1) {  # 奇数
      #colnames(x) <-paste('col', 1, sep = '')
      for (i in 1:n){
     x[i,j]=-2*ifelse(i <= n_threshold, 1, 0)+5*cos(theta[i]+pi*ifelse(i >= n_threshold, 1, 0))+epu[i]
    }}
    if (j %% 2 == 0) {
      #colnames(x) <-paste('col', 2, sep = '')
      for (i in 1:n){
      x[i,j]=5*ifelse(i <= n_threshold, 1, 0)+5*sin(theta[i]+pi*ifelse(i >= n_threshold, 1, 0))+epu[i]
    }}
  }
  x[,(true_p+1):p] = matrix(rnorm(n*(p-true_p),0,1),n,p-true_p)

  rownames(x) <- c(paste('row', c(rep(1,n_threshold),rep(2,n-n_threshold)), sep = '') )
  colnames(x) <- c(paste('col', rep(1,true_p), sep = ''),rep("noisy",p-true_p))
  return(x)
}
data <-data_gen_nonsph()
plot(data)



heatmap(data)
res3 = sparse.biADMM.speed(scale.X,  nu1, nu2, nu3,
                           gamma_1=gamma_1, gamma_2=gamma_2, gamma_3=gamma_3,
                           feature_weight = feature_weight,
                           m , phi,niter = 2000,tol = tol,output = 0)

heatmap(data)


#run cobra
scale.X <-data
bi.X.groups = bicluster.label(scale.X)
nGamma <- 8
gammaSeq <- 10**seq(0,5,length.out=nGamma)
## Generate solution path
phi <- 1
k <- 5#10
wts <- gkn_weights(scale.X,phi=phi,k_row=k,k_col=k)
w_row <- wts$w_row
w_col <- wts$w_col
E_row <- wts$E_row
E_col <- wts$E_col

## Connected Components of Row and Column Graphs
wts$nRowComp
wts$nColComp
#sol <- cobra(scale.X,E_row,E_col,w_row,w_col,gammaSeq)
sol <- cobra_validate(scale.X,E_row,E_col,w_row,w_col,gammaSeq,fraction=0.01)



# SpcvxBic
Identifies a bicluster of a data matrix, clustering observations and features simultaneously, and apply sparsity to features. 'SpcvxBic' can be used for high-dimensional data, identifies uninformative features. The method is based on Bi-ADMM algorithm. For a detailed introduction to this algorithm, please refer to the author's related papers.

# Note
Download the compressed file with the extension ".gz" to locally install the SpcvxBic package in RStudio.
main function is "SCB_ADMM_speed.R". Please use the corresponding Python file and place it in the working directory to smoothly call Python functions and accelerate program computations. "SCB_ADMM_speed_WS.R" can use previous grid search result {A,v,g,z} to accelerate calculate speed.

# Example
Full code of simulation were provided with supplyment file of paper "sparse convex bicluster".
Here we show the example code of simulation, compared with biADMM algorithm

#05/08/2024 at Beijing Normal University Zhuhai Campus
pac <- c("survival","plyr","ggplot2","reshape2","phyloseq",'dirmult',"microbiome","vegan","e1071","caret","pROC",
         "fossil","cvxclustr","cvxbiclustr","doParallel","foreach","mclust","Matrix","MASS","reticulate")
.tmp <- lapply(pac, library, character.only = T)
setwd('C:/Users/matebook 14/Desktop/~/Re_ sparse convex biclustering')

source("./R/data_gen.R")
source("./R/sylvester.R")
source("./R/SCB_ADMM.R")
source("./R/SCB_ADMM_speed.R")
source("./R/util.R")
source("./R/cluster assignment.R")
source('./R/prediction_validation_biclustering.R')
source('./R/bicluster.label.R')
source('./R/feature_selection_diagnostic.R')
source('./R/biADMM.R')
source("./R/SCB_ADMM_speed_WS.R")
#evaluate ARI of algo with mean and sd
MSD <- function(x) {
  if (is.list(x)) {
    x <- unlist(x)
  }
  mean_val <- round(mean(x),2)
  sd_val <- round(sd(x),2)
  paste("mean is", mean_val, ", sd is", sd_val)
}

feature.selection.diagnostic.ROC <-function(X,A,length.out =2000){
  ROC <-matrix(NA,length.out,2)
  tol.seq <-seq(min(A),max(A),length.out =length.out)
  for ( i in 1:length.out)   {
    ROC[i,] <-as.numeric(feature.selection.diagnostic(X=X,A=A,tol = tol.seq[i])[c(1,2)] )
  }
  return(ROC)
}

#define a new function diag sensetivity and specificity
feature.selection.diagnostic.quatile<-function(X, A,p,true_p,method =sd,tol = 1e-6){
  require(caret)
  true.feature <- colnames(X) != "noisy" # select not noisy feature
  tol = quantile(apply(A, 2, method),1-true_p/p) # use the label ratio
  estimated.feature <- apply(A, 2, method) > tol
  temp <- caret::confusionMatrix(factor(estimated.feature,levels=c("FALSE","TRUE")),
                                 factor(true.feature,levels=c("FALSE","TRUE")), positive = "TRUE",mode = "sens_spec")
  return(temp$byClass[c(1,2)])
}
#feature.selection.diagnostic.quatile (X=scale.X,A=best.bi.A.list[[ii]]* norm.X.f.vector[ii],p=p,true_p = true_p,method =mean)


nu1 = nu2 = nu3 = 1e+4 # 500
m = 5
phi = 1
tol = 5e-6

#data generation set-up

n = 240
p = 400 # 200
true_p = 60 #10 is good may 15
theta = 1 # 1 is good
theta_noise = 3# 1 is good,2^(1/2) is great
row_group = 4 #行聚成几类
col_group = 4 #列聚成几类，col对应特征
mu.lower = -5#-3
mu.upper = 5#3

#repetition
rep.num = 10

#check setting
X = data_gen(seed.cluster = 123, seed.data = 654, n=n, true_p=true_p, p=p,
             mu.lower=mu.lower, mu.upper=mu.upper,
             theta = theta, theta_noise=theta_noise, row_group = row_group, col_group = col_group)
dim(X)
#mu <-seq(-5,5,0.2);sample(mu, size = row_group*col_group, replace = T)

#heatmap example
plot.X = X
plot.X = X[,1:true_p] # true_p is the first p feature and noise otherwise

col_types = colnames(plot.X)
col_ty = as.numeric(factor(col_types))
row_types = rownames(plot.X)
row_ty = as.numeric(factor(row_types))
cols <- rainbow(4)
YlGnBu5 <- c('#ffffd9','#c7e9b4','#41b6c4','#225ea8','#081d58')
hmcols <- colorRampPalette(YlGnBu5)(256)

heatmap(plot.X,col=hmcols,labRow=NA,labCol=NA,
        ColSideCol=cols[col_ty],RowSideCol=cols[row_ty],
        main = 'Example')


#Initialize path parameters and structures
nGamma_1_2 <- 12
gammaSeq_1_2 <- exp(seq(1,7,length.out=nGamma_1_2)) #1,7
gammaSeq_1_2 <-gammaSeq_1_2[c(-1,-2,-3,-4)]
nGamma_1_2 <- length(gammaSeq_1_2)

nGamma_3 <- 12
gammaSeq_3 <- exp(seq(2,6,length.out=nGamma_3))#6
gammaSeq_3 <-gammaSeq_3[c(-1,-2,-3,-4)]
nGamma_3 <- length(gammaSeq_3)
gamma_12_3.grid <-expand.grid(gammaSeq_1_2,gammaSeq_3) # obtain a grid with gamma12_3

#get the right order to warm start, note that the previous results are great Warm Start
#if find the results (A,v,g,z) with very close parameters, iters will reduce significantly
gamma_12_3.grid <- gamma_12_3.grid[order(-gamma_12_3.grid$Var1, -gamma_12_3.grid$Var2), ]
gamma_12_3.grid$index <-seq_len(nrow(gamma_12_3.grid))
gamma_12_3.grid$ws_index <- seq_len(nrow(gamma_12_3.grid))-1

#get the gamma_3 max index and change ws_order
max3_index <-which(gamma_12_3.grid$index %% nGamma_3 == 0)+1
max3_index <-max3_index[-length(max3_index)] # delete last one
gamma_12_3.grid$ws_index[max3_index] <-gamma_12_3.grid$index[max3_index]-nGamma_3
gamma_12_3.grid
#with new gamma solution

#initialize vectors
best.bi.admm.adj.rand = numeric()
best.bi.admm.gamma.index = numeric()
best.val.bi.admm.adj.rand = numeric()
best.val.bi.admm.gamma.index = numeric()
best.bi.iters = numeric()

best.scb.admm.adj.rand = numeric()
best.scb.admm.gamma.index = numeric() #
best.val.scb.admm.adj.rand = numeric()
best.val.scb.admm.gamma.index = numeric()
best.scb.iters = numeric()

times <- numeric()

norm.X.f.vector <- numeric()

#initialize output lists
bi.admm.result.list = list()
scb.admm.result.list = list()
for (i in 1:rep.num) bi.admm.result.list[[i]]<-list()
for (i in 1:rep.num) scb.admm.result.list[[i]] <- list()
for (i in 1:rep.num) {
  for ( ii in 1:nrow(gamma_12_3.grid)) scb.admm.result.list[[i]][[ii]] <- list()
}

bi.admm.adj.rand.list = list()
bi.admm.validate.adj.rand.list = list()

best.bi.admm.gamma.list = list()
best.bi.admm.validate.gamma.list = list()

best.bi.A.list = list()
best.bi.admm.X.groups.list = list()

scb.admm.adj.rand.list = list()
scb.admm.validate.adj.rand.list = list()

best.scb.admm.gamma.list = list()
best.scb.admm.validate.gamma.list = list()

best.scb.A.list = list()
best.scb.admm.X.groups.list = list()

best.scb.fs.diagnostic.list = list()

#start simulations

for (ii in 1:rep.num){ # 重复50 次

  cat('\n',ii,'time repeat start')
  t0 <- Sys.time()
  X = data_gen(seed.cluster = 123*ii, seed.data = 654*ii, n=n, true_p=true_p, p=p,
               mu.lower=mu.lower, mu.upper=mu.upper,
               theta = theta, theta_noise=theta_noise, row_group = row_group, col_group = col_group)

  scale.X <- X
  norm.X.f <- norm(scale.X,'f')
  scale.X <- scale.X/norm.X.f
  norm.X.f.vector[ii] <- norm.X.f

  #generate a validation data given the same clustering structure
  val.X <- data_gen(seed.cluster = 123*ii, seed.data = 654*3*ii, n=n, true_p=true_p, p=p,
                    mu.lower=mu.lower, mu.upper=mu.upper,
                    theta = theta, theta_noise=theta_noise, row_group = row_group, col_group = col_group)
  scale.val.X <- val.X
  norm.scale.X.f <- norm(scale.val.X,'f')
  scale.val.X <- scale.val.X/norm.scale.X.f

  #create cluster labels by row and columns
  bi.X.groups = bicluster.label(scale.X)
  bi.val.X.groups = bicluster.label(scale.val.X)

  ######################################################
  #convex biclsutering
  for (para.index in 1:nGamma_1_2){
    cat(para.index)

    gamma_1 = gamma_2 = gammaSeq_1_2[para.index]
    #convex biclustering
    #t1 <- Sys.time()
    res_admm = biADMM.speed(scale.X, nu1, nu2, gamma_1, gamma_2, m=m, phi=phi,niter = 2000,
                            tol = tol,output=0)
    #t2 <- Sys.time()
    bi.MM = cluster_assign(scale.X, m, res_admm,phi = phi) # use scaled version

    bi.admm_group_row = bi.MM$clusters_row
    bi.admm_group_col = bi.MM$clusters_col

    #get the cluster label for bi.admm
    bi.admm.dat = scale.X
    rownames(bi.admm.dat) = bi.admm_group_row$cluster
    colnames(bi.admm.dat) = bi.admm_group_col$cluster

    bi.admm.X.groups = bicluster.label(bi.admm.dat)

    bi.admm.adj.rand = adjustedRandIndex(bi.admm.X.groups$num, bi.X.groups$num)
    bi.admm.validate = predict_bi_validate(scale.X,bi.admm_group_row$cluster,bi.admm_group_col$cluster,
                                           scale.val.X, bi.val.X.groups$num)
    bi.admm.validate.adj.rand = bi.admm.validate$adj.rand

    bi.admm.result.list[[ii]][[para.index]] <-list(A=res_admm$A, iters=res_admm$iters, clustering=bi.MM,adj.rand=bi.admm.adj.rand,validate.adj.rand =bi.admm.validate.adj.rand)
  }
  bi.admm.adj.rand.list[[ii]] = unlist( lapply(bi.admm.result.list[[ii]], function(data) {return(data$adj.rand)}) )
  bi.admm.validate.adj.rand.list[[ii]] = unlist( lapply(bi.admm.result.list[[ii]], function(data) {return(data$validate.adj.rand)}) )

  #the best gamma index and the best adj.rand
  best.bi.admm.adj.rand[ii] = max(bi.admm.adj.rand.list[[ii]] )
  best.bi.admm.gamma.index[ii] = which.max(bi.admm.adj.rand.list[[ii]] )

  best.val.bi.admm.adj.rand[ii] = max(bi.admm.validate.adj.rand.list[[ii]])
  best.val.bi.admm.gamma.index[ii] = which.max(bi.admm.validate.adj.rand.list[[ii]])

  best.bi.admm.gamma.list[[ii]] = data.frame(gamma1=gammaSeq_1_2[best.bi.admm.gamma.index[ii]],
                                             gamma2=gammaSeq_1_2[best.bi.admm.gamma.index[ii]])
  best.bi.admm.validate.gamma.list[[ii]] = data.frame(gamma1=gammaSeq_1_2[best.val.bi.admm.gamma.index[ii]],
                                                      gamma2=gammaSeq_1_2[best.val.bi.admm.gamma.index[ii]])
  #obtain the best convex biclustering
  best.bi.A.list[[ii]] = ( (bi.admm.result.list[[ii]])[[ best.val.bi.admm.gamma.index[ii] ]] )$A
  best.bi.iters[ii] = ( (bi.admm.result.list[[ii]])[[ best.val.bi.admm.gamma.index[ii] ]] )$iters

  best.bi.MM = ( (bi.admm.result.list[[ii]])[[ best.val.bi.admm.gamma.index[ii] ]] )$clustering
  best.bi.admm_group_row = best.bi.MM$clusters_row
  best.bi.admm_group_col = best.bi.MM$clusters_col

  #get the cluster label for bi.admm
  best.bi.admm.dat = scale.X
  rownames(best.bi.admm.dat) = best.bi.admm_group_row$cluster
  colnames(best.bi.admm.dat) = best.bi.admm_group_col$cluster

  best.bi.admm.X.groups.list[[ii]] = bicluster.label(best.bi.admm.dat)

  cat('\n bicluster finished, start Sparse\n')
  ##############################################
  #Sparse convex biclustering
  #print (best.bi.A.list[[ii]])
  feature_weight <- 1/ apply(best.bi.A.list[[ii]], 2, norm,"2")

    for (para.index in 1:nrow(gamma_12_3.grid)){
      cat(para.index)
      gamma_1 = gamma_2 = gamma_12_3.grid[para.index,1]
      gamma_3 = gamma_12_3.grid[para.index,2]
      start.index <-gamma_12_3.grid[para.index,4]
        if (para.index == 1){res1 = sparse.biADMM.speed(scale.X,  nu1, nu2, nu3,
                                                        gamma_1=gamma_1, gamma_2=gamma_2, gamma_3=gamma_3,
                                                        feature_weight = feature_weight,
                                                        m , phi,niter = 2000,tol = tol,output = 0)
        } else {res1 = sparse.biADMM.speed.WS(scale.X, A =scb.admm.result.list[[ii]][[start.index]]$A, nu1, nu2, nu3,v=scb.admm.result.list[[ii]][[start.index]]$v ,z=scb.admm.result.list[[ii]][[start.index]]$z ,g=scb.admm.result.list[[ii]][[start.index]]$g ,
                                              gamma_1=gamma_1, gamma_2=gamma_2, gamma_3=gamma_3,
                                              feature_weight = feature_weight,
                                              m , phi,niter = 2000,tol = tol,output = 0)
        }
      scb.MM = cluster_assign(scale.X, m, res1,phi = phi) # use scaled version and contain res1
      scb.admm_group_row = scb.MM$clusters_row
      scb.admm_group_col = scb.MM$clusters_col
      # get the cluster label for scb.admm
      scb.admm.dat = scale.X
      rownames(scb.admm.dat) = scb.admm_group_row$cluster
      colnames(scb.admm.dat) = scb.admm_group_col$cluster
      scb.admm.X.groups = bicluster.label(scb.admm.dat)
      scb.admm.adj.rand = adjustedRandIndex(scb.admm.X.groups$num, bi.X.groups$num)
      scb.admm.validate = predict_bi_validate(scale.X,scb.admm_group_row$cluster,scb.admm_group_col$cluster,
                                              scale.val.X, bi.val.X.groups$num)
      scb.admm.validate.adj.rand = scb.admm.validate$adj.rand
      scb.admm.result.list[[ii]][[para.index]] <-list(A=res1$A, v=res1$v, z=res1$z, g=res1$g, iters=res1$iters, clustering=scb.MM,adj.rand=scb.admm.adj.rand,validate.adj.rand =scb.admm.validate.adj.rand)
    }
  scb.admm.adj.rand.list[[ii]] = unlist( lapply(scb.admm.result.list[[ii]], function(data) {return(data$adj.rand)}) )
  scb.admm.validate.adj.rand.list[[ii]] = unlist( lapply(scb.admm.result.list[[ii]], function(data) {return(data$validate.adj.rand)}) )
  #the best gamma index and the best adj.rand
  best.scb.admm.adj.rand[ii] = max(scb.admm.adj.rand.list[[ii]] )
  best.scb.admm.gamma.index[ii] = which.max(scb.admm.adj.rand.list[[ii]] )
  best.val.scb.admm.adj.rand[ii] = max(scb.admm.validate.adj.rand.list[[ii]])
  best.val.scb.admm.gamma.index[ii] = which.max(scb.admm.validate.adj.rand.list[[ii]])

  best.scb.admm.gamma.list[[ii]] = data.frame(gamma1=gamma_12_3.grid[best.scb.admm.gamma.index[ii],1],
                                              gamma2=gamma_12_3.grid[best.scb.admm.gamma.index[ii],1],
                                              gamma3=gamma_12_3.grid[best.scb.admm.gamma.index[ii],2])
  best.scb.admm.validate.gamma.list[[ii]] = data.frame(gamma1=gamma_12_3.grid[best.val.scb.admm.gamma.index[ii],1],
                                                       gamma2=gamma_12_3.grid[best.val.scb.admm.gamma.index[ii],1],
                                                       gamma3=gamma_12_3.grid[best.val.scb.admm.gamma.index[ii],2])
  #obtain the best convex scbclustering
  best.scb.A.list[[ii]] = ( (scb.admm.result.list[[ii]])[[ best.val.scb.admm.gamma.index[ii] ]] )$A
  best.scb.iters[ii] = ( (scb.admm.result.list[[ii]])[[ best.val.scb.admm.gamma.index[ii] ]] )$iters
  best.scb.MM = ( (scb.admm.result.list[[ii]])[[ best.val.scb.admm.gamma.index[ii] ]] )$clustering
  best.scb.admm_group_row = best.scb.MM$clusters_row
  best.scb.admm_group_col = best.scb.MM$clusters_col

  #get the cluster label for scb.admm
  best.scb.admm.dat = scale.X
  rownames(best.scb.admm.dat) = best.scb.admm_group_row$cluster
  colnames(best.scb.admm.dat) = best.scb.admm_group_col$cluster

  best.scb.admm.X.groups = bicluster.label(best.scb.admm.dat)

  times[ii] = Sys.time() - t0
  cat( 'total times of ',ii,'repeat is',times[ii])
  #delete the results of ii loop for saving storage
  scb.admm.result.list[[ii]]<-NULL
  #delete v,z,g in this loop for storage burden
  #for (para.index in 1:nrow(gamma_12_3.grid) ){
  #scb.admm.result.list[[ii]][[para.index]]$v <-NULL
  #scb.admm.result.list[[ii]][[para.index]]$z <-NULL
  #scb.admm.result.list[[ii]][[para.index]]$g <-NULL
  #}
}

table(best.scb.admm.gamma.index)
MSD(best.bi.admm.adj.rand)
MSD(best.scb.admm.adj.rand)

MSD(best.val.bi.admm.adj.rand)
MSD(best.val.scb.admm.adj.rand)

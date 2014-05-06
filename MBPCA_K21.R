# MBPCA for K21 data
# written by Chen Yue on 05/01/2014

# load data

load("adj_proc_xc_rest_sfnfwav_craddock_0128.Rdata")

# The third data is not good
preDat <- graphs[-3,,,]

node_pick <- 1:128
pickDat <- preDat[,,node_pick, node_pick]
XX <- NULL
for(i in 1:dim(preDat)[1]){
  for(j in 1:2){
    currMat <- pickDat[i,j,,]
    XX <- rbind(XX, currMat[lower.tri(currMat)])
  }
}
J <- rep(2, times=20)
set.seed(0822)
X <- matrix(rbinom(prod(dim(XX)), 1, c(abs(XX))), dim(XX)[1], dim(XX)[2])
source('MBPCA.R')
fit1 <- MBPCA(X,J,k_between=5, k_within=5, n_iter=30)
save(fit1, file="MBPCA_K21.rda")

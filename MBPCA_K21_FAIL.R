# MBPCA for K21 data
# written by Chen Yue on 05/01/2014

# load data

load("C:/Users/Chen/Dropbox/Research/PCA_working_group/code/Data_Kirby21/adj_proc_xc_rest_sfnfwav_craddock_0128.Rdata")

# The third data is not good
preDat <- graphs[-3,,,]

node_pick <- 1:40
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
source('C:/Users/Chen/Documents/GitHub/Graphical_ICC/MBPCA.R')
fit1 <- MBPCA(X,J,k_between=2, k_within=2, n_iter=30)
plotDat <- data.frame(latent=c(abs(XX)), recon=c(fit1$recon_map))
# dim(svd(fit1$pc_between)$u)
p1 <- ggplot(plotDat, aes(x=latent, y=recon))
p1 <- p1 + layer(geom="point", colour='darkblue', alpha=0.3, size=1)
p1 <- p1 + xlab('Latent Probability') + ylab('Reconstructed Probability') + ggtitle('Simulation Results')
p1

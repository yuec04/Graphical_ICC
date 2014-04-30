# MBPCA Project
# 2/05/2014

# Generate a test data
D <- 8

# define the true mu
mu_true <- rep(0, D)

# define the true theta, theta is dimension D by K1
# theta is the between subject pc
K1 <- 2
theta1_true <- c(rep(1,D/2),rep(-1,D/2))
theta2_true <- c(rep(1,D/4),rep(-1,D/4),rep(1,D/4),rep(-1,D/4))
theta_true <- cbind(theta1_true,theta2_true)

# define the true psi, psi is dimension D by K2
# psi is the within subject pc
K2 <- 2
psi1_true <- rep(.2,D)
psi2_true <- seq(-.5,.5,length.out=D)
psi_true <- cbind(psi1_true, psi2_true)

# generate standard normal scores,
# between subject score V: dimension is I by K1
# within subject score U: dimension is (sum Ji) by K2
I <- 100
set.seed(0822)
J <- floor(runif(I,5,10))
V <- matrix(rnorm(I*K1), I, K1)
U <- matrix(rnorm(sum(J)*K2), sum(J), K2)
V_expand <- t(matrix(unlist(apply(cbind(V,J), 1, function(x){
  matrix(rep(x[1:K1], times=x[K1+1]),K1,x[K1+1])
})),2, sum(J)))

# generate the latent variable dimension is sum(Ji) by D
bet_sub <- t(theta_true %*% t(V_expand))
wit_sub <- t(psi_true %*% t(U))
mu_expand <- matrix(rep(mu_true, each=sum(J)),sum(J),D)
sigma_f <- function(x){
  (1+exp(-x))^(-1)
}
lambda_f <- function(x){
  (0.5-(1+exp(-x))^(-1))/2/x
}
latent <- sigma_f(mu_expand+wit_sub+bet_sub)

# generate data X, dimension sum(Ji) by D
X <- matrix(rbinom(sum(J)*D, 1, as.vector(latent)), sum(J), D)


#############################################
#############################################
# one step algorithm
#############################################
#############################################
# Step 0: Dimension parameter
# function input: X, J, Kb, Kw
I <- length(J)
if(nrow(X)!=sum(J)){
  cat('The dimension does not match!')
}
D <- ncol(X)
Kb <- 2 # number of pc between
Kw <- 2 # number of pc within

# Step 1: Initialize
# initialize mu, theta, psi, xi
mu_old <- rep(.1, D)
theta_old <- matrix(rep(.1, D*Kb),D, Kb)
psi_old <- matrix(rep(.1, D*Kw), D, Kw)
set.seed(0822)
xi <- matrix(rnorm(D*sum(J)), sum(J),D)

# Step 2: Posterior mean and variance step
# test i=1
i <- 1
Ji <- which(unlist(mapply(rep,1:I,J))==i)
xi_i <- xi[Ji,]
X_i <- X[Ji,]

H_i <- diag(1,Kb)-2*matrix(rowSums(apply(rbind(lambda_f(xi_i),t(theta_old)), 2, function(x){
  as.vector(sum(x[1:J[i]])*x[-c(1:J[i])]%o%x[-c(1:J[i])])
} )),Kb,Kb)

G_i <- NULL
B_i <- NULL
g_i <- NULL
for(j in 1:J[i]){
  G_ij <- as.vector(diag(1,Kw))- 2*rowSums(apply(rbind(lambda_f(xi_i[j,]),t(theta_old)), 2, function(x){
    as.vector(x[1]*x[-1]%o%x[-1])
  }))
  G_i <- rbind(G_i, G_ij)
  B_ij <- rowSums(apply(rbind(lambda_f(xi_i[j,]),t(theta_old),t(psi_old)), 2, function(x){
    as.vector(x[1]*x[2:(1+Kb)]%o%x[(2+Kb):(1+Kb+Kw)])
  }))
  B_i <- rbind(B_i, B_ij)
  g_ij <- rowSums(apply(cbind(X_i[j,]-.5+2*lambda_f(xi_i[j,])*mu_old, theta_old), 1, function(x){
    x[1]*x[-1]
  }
                ))
  g_i <- rbind(g_i, g_ij)
}
h_i <- rowSums(apply(cbind(colSums(X_i - .5 + 2 * t(t(lambda_f(xi_i)) * mu_old)), theta_old ), 1, function(x){
  x[1]*x[-1]
}))

# Now we are calculating the inverse of C_i
diag_paste <- function(A, B){
  if(!is.null(A)){
    rbind(cbind(A, matrix(0, dim(A)[1], dim(B)[2])), cbind(matrix(0, dim(B)[1], dim(A)[2]), B))
  } else {
    return(B)
  }
}
D_inv <- NULL
D <- NULL
B <- NULL
m1_i <- h_i
for (j in 1:J[i]){
  D_inv <- diag_paste(D_inv, solve(matrix(G_i[j,],Kw, Kw)))
  D <- diag_paste(D, matrix(G_i[j,],Kw, Kw))
  B <- cbind(B, matrix(B_i[j,], Kb, Kw))  
  m1_i <- c(m1_i, g_i[j,])
}
Cinv_tl <- solve(H_i- B %*% D_inv %*% t(B))
Cinv_tr <- -Cinv_tl%*%B%*%D_inv
Cinv_br <- D_inv + D_inv%*%t(B)%*%Cinv_tl%*%B%*%D
m_i <- rbind(cbind(Cinv_tl, Cinv_tr), cbind(t(Cinv_tr),Cinv_br)) %*% matrix(m1_i, length(m1_i),1)

# Now updating xi_i
VVi <- Cinv_tl + m_i[1:Kb] %o% m_i[1:Kb]
for (j in 1:J[i]){
  j_index <- (Kw*(j-1)+1):(Kw*j)
  UUij <- Cinv_br[j_index,j_index] + m_i[Kb+j_index] %o% m_i[Kb+j_index]
  UVij <- Cinv_tr[, j_index] + m_i[1:Kb] %o% m_i[Kb+j_index]
  xi_i[j,] <- (diag(theta_old %*% VVi %*% t(theta_old)) + diag(psi_old %*% UUij %*% t(psi_old)) + 
    2* diag(theta_old %*% UVij %*% t(psi_old)) + c(2*mu_old*theta_old%*%matrix(m_i[1:Kb],Kb,1)) + 
    c(2* mu_old * psi_old %*% matrix(m_i[Kb+j_index], Kw, 1)) + mu_old^2)^.5
}


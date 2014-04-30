# MBPCA Project
# 2/05/2014

# Generate a test data
D <- 100

# define the true mu
mu_true <- rep(0, D)

# define the true theta, theta is dimension D by K1
# theta is the between subject pc
K1 <- 2
theta1_true <- c(rep(1,D/2),rep(-1,D/2))
theta2_true <- c(rep(.5,D/4),rep(-.5,D/4),rep(.5,D/4),rep(-.5,D/4))
theta_true <- cbind(theta1_true,theta2_true)

# define the true psi, psi is dimension D by K2
# psi is the within subject pc
K2 <- 2
psi1_true <- rep(0.4,D)
psi2_true <- seq(-.3,.3,length.out=D)
psi_true <- cbind(psi1_true, psi2_true)

# generate standard normal scores,
# between subject score V: dimension is I by K1
# within subject score U: dimension is (sum Ji) by K2
I <- 500
set.seed(0822)
J <- rep(2, I)
V <- matrix(rnorm(I*K1), I, K1)
U <- matrix(rnorm(sum(J)*K2), sum(J), K2)
V_expand <- t(matrix(unlist(apply(cbind(V,J), 1, function(x){
  matrix(rep(x[1:K1], times=x[K1+1]),K1,x[K1+1])
})), K1, sum(J)))

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

# test 
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
set.seed(0822)
theta_old <- matrix(rnorm(D*Kb, 0, 1),D, Kb)
psi_old <- matrix(rnorm(D*Kw, 0, 0.1), D, Kw)
xi <- matrix(rnorm(D*sum(J)), sum(J),D)
EM_iter <- 0
diag_paste <- function(A, B){
  if(!is.null(A)){
    rbind(cbind(A, matrix(0, dim(A)[1], dim(B)[2])), cbind(matrix(0, dim(B)[1], dim(A)[2]), B))
  } else {
    return(B)
  }
}

while(EM_iter < 30){
  T1 <- 0
  T2 <- 0
  post <- NULL
  for (e_step_rep in 1:2){
    for (i in 1:I){
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
        g_ij <- rowSums(apply(cbind(X_i[j,]-.5+2*lambda_f(xi_i[j,])*mu_old, psi_old), 1, function(x){
          x[1]*x[-1]
        }
        ))
        g_i <- rbind(g_i, g_ij)
      }
      h_i <- rowSums(apply(cbind(colSums(X_i - .5 + 2 * t(t(lambda_f(xi_i)) * mu_old)), theta_old ), 1, function(x){
        x[1]*x[-1]
      }))
      
      # Now we are calculating the inverse of C_i
      DD_inv <- NULL
      B <- NULL
      m1_i <- h_i
      for (j in 1:J[i]){
        DD_inv <- diag_paste(DD_inv, solve(matrix(G_i[j,],Kw, Kw)))
        B <- cbind(B, matrix(B_i[j,], Kb, Kw))  
        m1_i <- c(m1_i, g_i[j,])
      }
      Cinv_tl <- solve(H_i- B %*% DD_inv %*% t(B))
      Cinv_tr <- -Cinv_tl%*%B%*%DD_inv
      Cinv_br <- DD_inv - DD_inv%*%t(B)%*%Cinv_tr
      m_i <- rbind(cbind(Cinv_tl, Cinv_tr), cbind(t(Cinv_tr),Cinv_br)) %*% matrix(m1_i, length(m1_i),1)
      
      # Step 3: Updating xi_i
      if(e_step_rep!=2){
        VVi <- Cinv_tl + m_i[1:Kb] %o% m_i[1:Kb]
        for (j in 1:J[i]){
          j_index <- (Kw*(j-1)+1):(Kw*j)
          UUij <- Cinv_br[j_index,j_index] + m_i[Kb+j_index] %o% m_i[Kb+j_index]
          UVij <- Cinv_tr[, j_index] + m_i[1:Kb] %o% m_i[Kb+j_index]
          xi_i[j,] <- (diag(theta_old %*% VVi %*% t(theta_old)) + diag(psi_old %*% UUij %*% t(psi_old)) + 
                         2* diag(theta_old %*% UVij %*% t(psi_old)) + c(2*mu_old*(theta_old%*%matrix(m_i[1:Kb],Kb,1))) + 
                         c(2* mu_old * (psi_old %*% matrix(m_i[Kb+j_index], Kw, 1))) + mu_old^2)^.5
        }
      } else {
        VVi <- Cinv_tl + m_i[1:Kb] %o% m_i[1:Kb]
        for (j in 1:J[i]){
          j_index <- (Kw*(j-1)+1):(Kw*j)
          UUij <- Cinv_br[j_index,j_index] + m_i[Kb+j_index] %o% m_i[Kb+j_index]
          UVij <- Cinv_tr[, j_index] + m_i[1:Kb] %o% m_i[Kb+j_index]
          xi_i[j,] <- (diag(theta_old %*% VVi %*% t(theta_old)) + diag(psi_old %*% UUij %*% t(psi_old)) + 
                         2* diag(theta_old %*% UVij %*% t(psi_old)) + c(2*mu_old*(theta_old%*%matrix(m_i[1:Kb],Kb,1))) + 
                         c(2* mu_old * (psi_old %*% matrix(m_i[Kb+j_index], Kw, 1))) + mu_old^2)^.5
          curr_mat <- rbind(cbind(UUij, t(UVij), m_i[Kb+j_index]), cbind(UVij, VVi, m_i[1:Kb]), c(m_i[Kb+j_index], m_i[1:Kb], 1))
          curr_mat_D <- c(curr_mat) %o% lambda_f(xi_i[j,])
          curr_vec <- curr_mat[,Kb+Kw+1]
          curr_vec_D <- curr_vec %o% (X_i[j,]-0.5)
          T1 <- T1 + curr_mat_D
          T2 <- T2 + curr_vec_D
        }
        post <- rbind(post, c(m_i))
      }
      xi[Ji,] <- xi_i
      
    }
  }
  
  # M-step
  
  phi_new <- NULL
  for (d in 1:D){
    phi_est <- c(-solve(2*matrix(T1[,d], Kb+Kw+1, Kb+Kw+1))%*%T2[,d])
    phi_new <- rbind(phi_new, phi_est)
  }
  mu_new <- c(phi_new[,Kb+Kw+1])
  theta_new <- phi_new[,((Kw+1):(Kw+Kb))]
  psi_new <- phi_new[,(1:Kw)]
  
  # get the error
  error <- sum((mu_new-mu_old)^2)+sum((theta_new-theta_old)^2) + sum((psi_new-psi_old)^2)
  cat("EM Iteration No.",EM_iter,'.\n') 
  cat("The current error is ", error, '.\n')
  
  # Update the parameter estimate
  mu_old <- mu_new
  theta_old <- theta_new
  psi_old <- psi_new
  EM_iter <- EM_iter + 1
}

# test the generated latent variable
mu_est_expand <- matrix(rep(mu_new, each=sum(J)), sum(J), D)
V_est_expand <- t(matrix(unlist(apply(cbind(post[,1:Kb],J), 1, function(x){
  matrix(rep(x[1:K1], times=x[K1+1]),K1,x[K1+1])
})), Kb, sum(J)))
U_est_expand <- NULL
for(i in 1:I){
  U_est_expand <- rbind(U_est_expand, t(matrix(post[i,-c(1:Kb)], Kw, J[i])))
}
bet_est_expand <- V_est_expand %*% t(theta_new)
wit_est_expand <- U_est_expand %*% t(psi_new)
recon <- mu_est_expand + bet_est_expand + wit_est_expand

# Simulation result illustration
# Plot the reconstruction probability versus the true probability

library(ggplot2)

qplot(c(latent), c(sigma_f(recon)), xlab='True Probability', 
      ylab='Estimated probability', main='Simulation results', col=2
      )

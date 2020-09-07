## SPMIX Data Simulation ##

# Set working directory
setwd("~/Documents/SPMIX")

# Required libraries
library("mvtnorm")
library("reticulate")

# Python module inclusion
use_python("/bin/python3")
u <- import_from_path("utilities", path = "./inc/")

# Dimensions and fixed parameters
I <- 4; H <- 3; N <- 100
mu0 <- 0
a <- 2; b <- 2
lambda <- 0.1
m_tilde <- rep(0,I*(H-1))
G <- diag(0,I,I); G[1,2] <- 1; G[2,3] <- 1; G[3,4] <- 1; G <- G + t(G)
set.seed(230196)

# Simulating rho
rho <- rbeta(1, shape1 = 1, shape2 = 1)

# Simulating tau kernels
tau <- matrix(0, nrow = H, ncol = 2)
for (h in 1:H) {
  tau[h,2] <- 1/rgamma(1,a,b)
  tau[h,1] <- rnorm(1,mu0,sqrt(tau[h,2]/lambda))
}

# Simulating Sigma
sigma <- 1/rgamma(1,a,b)
Sigma <- sigma*diag(nrow = H-1, ncol = H-1)

# Building F matrix
F_matrix <- matrix(0,I,I)
for (i in 1:I) {
  F_matrix[i,i] <- rho*sum(G[i,]+1-rho)
}

# Simulating weights matrix
w_tilde <- rmvnorm(1, mean=m_tilde, sigma = solve((F_matrix - rho*G) %x% solve(Sigma)))
w_tilde_mat <- matrix(w_tilde, nrow = I, ncol = H-1, byrow = T)
w_mat <- u$inv_alr(w_tilde_mat)

# Simulating data
y <- matrix(0, nrow = N, ncol = I)
for (i in 1:I) {
  for (h in 1:h) {
    y[,i] <- y[,i] + w_mat[i,h]*rnorm(N,tau[h,1],sqrt(tau[h,2]))
  }
}
colnames(y) <- c("Area1","Area2","Area3","Area4")

# Simulating allocation matrix
S <- matrix(0, nrow = N, ncol = I)
for (i in 1:I) {
  S[,i] <- sample(1:H, N, replace = T, prob = w_mat[i,])
}
colnames(S) <- colnames(y)

# Save data image for importation
rm(list='u')
save.image(file = "SPMIX_simulated_data.Rdata")

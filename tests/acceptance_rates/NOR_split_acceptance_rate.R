## First test on acceptance rates computation according to the NORETS approach ##

# Set working directory
setwd("~/Documents/SPMIX/")

# Required Libraries
library("mvtnorm")
library("reticulate")
library("madness")

# Python module inclusion
use_python("/bin/python3")
u <- import_from_path("utilities", path = "./inc")

# Data importation
load("./data/SPMIX_simulated_data.Rdata")

# Dimension growth with Nor fashion ---------------------------------------

###########################################################################
# Determination of proposals parameters (LONG EXECUTION) ------------------
input_list <- list(
  "data" = y,
  "w_tilde" = w_tilde,
  "tau" = tau,
  "h_dim" = as.integer(H+1)
)

param_list <- list(
  "a" = a,
  "b" = b,
  "mu0" = mu0,
  "lambda" = lambda,
  "m_tilde" = matrix(rep(0,I*H),nrow=1),
  "S_tilde" = solve((F_matrix - rho*G) %x% solve(sigma*diag(nrow = H, ncol = H)))
)

opt_options <- list(
  "max_iter" = 10L,
  "tol" = 1e-8,
  "w_0" = matrix(c(0,0,0,0), nrow = 1),
  "tau_0" = matrix(c(0,1), nrow = 1)
)

opt_options_cont <- list(
  "max_iter" = 10L,
  "tol" = 1e-8,
  "w_0" = matrix(tmp$post_mode[1,1:4], nrow=1),
  "tau_0" = matrix(tmp$post_mode[1,5:6], nrow=1)
)

toc = Sys.time()
tmp <- u$norets_optimal_proposal(input_list, param_list, opt_options_cont)
tic = Sys.time()
print(paste0("Iter Time: ", tic-toc))
###########################################################################

# To speed-up, the data has already been stored
load("./data/opt_proposal_output.Rdata")

# Sampling from the approximated posterior
# Symmetry correction
tmp$post_variance[lower.tri(tmp$post_variance)] <- 0; 
diag(tmp$post_variance) <- 0.5*diag(tmp$post_variance)
post_variance <- solve(tmp$post_variance + t(tmp$post_variance))

# Sampling (block-independent, even if the hessian is full)
set.seed(230196)
w_tilde_toadd <- t(rmvnorm(1, tmp$post_mode[1:I], post_variance[1:I,1:I]))
tau_toadd <- c(rnorm(1,tmp$post_mode[I+1], post_variance[I+1,I+1]),
               rgamma(1,tmp$post_mode[I+2], post_variance[I+2,I+2]))

# Building new entities
H_new <- as.integer(H+1)
w_tilde_mat_new <- cbind(w_tilde_mat, w_tilde_toadd)
w_tilde_new <- rbind(c(t(w_tilde_mat_new)))
tau_new <- rbind(tau,tau_toadd); row.names(tau_new) <- NULL
m_tilde_new <- rep(0,dim(w_tilde_new)[2])
Sigma_new <- sigma*diag(nrow = H_new-1, ncol = H_new-1)

# Computing acceptance rate
param_num = list(
  "a" = a,
  "b" = b,
  "mu0" = mu0,
  "lambda" = lambda,
  "m_tilde" = m_tilde_new,
  "S_tilde" = solve((F_matrix - rho*G) %x% solve(Sigma_new))
)

param_den = list(
  "a" = a,
  "b" = b,
  "mu0" = mu0,
  "lambda" = lambda,
  "m_tilde" = m_tilde,
  "S_tilde" = solve((F_matrix - rho*G) %x% solve(Sigma))
)

# Computing loglikelihood (to stabilize the final value)
L_num <- u$spmix_loglikelihood(data = y, w_tilde = w_tilde_new, tau = tau_new,
                               h_dim = as.integer(H_new), param = param_num)
L_den <- u$spmix_loglikelihood(data = y, w_tilde = w_tilde, tau = tau,
                               h_dim = as.integer(H), param = param_den)

# Adding the term coming from approximated posterior
L_den <- L_den + dmvnorm(t(w_tilde_toadd), tmp$post_mode[1:I], post_variance[1:I,1:I], log = T) +
  dnorm(tau_toadd[1], tmp$post_mode[I+1], post_variance[I+1,I+1], log = T) + 
  dgamma(tau_toadd[2], tmp$post_mode[I+2], post_variance[I+2,I+2], log = T)
# Computing alpha
alpha <- exp(L_num - L_den)

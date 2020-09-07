## First test on acceptance rates computation according to the RG approach ##

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

# Split Move with R&J Fashion ---------------------------------------------

# Choose the component to split
set.seed(230196)
to_split <- 1L # in general is simulated as: sample(1:H, 1, prob = rep(1/H,H))

# Sample the u vector
u1 <- rbeta(I,2,2)
u2 <- rbeta(1,2,2)
u3 <- rbeta(1,1,1)

# Compute new weights carrying the derivatives through madness
x <- rbind(c(w_mat[,to_split], tau[to_split,],u1,u2,u3)); x <- madness(x)
x_new <- cbind(x[,1:I]*x[,(I+3):(I+6)],
               x[,1:I]*(1-x[,(I+3):(I+6)]),
               x[,I+1]+x[,2*I+3]*sqrt(x[,I+2]*sum(x[,1:I]*(1-x[,(I+3):(I+6)]))/sum(x[,1:I]*x[,(I+3):(I+6)])),
               x[,I+1]-x[,2*I+3]*sqrt(x[,I+2]*sum(x[,1:I]*x[,(I+3):(I+6)])/sum(x[,1:I]*(1-x[,(I+3):(I+6)]))),
               x[,2*I+4]*(1-x[,2*I+3]^2)*x[,I+2]*(sum(x[,1:I])/sum(x[,1:I]*x[,(I+3):(I+6)])),
               (1-x[,2*I+4])*(1-x[,2*I+3]^2)*x[,I+2]*(sum(x[,1:I])/sum(x[,1:I]*(1-x[,(I+3):(I+6)]))) )

# Extracting Jacobian and re-ordering things
Jac <- abs(det(dvdx(x_new)))
w_mat_split <- matrix(val(x_new)[1:(2*I)], nrow = I, ncol = 2, byrow = F)
tau_split <- matrix(val(x_new)[(2*I+1):(2*I+4)], nrow = 2, ncol = 2, byrow = F)

# Creating new parameters to compute the likelihood ratio
H_new <- as.integer(H+1)
w_mat_new <- cbind(w_mat_split, w_mat[,-to_split])
w_tilde_mat_new <- u$alr(w_mat_new)
w_tilde_new <- rbind(c(t(w_tilde_mat_new)))
tau_new <- rbind(tau_split, tau[-to_split,])
m_tilde_new <- rep(0,dim(w_tilde_new)[2])
Sigma_new <- sigma*diag(nrow = H_new-1, ncol = H_new-1)

# TODO: scrivere la mappa per la mossa di merge (serve?)

# Computation of the acceptance rate (R&G approach)
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

# Adding contribution from proposal and finalizing alpha
L_num <- L_num + log(0.5)
L_den <- L_den + log(0.5)
L <- exp(L_num-L_den)
alpha <- (L*Jac)/(prod(dbeta(u1,2,2))*dbeta(u2,2,2)*dbeta(u3,1,1))

# Since data are simulated from a mixture of three weights, it makes sense that
# in the end alpha is practically zero even though I expect that this alpha is
# in general always quite low.

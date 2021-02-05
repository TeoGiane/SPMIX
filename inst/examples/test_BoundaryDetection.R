## Test for boundary detection for the Reversible Jump Sampler ##
# The scenario is described in Sec. 6.1 of Beraha et al. (2020)

# Required libraries
library("SPMIX")
library("ggplot2")

###########################################################################
# Data Generation ---------------------------------------------------------

# Generating data
I <- 6; Ns <- rep(1000,I); set.seed(230196)
data <- list(); labels <- sprintf("g%d", 1:I)
data[[1]] <- metRology::rt.scaled(n=Ns[1], df=6, mean=-4, sd=1)
data[[2]] <- metRology::rt.scaled(n=Ns[2], df=6, mean=-4, sd=1)
data[[3]] <- sn::rsn(n=Ns[3], xi=4, omega=4, alpha=1); attributes(data[[3]]) <- NULL
data[[4]] <- sn::rsn(n=Ns[4], xi=4, omega=4, alpha=1); attributes(data[[4]]) <- NULL
data[[5]] <- rchisq(n=Ns[5], df=3)
data[[6]] <- rchisq(n=Ns[6], df=3)
names(data) <- labels

# Setting initial W
W <- matrix(1,I,I)

###########################################################################

###########################################################################
# Sampler execution -------------------------------------------------------

# Setting MCMC parameters
burnin = 5000
niter = 5000
thin = 2

# Grab input filenames
params_filename = system.file("input_files/rjsampler_params.asciipb", package = "SPMIX")

# Run Spatial sampler
out <- SPMIXSampler(burnin, niter, thin, data, W, params_filename, type = "rjmcmc")

###########################################################################

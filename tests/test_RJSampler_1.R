## First test for the Reversible Jump Sampler ##

# Required libraries
library("SPMIX")

# Generate data (1 location, mixture of 3 normals)
set.seed(230196)
ngroups <- 1; ncomponents <- 3; N <- 1000
means <- c(-4,1,5); std_devs <- c(1,1,1); weights <- c(1/6,3/6,2/6)
cluster_alloc <- sample(1:ncomponents, prob = weights, size = N, replace = T)
data <- list(); data[[1]] <- rnorm(N, mean = means[cluster_alloc], sd = std_devs[cluster_alloc])

# Generate Matrix W
W <- matrix(0, nrow = 1, ncol = 1, byrow = T)

# Run Sampler
# Setting MCMC parameters
burnin = 5000
niter = 5000
thin = 2

# Grab input filenames
params_filename = system.file("input_files/rjsampler_params.asciipb", package = "SPMIX")

# Run Spatial sampler
out <- SPMIXSampler(burnin, niter, thin, data, W, params_filename, type = "rjmcmc")

# Analyses
chains <- sapply(out, function(x) DeserializeSPMIXProto("UnivariateState",x))
H_chain <- sapply(chains, function(x) x$num_components)
means_chain <- lapply(chains, function(x) sapply(x$atoms, function(x) x$mean))
stdev_chain <- lapply(chains, function(x) sapply(x$atoms, function(x) x$stdev))
weights_chain <- lapply(chains, function(x) x$groupParams[[1]]$weights)

# Barplot for the number of components
x11(height = 4, width = 8.27); barplot(table(H_chain)/length(chains))
title("Posterior of H")

# Plotting the average density over iterations and compare with true curve
# Computing estimated density
data_range <- sapply(data, range)
est_dens <- ComputeDensities(chains, 500, data_range)

# Computing true density
x_grid <- seq(data_range[1,1], data_range[2,1], length.out = 500)
xgrid_expand <- t(rbind(replicate(ncomponents, x_grid, simplify = "matrix")))
true_dens <- t(as.matrix(weights)) %*% dnorm(xgrid_expand,means,std_devs)

# Visualization (first sketch)
x11(height = 4, width = 8.27)
plot(x_grid, true_dens, type = 'l', col='blue', lwd=2, lty=1,xlab = "Grid", ylab = "Density")
lines(x_grid, est_dens[[1]], col='darkorange', lwd=2, lty=1)
legend("topleft", legend = c("True", "Estimated"), col=c("blue","darkorange"), lty = 1, lwd = 2)
title("Area 1 - Density estimation")

# Required libraries
library("SPMIX")

# Import proto files
readProtoFiles(system.file("proto/univariate_mixture_state.proto", package = "SPMIX"))

# Generate data (1 location, mixture of 3 normals)
set.seed(230196)
ngroups <- 1; ncomponents <- 3; N <- 1000
means <- c(-4,0,4); std_devs <- c(1,1,1); weights <- c(rep(1/3,3))
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
out <- SPMIX_sampler(burnin, niter, thin, data, W, params_filename, type = "rjmcmc", display_progress = TRUE)
save(out, file = "RJTest1_output_10k.dat")

# Analyses
chains <- sapply(out, function(x) unserialize_proto("UnivariateState",x))
H_chain <- sapply(chains, function(x) x$num_components)
means_chain <- sapply(chains, function(x) sapply(x$atoms, function(x) x$mean))
stdev_chain <- sapply(chains, function(x) sapply(x$atoms, function(x) x$stdev))
weights_chain <- sapply(chains, function(x) x$groupParams[[1]]$weights)

# Barplot for the number of components
x11(height = 4, width = 8.27); barplot(table(H_chain))
title("Posterior of H")

# Plotting the average density over iterations and compare with true curve
# Computing average density
x_grid <- seq(-5,5,length.out = 500)
est_dens <- rep(0,length(x_grid))
for (i in 1:length(chains)) {
  xgrid_expand <- t(rbind(replicate(H_chain[i], x_grid, simplify = "matrix")))
  est_dens <- est_dens + t(as.matrix(weights_chain[[i]])) %*% dnorm(xgrid_expand,
                                                                    means_chain[[i]],stdev_chain[[i]])
}
est_dens <- est_dens / length(chains)

# Computing true density
xgrid_expand <- t(rbind(replicate(ncomponents, x_grid, simplify = "matrix")))
true_dens <- t(as.matrix(weights)) %*% dnorm(xgrid_expand,means,std_devs)

# Visualization (first sketch)
x11(height = 4, width = 8.27)
plot(x_grid, est_dens, type = 'l', col='darkorange', lwd=2, lty=1,
     xlab = "Grid", ylab = "Density")
lines(x_grid, true_dens, col='blue', lwd=2, lty=1)
legend("bottomright", legend = c("Estimated", "True"), col=c("darkorange","blue"), lty = 1, lwd = 2, )
title("Area 1 - Density estimation")

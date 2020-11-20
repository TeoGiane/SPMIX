## First test for the RJSampler ##
# The scenario is described in Sec. 6.2 of Beraha et al. (2020)

# Required libraries
library("SPMIX")
library("ggplot2")
library("sp") # also rgdal, rgeos, maptools have been installed

###########################################################################
# Data Generation ---------------------------------------------------------

# Generate spatial grid from coordinates
set.seed(230196)
numGroups <- 9; numComponents <- 3
coords <- expand.grid(x = seq(0,1,length.out = sqrt(numGroups)),
                      y = seq(0,1,length.out = sqrt(numGroups)))
xbar <- ybar <- 0.5
spatial_grid <- as.SpatialPolygons.GridTopology(SpatialPixels(SpatialPoints(coords))@grid)
labels <- row.names(coordinates(spatial_grid))
rm(list='coords')

# Generating weights matrices
trasformed_weights <- matrix(0,nrow=numGroups,ncol=numComponents)
weights <- matrix(0,nrow=numGroups,ncol=numComponents)
for (i in 1:numGroups) {
  center <- coordinates(spatial_grid)[i,]
  trasformed_weights[i,1] <- 3*(center[1]-xbar)+3*(center[2]-ybar)
  trasformed_weights[i,2] <- -3*(center[1]-xbar)-3*(center[2]-ybar)
  weights[i,] <- inv_alr(trasformed_weights[i,], padded_zero = TRUE)
}
rm(list=c('i','center'))
row.names(trasformed_weights) <- labels
row.names(weights) <- labels

# Generate data
means <- c(-5,0,5)
sds <- c(1,1,1)
Ns <- rep(50, numGroups)
data <- list()
for (i in 1:numGroups) {
  cluster_alloc <- sample(1:numComponents, prob = weights[i,], size = Ns[i], replace = T)
  data[[i]] <- rnorm(Ns[i], mean = means[cluster_alloc], sd = sds[cluster_alloc])
}
names(data) <- labels
rm(list=c('cluster_alloc','i'))

# Build adjacency matrix (FIND A BETTER WAY)
W <- matrix(0, numGroups, numGroups); colnames(W) <- rownames(W) <- labels
check_prox <- function(x,y) {
  dist <- 1/(sqrt(numGroups)-1)
  return(ifelse((x[1]==y[1] & abs(x[2]-y[2])<=dist) | (x[2]==y[2] & abs(x[1]-y[1])<=dist), 1, 0))
}
for (i in 1:numGroups) {
  for (j in 1:numGroups) {
    x_coor <- coordinates(spatial_grid)[labels[i],]
    y_coor <- coordinates(spatial_grid)[labels[j],]
    W[i,j] <- check_prox(x_coor, y_coor)
  }
  W[i,i] <- 0
}
rm(list=c('i','j','x_coor','y_coor','check_prox'))

###########################################################################

###########################################################################
# Sampler Execution -------------------------------------------------------

# Setting MCMC parameters
burnin = 5000
niter = 5000
thin = 2

# Grab input filenames
params_filename = system.file("input_files/rjsampler_params.asciipb", package = "SPMIX")
# options_filename = system.file("input_files/optimization_options.asciipb", package = "SPMIX")

# Run Spatial sampler
out <- SPMIX_sampler(burnin, niter, thin, data, W, params_filename, type = "rjmcmc")
save('out', file = "SPMIXoutput10K.dat")

# Data analysis on chains (ONGOING)
chains <- sapply(out, function(x) unserialize_proto("UnivariateState",x))
H_chain <- sapply(chains, function(x) x$num_components)
means_chain <- sapply(chains, function(x) sapply(x$atoms, function(x) x$mean))
stdev_chain <- sapply(chains, function(x) sapply(x$atoms, function(x) x$stdev))
weights_chain <- sapply(chains, function(x) x$groupParams[[1]]$weights)
x_grid <- seq(-6,6, length.out = 500)
est_dens <- seq(0, length.out = length(x_grid))
for (i in 1:length(chains)) {
  for (h in 1:H_chain[i]) {
    est_dens <- est_dens + weights_chain[[i]][h]*dnorm(x_grid,means_chain[[i]][h], stdev_chain[[i]][h])
  }
}
est_dens <- est_dens/length(chains)

###########################################################################

###########################################################################
# Visualization -----------------------------------------------------------

# Generate dataframe for plotting
# df <- fortify(spatial_grid)
# for (i in 1:length(labels)) {
#   df[which(df$id==labels[i]),'w_i1'] <- weights[i,1]
#   df[which(df$id==labels[i]),'w_i2'] <- weights[i,2]
# }
# plot1 <- ggplot(data = df, aes(x=long, y=lat, group=group, fill=w_i1)) +
#   geom_polygon() + scale_fill_continuous(type = "gradient")
# plot2 <- ggplot(data = df, aes(x=long, y=lat, group=group, fill=w_i2)) +
#   geom_polygon() + scale_fill_continuous(type = "gradient")
# x11(height = 4, width = 8.27); gridExtra::grid.arrange(plot1, plot2, ncol=2)

###########################################################################

## Second test for the RJSampler ##
# The scenario is described in Sec. 6.2 of Beraha et al. (2020)

# Required libraries
library("SPMIX")
library("ggplot2")
library("sp") # also rgdal, rgeos, maptools have been installed

# Helper functions (INSERT INTO SPMIX PACKAGE??)
buildLattice <- function(coords) {
  numGroups <- dim(coords)[1]

  # Build spatial grid from coordinates
  spPixels <- sp::SpatialPixels(sp::SpatialPoints(coords))
  spatialGrid <- sp::as.SpatialPolygons.GridTopology(spPixels@grid)

  # Compute proximity matrix
  nblist <- sp::gridIndex2nb(spPixels, maxdist = 1)
  W <- matrix(0,numGroups,numGroups)
  for (i in 1:numGroups) {
    row_ind <- i; col_inds <- unlist(nblist[[i]]); names(col_inds) <- NULL
    W[row_ind,col_inds] <- 1; W[i,i] <- 0
  }

  # Return output list
  out <- list(); out[[1]] <- spatialGrid; out[[2]] <- W
  names(out) <- c("spatialGrid", "proxMatrix")
  return(out)
}
posteriorDensities <- function(unserialized_chains, ranges, names = NULL) {

  # Elicit numGroups
  numGroups <- length(unserialized_chains[[1]]$groupParams)

  # Extract necessary chains
  H_chain <- sapply(unserialized_chains, function(x) x$num_components)
  means_chain <- lapply(unserialized_chains, function(x) sapply(x$atoms, function(x) x$mean))
  stdev_chain <- lapply(unserialized_chains, function(x) sapply(x$atoms, function(x) x$stdev))

  # Compute Estimated densities
  estimated_densities <- list()
  for (i in 1:numGroups) {
    weights_chain <- lapply(unserialized_chains, function(x) x$groupParams[[i]]$weights)
    x_grid <- seq(ranges[1,i], ranges[2,i], length.out = 500)
    est_dens <- rep(0,length(x_grid))
    for (j in 1:length(unserialized_chains)) {
      xgrid_expand <- t(rbind(replicate(H_chain[j], x_grid, simplify = "matrix")))
      est_dens <- est_dens + t(as.matrix(weights_chain[[j]])) %*% dnorm(xgrid_expand,
                                                                        means_chain[[j]],stdev_chain[[j]])
    }
    est_dens <- est_dens/length(chains)
    estimated_densities[[i]] <- est_dens
  }

  # Set names if passed
  if (!is.null(names)) {
    names(estimated_densities) <- labels
  }

  # Return list
  return(estimated_densities)
}
dataRecast <- function(x, y, labels = row.names(y)) {
  rows <- dim(y)[1]; cols <- dim(y)[2]
  out <- data.frame()
  for (i in 1:rows) {
    tmp <- data.frame(x, y[i,], as.factor(rep(labels[i],cols)));
    names(tmp) <- c("Grid", "Value", "Density"); row.names(tmp) <- NULL
    out <- rbind(out,tmp)
  }
  return(out)
}

###########################################################################
# Data Generation ---------------------------------------------------------

# Build spatial grid and proximity matrix from coordinates
numGroups <- 9; numComponents <- 3
coords <- expand.grid(x = seq(0,1,length.out = sqrt(numGroups)),
                      y = seq(0,1,length.out = sqrt(numGroups)))
xbar <- ybar <- 0.5

spatial_grid <- buildLattice(coords)$spatialGrid
W <- buildLattice(coords)$proxMatrix
labels <- row.names(coordinates(spatial_grid))
rm(list='coords')

# Generating weights matrices
trasformed_weights <- matrix(0,nrow=numGroups,ncol=numComponents-1)
weights <- matrix(0,nrow=numGroups,ncol=numComponents)
for (i in 1:numGroups) {
  center <- coordinates(spatial_grid)[i,]
  trasformed_weights[i,1] <- 3*(center[1]-xbar)+3*(center[2]-ybar)
  trasformed_weights[i,2] <- -3*(center[1]-xbar)-3*(center[2]-ybar)
  weights[i,] <- InvAlr(trasformed_weights[i,])
}
rm(list=c('i','center'))
row.names(trasformed_weights) <- labels
row.names(weights) <- labels

# Generate data
set.seed(230196)
means <- c(-5,0,5); sds <- c(1,1,1); Ns <- rep(50, numGroups)
data <- list()
for (i in 1:numGroups) {
  cluster_alloc <- sample(1:numComponents, prob = weights[i,], size = Ns[i], replace = T)
  data[[i]] <- rnorm(Ns[i], mean = means[cluster_alloc], sd = sds[cluster_alloc])
}
names(data) <- labels
rm(list=c('cluster_alloc','i'))

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
out <- SPMIXSampler(burnin, niter, thin, data, W, params_filename, type = "rjmcmc")

# Data analysis on chains (ONGOING)
chains <- sapply(out, function(x) UnserializeSPMIXProto("UnivariateState",x))
H_chain <- sapply(chains, function(x) x$num_components)
means_chain <- lapply(chains, function(x) sapply(x$atoms, function(x) x$mean))
stdev_chain <- lapply(chains, function(x) sapply(x$atoms, function(x) x$stdev))

# Computing estimated densities
data_ranges <- sapply(data, range)
estimated_densities <- posteriorDensities(chains, ranges = data_ranges, names = labels)

# Computing true densities for comparison
true_densities <- list()
for (i in 1:numGroups) {
  data_range <- range(data[[i]])
  x_grid <- seq(data_range[1], data_range[2], length.out = 500)
  xgrid_expand <- t(rbind(replicate(numComponents, x_grid, simplify = "matrix")))
  true_dens <- t(as.matrix(weights[i,])) %*% dnorm(xgrid_expand, means, sds)
  true_densities[[i]] <- true_dens
}
names(true_densities) <- labels

###########################################################################

###########################################################################
# Visualization -----------------------------------------------------------

# Weights distribution on the spatial grid - Plot
df <- fortify(spatial_grid)
for (i in 1:length(labels)) {
  df[which(df$id==labels[i]),'w_i1'] <- weights[i,1]
  df[which(df$id==labels[i]),'w_i2'] <- weights[i,2]
}
plot_weightsi1 <- ggplot(data = df, aes(x=long, y=lat, group=group, fill=w_i1)) +
  geom_polygon() + scale_fill_continuous(type = "gradient")
plot_weightsi2 <- ggplot(data = df, aes(x=long, y=lat, group=group, fill=w_i2)) +
  geom_polygon() + scale_fill_continuous(type = "gradient")
rm(list='df')

# Posterior of H - Barplot
df <- as.data.frame(table(H_chain)/length(H_chain)); names(df) <- c("NumComponents", "Prob.")
plot_postH <- ggplot(data = df, aes(x=NumComponents, y=Prob.)) +
  geom_bar(stat="identity", color="steelblue", fill="lightblue") +#"white") +
  theme(plot.title = element_text(face="bold", hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  ggtitle("Posterior of H")
rm(list='df')

# Posterior of H - Traceplot
df <- data.frame("Iteration"=1:niter, "LowPoints"=H_chain-0.3, "UpPoints"=H_chain+0.3)
plot_traceH <- ggplot(data=df, aes(x=Iteration, y=LowPoints, xend=Iteration, yend=UpPoints)) +
  ylim(range(df[,-1])) + ylab("NumComponents") + geom_segment(lwd=0.1) +
  theme(plot.title = element_text(face="bold", hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  ggtitle("Traceplot of H")
rm(list='df')

# Comparison plots between estimated and true densities in i-th area
plots_area <- list()
for (i in 1:numGroups) {
  x <- seq(data_ranges[1,i], data_ranges[2,i], length.out = 500)
  y <- rbind(estimated_densities[[i]], true_densities[[i]]); row.names(y) <- c("Estimated", "True")
  tmp <- ggplot(data = dataRecast(x,y), aes(x=Grid, y=Value, group=Density, col=Density)) +
    geom_line(lwd=1) + theme(plot.title = element_text(face="bold", hjust = 0.5)) +
    ggtitle(paste0("Area ", i))
  plots_area[[i]] <- tmp; rm(list=c('x','y','tmp'))
}
names(plots_area) <- labels

# Printing plots
x11(height = 2, width = 8.27); gridExtra::grid.arrange(plot_weightsi1, plot_weightsi2, plot_postH, ncol=3)
x11(height = 8.27, width = 8.27); gridExtra::grid.arrange(grobs=plots_area, nrow=3, ncol=3)
x11(height = 4, width = 8.27); plot_traceH

###########################################################################

## First test for the RJSampler ##
# The scenario is described in Sec. 6.2 of Beraha et al. (2020)

# Required libraries
library("SPMIX")
library("ggplot2")
library("sp") # also rgdal, rgeos, maptools have been installed

# Getting useful descriptors
readProtoFiles(file = system.file("proto/sampler_params.proto", package = "SPMIX"))
readProtoFiles(files = system.file("proto/optimization_options.proto", package = "SPMIX"))

###########################################################################
# Data Generation ---------------------------------------------------------

# Generate spatial grid from coordinates
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

###########################################################################

###########################################################################
# Sampler Execution -------------------------------------------------------

# Setting MCMC parameters
burnin = 10000
niter = 10000
thin = 5

# Grab input filenames
data_filename = system.file("input_files/datasets/scenario0/rep0.csv", package = "SPMIX")
w_filename = system.file("input_files/prox_matrix.csv", package = "SPMIX")
params_filename = system.file("input_files/rjsampler_params.asciipb", package = "SPMIX")
options_filename = system.file("input_files/optimization_options.asciipb", package = "SPMIX")

# building input objects
data_obj <- readDataFromCSV(data_filename)
w_obj <- readMatrixFromCSV(w_filename)
params_obj <- readASCII(SamplerParams, file(params_filename))
options_obj <- readASCII(OptimOptions, file(options_filename))

# Run test utility
RJsampler_test(data_obj, w_obj, params_obj, options_obj)

###########################################################################

###########################################################################
# Visualization -----------------------------------------------------------

# Generate dataframe for plotting
df <- fortify(spatial_grid)
for (i in 1:length(labels)) {
  df[which(df$id==labels[i]),'w_i1'] <- weights[i,1]
  df[which(df$id==labels[i]),'w_i2'] <- weights[i,2]
}
plot1 <- ggplot(data = df, aes(x=long, y=lat, group=group, fill=w_i1)) +
  geom_polygon() + scale_fill_continuous(type = "gradient")
plot2 <- ggplot(data = df, aes(x=long, y=lat, group=group, fill=w_i2)) +
  geom_polygon() + scale_fill_continuous(type = "gradient")
x11(height = 4, width = 8.27); gridExtra::grid.arrange(plot1, plot2, ncol=2)

###########################################################################

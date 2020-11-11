# Required libraries
library("SPMIX")

# Import proto files
readProtoFiles(system.file("proto/univariate_mixture_state.proto", package = "SPMIX"))

# Generate data (1 location, mixture of 3 normals)
set.seed(230196)
means <- c(-5,1,4)
std_devs <- c(1,1,1)
weights <- c(0.2,0.5,0.3)
transformed_weights <- alr(weights)
data <- list(); data[[1]] <- rep(0,1000)
for (i in 1:length(means)) {
  data[[1]] <- data[[1]] + weights[i]*rnorm(1000, means[i], std_devs[i])
}

# Generate Matrix W
W <- matrix(0,nrow = 1, ncol = 1, byrow = T)

# Read Params
params_filename <- system.file("input_files/sampler_params.asciipb", package = "SPMIX")
params <- readASCII(SamplerParams, file(params_filename))

# Read options
options_filename <- system.file("input_files/optimization_options.asciipb", package = "SPMIX")
options <- readASCII(OptimOptions, file(options_filename))

###########################################################################
# Case 1 - A reduced state, Increase move ---------------------------------
# Generate the state
state <- new(UnivariateState)
state$num_components <- 2
state$atoms <- list(new(UnivariateMixtureAtom, mean = means[1], stdev = std_devs[1]),
                    new(UnivariateMixtureAtom, mean = means[2], stdev = std_devs[2]))
state$groupParams <- new(UnivariateState.GroupParams, weights = inv_alr(alr(weights[1:2])))
state$rho <- 0.9
state$Sigma <- new(EigenMatrix, rows=1, cols=1, data=1.2)

# Run test
cat("IncreaseMove_test ...\n\n")
IncreaseMove_test(data, W, params, state, options)

###########################################################################

###########################################################################
# Case 2 - The right state, Increase move ---------------------------------
# Generate the state
state <- new(UnivariateState)
state$num_components <- 3
state$atoms <- list(new(UnivariateMixtureAtom, mean = means[1], stdev = std_devs[1]),
                    new(UnivariateMixtureAtom, mean = means[2], stdev = std_devs[2]),
                    new(UnivariateMixtureAtom, mean = means[3], stdev = std_devs[3]))
state$groupParams <- new(UnivariateState.GroupParams, weights = weights)
state$rho <- 0.9
state$Sigma <- new(EigenMatrix, rows=2, cols=2, data=c(1.2,0,0,1.2))

# Run test
cat("IncreaseMove_test ...\n\n")
IncreaseMove_test(data, W, params, state, options)
cat("ReduceMove_test ...\n\n")
ReduceMove_test(data, W, params, state, options)

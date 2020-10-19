# Test Environment for teh RJSampler
library("SPMIX")

# Getting useful descriptors
readProtoFiles(file = system.file("proto/sampler_params.proto", package = "SPMIX"))

# Setting MCMC parameters
burnin = 10000
niter = 10000
thin = 5

# Grab input filenames
data_filename = system.file("input_files/datasets/scenario0/rep0.csv", package = "SPMIX")
w_filename = system.file("input_files/prox_matrix.csv", package = "SPMIX")
params_filename = system.file("input_files/rjsampler_params.asciipb", package = "SPMIX")

# building input objects
data_obj <- readDataFromCSV(data_filename)
w_obj <- readMatrixFromCSV(w_filename)
params_obj <- readASCII(SamplerParams, file(params_filename))

# Run test utility
RJsampler_test(data_obj, w_obj, params_obj)

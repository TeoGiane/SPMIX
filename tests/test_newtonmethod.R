library("SPMIX")

# # New Message::NewtonOptions from scratch
# readProtoFiles(system.file("/proto/newton_options.proto", package = "SPMIX"))
# options <- new(NewtonOptions)
# options$max_iter <- 100
# options$tol <- 1e-8;
# write(toString(options),
#       file = "/home/m_gianella/Documents/R-Packages/SPMIX_stuff/InputFiles/newton_options.asciipb")

data_filename <- system.file("input_files/datasets/scenario0/rep0.csv", package = "SPMIX")
w_filename = system.file("input_files/prox_matrix.csv", package = "SPMIX")
params_filename <- system.file("input_files/sampler_params.asciipb", package = "SPMIX")
options_filename <- system.file("input_files/optimization_options.asciipb", package = "SPMIX")

data <- readDataFromCSV(data_filename)
W <- readMatrixFromCSV(w_filename)
params <- RProtoBuf::readASCII(SamplerParams, file(params_filename))
options <- RProtoBuf::readASCII(OptimOptions, file(options_filename))

load(system.file("input_files/output_chains_serialized.dat", package = "SPMIX"))
state <- out[[1]]; rm(list='out'); state <- unserialize_proto("UnivariateState", state)

newton_opt_test(state, data, W, params, options)

library("SPMIX")

data_filename <- system.file("input_files/datasets/scenario0/rep0.csv", package = "SPMIX")
w_filename = system.file("input_files/prox_matrix.csv", package = "SPMIX")
params_filename <- system.file("input_files/sampler_params.asciipb", package = "SPMIX")

data <- readDataFromCSV(data_filename)
W <- readMatrixFromCSV(w_filename)
params <- RProtoBuf::readASCII(SamplerParams, file(params_filename))

load(system.file("input_files/output_chains_serialized.dat", package = "SPMIX"))
state <- out[[1]]; rm(list='out'); state <- unserialize_proto("UnivariateState", state)

params_str <- toString(params); state_str <- toString(state);
test_spmixloglikelihood(data, W, params_str, state_str)

# First attempt of simulation
library("SPMIX")

# Setting MCMC parameters
burnin = 10000
niter = 10000
thin = 5

# Setting input filenames (ALSO FAKE TO CHECK EXCETIONS HANDLING)
data_filename = system.file("input_files/datasets/scenario0/rep0.csv", package = "SPMIX")
w_filename = system.file("input_files/prox_matrix.csv", package = "SPMIX")
params_filename = system.file("input_files/sampler_params.asciipb", package = "SPMIX")

data_filename_fake = "/home/m_gianella/rep0.csv"
w_filename_fake = "/home/m_gianella/prox_matrix.csv"
params_filename_fake = "/home/m_gianella/sampler_params.asciipb"

# Setting type-ready inputs (to check proper behaviour of the sampler)
data_obj <- ReadDataFromCSV(data_filename)
w_obj <- ReadMatrixFromCSV(w_filename)
readProtoFiles(file = system.file("proto/sampler_params.proto", package = "SPMIX"))
params_obj <- RProtoBuf::readASCII(SamplerParams, file(params_filename))

# Executing Sampler
out <- SPMIXSampler(burnin, niter, thin, data_filename, w_filename, params_filename, type = "no_rjmcmc")

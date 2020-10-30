## This is an example script to understand how to proper set parameters using RProtoBuf package ##

# Required libraries (SPMIX load RProtoBuf automatically)
library("SPMIX")

# Read proto file
readProtoFiles(file = system.file("proto/sampler_params.proto", package = "SPMIX"))

# Instanciate a new SamplerParams object
params <- new(SamplerParams)

# Inspect the fields of params
names(params)

# Set values for the fields (MANUALLY)
params$num_components <- 5
params$p0_params <- new(SamplerParams.NormalGammaParams, mu0 = 0, a = 2, b = 2, lam_ = 0.1)
params$rho_params <- new(SamplerParams.BetaParams, a = 1, b = 5)
params$sigma_params$inv_wishart$nu <- 100; params$sigma_params$inv_wishart$identity <- TRUE
params$mtilde_sigmasq <- 5

# Directly define a SamplerParams message instanciated from an ASCII file (AUTOMATICALLY)
params_fromfile <- readASCII(SamplerParams,
                             file(system.file("input_files/sampler_params.asciipb", package = "SPMIX")))
all.equal(params, params_fromfile) # TRUE, clearly!

# Reset descriptor (you need to reload proto files after this call)
resetDescriptorPool()

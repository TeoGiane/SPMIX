library("SPMIX")

# Setting input filenames
data_filename = system.file("input_files/datasets/scenario0/rep0.csv", package = "SPMIX")
w_filename = system.file("input_files/prox_matrix.csv", package = "SPMIX")
paramsrj_filename = system.file("input_files/rjsampler_params.asciipb", package = "SPMIX")
paramsfix_filename = system.file("input_files/sampler_params.asciipb", package = "SPMIX")
options_filename = system.file("input_files/optimization_options.asciipb", package = "SPMIX")

# Parsing files
data_obj <- ReadDataFromCSV(data_filename)
w_obj <- ReadMatrixFromCSV(w_filename)
readProtoFiles(file = system.file("proto/sampler_params.proto", package = "SPMIX"))
paramsrj_obj <- RProtoBuf::readASCII(SamplerParams, file(paramsrj_filename))
paramsfix_obj <- RProtoBuf::readASCII(SamplerParams, file(paramsfix_filename))
readProtoFiles(file = system.file("proto/optimization_options.proto", package = "SPMIX"))
options_obj <- RProtoBuf::readASCII(OptimOptions, file(options_filename))

# Running test function
Samplers_test(data_obj, w_obj, paramsrj_obj, paramsfix_obj, options_obj)

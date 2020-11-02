library("SPMIX")

# Getting options from file
options_filename <- system.file("input_files/optimization_options.asciipb", package = "SPMIX")
options <- RProtoBuf::readASCII(OptimOptions, file(options_filename))

# Running test
grad_ascent_test(options)

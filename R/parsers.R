## Internals - Parsers for samplers input to check and manage data types ##

###########################################################################
# Maybe Print to File -----------------------------------------------------
# If maybe_proto is a file, returns the file name. If maybe_proto is a string
# representing a message, prints the message to a file and returns its path.
maybe_print_to_file <- function(maybe_proto, proto_name = NULL, out_dir = NULL) {
  if(file.exists(maybe_proto)){
    return(maybe_proto)
  }
  proto_file = sprintf("%s/%s.asciipb", out_dir, proto_name)
  write(maybe_proto, file = proto_file)
  return(proto_file)
}

###########################################################################


###########################################################################
# Data Parser -------------------------------------------------------------
parseData <- function(data) {
  # Checking if data is given or needs to be read from file
  if(typeof(data) == "character") {
    cat("Data are provided as a path to a csv file\n")
    data_in <- ReadDataFromCSV(data)
    # Return if the input filepath does not exist
    if (all(is.na(data_in)))
      stop("Input file for 'data' does not exist.")
  } else if ( typeof(data)=="list" && prod(sapply(data, function(x) return(typeof(x)=="double" && length(x)>0))) ) {
    cat("Data are provided as a list of numeric vectors\n")
    data_in <- data
  } else {
    stop("Input parameter 'data' is of unknown type.")
  }
  # Return the parsed data structure for samplers
  return(data_in)
}

###########################################################################

###########################################################################
# W parser ----------------------------------------------------------------
parseW <- function(W) {
  # Checking if W is given or needs to be read from file
  if(typeof(W) == "character") {
    cat("Proximity Matrix is provided as a path to a csv file\n")
    W_in <- SPMIX::ReadMatrixFromCSV(W);
    # Return if the input filepath does not exist
    if (all(is.na(W_in)))
      stop("Input file for 'W' does not exist.")
  } else if ( typeof(W)=="double" && any(is(W)=="matrix") ) {
    cat("Proximity Matrix is provided as a matrix of double\n")
    W_in <- W
  } else {
    stop("Input parameter 'W' is of unknown type.")
  }
  # Return the parsed W structure for samplers
  return(W_in)
}

###########################################################################

###########################################################################
# Parameters parser -------------------------------------------------------
parseParams <- function(params, out_dir = NULL) {
  # Create file where to store serialized params
  serialized_params_file <- sprintf("%s/sampler_params.bin", out_dir)
  # Checking if params is given or needs to be read from file
  if(typeof(params) == "character") {
    cat("Hyperparameters are provided as a path to an asciipb file\n")
    # Check if file exists
    if(!file.exists(params))
      stop("Input file does not exist.")
    # Read ASCII file
    cat("readParamsfromASCII ... ")
    RProtoBuf::readProtoFiles(file = system.file("proto/sampler_params.proto", package = "SPMIX"))
    parsed_params <- RProtoBuf::readASCII(spmix.SamplerParams, file(params))
    mcmc_type <- ifelse(parsed_params$num_components$has("shifted_poisson_prior"), "rjmcmc", "no_rjmcmc")
    RProtoBuf::serialize(parsed_params, serialized_params_file)
    cat("done!\n")
  } else if ( is(params)=="Message" && params@type=="spmix.SamplerParams" ) {
    cat("Hyperparameters are provided as an RProtoBuf::Message\n")
    mcmc_type <- ifelse(params$num_components$has("shifted_poisson_prior"), "rjmcmc", "no_rjmcmc")
    RProtoBuf::serialize(params, serialized_params_file)
  } else {
    stop("Input parameter 'params' is of unknown type.")
  }
  # Return the serialized params file path for samplers
  returned_list <- list("filepath" = serialized_params_file, "mcmc_type" = mcmc_type)
  return(returned_list)
}

###########################################################################

###########################################################################
# Option parser -----------------------------------------------------------
parseOptions <- function(options, out_dir = NULL) {
  # Create file where to store serialized options
  serialized_options_file <- sprintf("%s/optim_options.bin", out_dir)
  # Checking if options is NULL, given or needs to be read from file
  if (is.null(options)) {
    cat("Optimization Options required but not given: setting default values ... ")
    RProtoBuf::readProtoFiles(file = system.file("proto/optimization_options.proto", package = "SPMIX"))
    RProtoBuf::serialize(RProtoBuf::new(spmix.OptimOptions, max_iter = 20, tol = 1e-6, jump_every = 1), serialized_options_file)
    cat("done!\n")
  } else if(typeof(options) == "character") {
    cat("Optimization Options are provided as a path to an asciipb file\n")
    # Check if file exists
    if(!file.exists(options))
      stop("Input file does not exist.")
    # Read ASCII file
    cat("readOptimOptionsfromASCII ... ")
    RProtoBuf::readProtoFiles(file = system.file("proto/optimization_options.proto", package = "SPMIX"))
    RProtoBuf::serialize(RProtoBuf::readASCII(spmix.OptimOptions, file(options)), serialized_options_file)
    cat("done!\n")
  } else if ( is(options)=="Message" && options@type=="spmix.OptimOptions" ) {
    cat("Optimization Options are provided as an RProtoBuf::Message\n")
    RProtoBuf::serialize(options, serialized_options_file)
  } else {
    stop("Input parameter 'options' is of unknown type.")
  }
  # Return the serialized options file path for samplers
  return(serialized_options_file)
}

###########################################################################

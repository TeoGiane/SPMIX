## Internals - Parsers for samplers input to check and manage data types ##

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
parseParams <- function(params) {
  # Checking if params is given or needs to be read from file
  if(typeof(params) == "character") {
    cat("Hyperparameters are provided as a path to an asciipb file\n")
    # Check if file exists
    if(!file.exists(params))
      stop("Input file does not exist.")
    # Read ASCII file
    cat("readParamsfromASCII ... ")
    RProtoBuf::readProtoFiles(file = system.file("proto/sampler_params.proto", package = "SPMIX"))
    params_in <- RProtoBuf::toString(RProtoBuf::readASCII(SamplerParams, file(params)))
    cat("done!\n")
  } else if ( is(params)=="Message" && params@type=="SamplerParams" ) {
    cat("Hyperparameters are provided as an RProtoBuf::Message\n")
    params_in <- RProtoBuf::toString(params)
  } else {
    stop("Input parameter 'params' is of unknown type.")
  }

  # Return the parsed params structure for samplers
  return(params_in)
}

###########################################################################

###########################################################################
# Option parser -----------------------------------------------------------
parseOptions <- function(options) {
  # Checking if options is NULL, given or needs to be read from file
  if (is.null(options)) {
    cat("Optimization Options required but not given: setting default values ... ")
    RProtoBuf::readProtoFiles(file = system.file("proto/optimization_options.proto", package = "SPMIX"))
    options_in <- RProtoBuf::toString(RProtoBuf::new(OptimOptions, max_iter = 1000, tol = 1e-6))
    cat("done!\n")
  } else if(typeof(options) == "character") {
    cat("Optimization Options are provided as a path to an asciipb file\n")
    # Check if file exists
    if(!file.exists(options)){}
    stop("Input file does not exist.")
    # Read ASCII file
    cat("readOptimOptionsfromASCII ... ")
    RProtoBuf::readProtoFiles(file = system.file("proto/optimization_options.proto", package = "SPMIX"))
    options_in <- RProtoBuf::toString(RProtoBuf::readASCII(OptimOptions, file(options)))
    cat("done!\n")
  } else if ( is(options)=="Message" && options@type=="OptimOptions" ) {
    cat("Optimization Options are provided as an RProtoBuf::Message\n")
    options_in <- RProtoBuf::toString(options)
  } else {
    stop("Input parameter 'options' is of unknown type.")
  }

  # Return the parsed options structure for samplers
  return(options_in)
}

###########################################################################

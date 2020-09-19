#' Spatial Mixture Model Sampler
#'
#' This function provides the sampling algorithm for Spatial Mixture Models data. A seguire, se tutto funziona,
#' una bella documentazione.
#' @param burnin blabla
#' @param niter blabla
#' @param thin blabla
#' @param data blabla
#' @param W blabla
#' @param params blabla
#' @param cov blabla
#' @param display_progress blabla
#' @return bleble
#'
#' @export
SPMIX_sampler <- function(burnin, niter, thin, data, W, params, cov = list(), display_progress = TRUE) {

  # Checking if data is given or needs to be read from file
  if(typeof(data) == "character") {
    cat("Data are provided as a path to a csv file\n")
    data_in <- readDataFromCSV(data)
    # Return if the input filepath does not exist
    if (all(is.na(data_in)))
      return()
  } else if ( typeof(data)=="list" && sapply(data, function(x) return(typeof(x)=="double" && length(x)>0)) ) {
    cat("Data are provided as a list of numeric vectors\n")
    data_in <- data
  } else {
    cat("ERROR: input parameter 'data' is of unknown type\n")
    return()
  }

  # Checking if W is given or needs to be read from file
  if(typeof(W) == "character") {
    cat("Proximity Matrix is provided as a path to a csv file\n")
    W_in <- SPMIX::readMatrixFromCSV(W);
    # Return if the input filepath does not exist
    if (all(is.na(W_in)))
      return()
  } else if ( typeof(W)=="double" && any(is(W)=="matrix") ) {
    cat("Proximity Matrix is provided as a matrix of double\n")
    W_in <- W
  } else {
    cat("ERROR: input parameter 'W' is of unknown type\n")
    return()
  }

  # Checking if params is given or needs to be read from file
  if(typeof(params) == "character") {
    cat("Hyperparameters are provided as a path to an asciipb file\n")

    # Check if file exists
    if(!file.exists(params))
      stop("Input file does not exist.")

    # Read ASCII file
    cat("readParamsfromASCII ... ")
    RProtoBuf::readProtoFiles(file = system.file("proto/sampler_params.proto", package = "SPMIX"))
    params_in <- RProtoBuf::readASCII(SamplerParams, file(params))
    cat("done!\n")
  } else if ( is(params)=="Message" && params@type=="SamplerParams" ) {
    cat("Hyperparameters are provided as an RProtoBuf::Message\n")
    params_in <- params
  } else {
    cat("ERROR: input parameter 'params' is of unknown type\n")
    return()
  }

 # Calling the Sampler
 output <- SPMIX:::runSpatialSampler(burnin, niter, thin, data_in, W_in, params_in, cov, display_progress)
 return(output)
}

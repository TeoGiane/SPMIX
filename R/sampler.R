#' Spatial Mixture Model Sampler
#'
#' Runs the Gibbs sampler for the SPMIX model for a total of burnin + niter iterations,
#' discarding the first 'burnin' ones and keeping in memory only one every 'thin' iterations.
#'
#' @usage SPMIX_sampler(burnin, niter, thin, data, W, params,
#'     cov = list(), display_progress = TRUE)
#'
#' @param burnin Integer, the number of steps of the burnin phase.
#' @param niter Integer, the number of steps to run *after* the burnin,
#' @param thin Integer, is the thinning parameter, hence only one every 'thin' iterations will be saved
#' during the sampling phase.
#' @param data The data that needs to be fitted by the model. Data can be passed either as string, thus indicating
#' the path to a \code{.csv} file representing the data, or as a list of vectors. In this case, data[[i]] represents
#' the vector of data in location i.
#' @param W The proximity matrix between different locations. W can be passed as string, i.e. a string containing
#' the path to a \code{.csv} file storing the proximity matrix, or as a simple R matrix.
#' @param params The sampler parameters that needs to be provided as input. params can be passed as a string
#' containing the path to an \code{.asciipb} file or as an \code{S4::Message} class of type SamplerParams, i.e. a
#' Google Protocol Buffer Message available in the package and interfaced to R through \code{\link{RProtoBuf}}.
#' @param cov a list of vectors that represents covariates. As default value, is an empty list.
#' In case this input parameter is not empty, the sampler performs a regression on these covariates.
#' @param display_progress Boolean, it allows you display on the console the progress bar during burn-in
#' and sampling phase. As default value, is set to TRUE.
#'
#' @return A list of raw vectors where the i-th element is the i-th saved draw In order to reduce the space
#' occupied by these draws, data are serialized through Google Protocol Buffers serialization procedure.
#' Each state can be easily deserialized in R using the \code{\link{unserialize_proto}} function in this package.
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

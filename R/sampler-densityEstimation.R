#' \code{SPMIX} Sampler for Density Estimation
#'
#' \loadmathjax Runs the Gibbs samplers for the \code{SPMIX} model for a total of burnin + niter iterations,
#' discarding the first 'burnin' ones and keeping in memory only one every 'thin' iterations.
#'
#' @param burnin Integer, the number of steps of the burnin phase.
#' @param niter Integer, the number of steps to run *after* the burnin.
#' @param thin Integer, is the thinning parameter, hence only one every 'thin' iterations will be saved
#' during the sampling phase.
#' @param data The data that needs to be fitted by the model. Data can be passed either as string, thus indicating
#' the path to a \code{.csv} file representing the data, or as a list of vectors. In this case, the
#' \mjseqn{i}-th element of the data list represents the vector of data assigned to the \mjseqn{i}-th location.
#' @param W The proximity matrix between different locations. It can be passed as string, i.e. a string containing
#' the path to a \code{.csv} file storing the proximity matrix, or as a simple R matrix.
#' @param params The sampler parameters that needs to be provided as input. params can be passed as a string
#' containing the path to an \code{.asciipb} file or as an \code{S4::Message} class of type SamplerParams, i.e. a
#' Google Protocol Buffer Message available in the package and interfaced to R through \code{\link{RProtoBuf}}.
#' @param type A string identifying the type of sampler to run. If type is "rjmcmc", the algorithm will run
#' the spatial mixture sampler putting a prior on the number of components \mjseqn{H}.
#' The default value is "no_rjmcmc", which samples from the spatial mixture model with a fixed number
#' of components.
#' @param options The sampler optimization options used in the execution of the reversible jump sampler.
#' Default value is set to \code{NULL} and in case type is "rjmcmc", a \code{S4::Message} object of type
#' OptimOptions is istanciated with default values. In order to override the default values, options can be
#' passed either as a string, representing the path to an \code{.asciipb} file, or as a \code{S4::Message}
#' object of type OptimOptions generated via \code{\link{RProtoBuf}}.
#' @param display_progress Boolean, it allows you display on the console the progress bar during burn-in
#' and sampling phase. As default value, it is set to TRUE.
#'
#' @return A list of raw vectors where the \mjseqn{i}-th element is the \mjseqn{i}-th saved draw.
#' In order to reduce the space occupied by these draws, data are serialized using Google Protocol Buffers.
#' Each state can be easily deserialized in R using the \code{\link{DeserializeSPMIXProto}} function of this package.
#'
#' @export
Sampler.DensityEstimation <- function(burnin, niter, thin, data, W, params,
                                      type = "no_rjmcmc", options = NULL, display_progress = TRUE) {

  # Create .asciipb files in temporary directory if needed
  out_dir = tempdir()
  if(typeof(params) == "character") { params <- maybe_print_to_file(params, 'sampler_params', out_dir) }
  if(typeof(options) == "character") { options <- maybe_print_to_file(options, 'optim_options', out_dir) }

  # Check and parse of input members
  data_in <- parseData(data)
  W_in <- parseW(W)
  params_in <- parseParams(params)

  # Check sampler type and run
  if (type == "no_rjmcmc") {
    output <- SPMIX:::runSpatialSampler(burnin, niter, thin, data_in, W_in, params_in, list(), FALSE, display_progress)
  } else if (type == "rjmcmc") {
    options_in <- parseOptions(options)
    output <- SPMIX:::runSpatialRJSampler(burnin, niter, thin, data_in, W_in, params_in, list(), options_in, FALSE, display_progress)
  } else {
    stop("Input parameter 'type' is of unknown type.")
  }

  # Remove temporary files amd return
  unlink(paste0(out_dir,"/*.asciipb"))
  return(output)
}

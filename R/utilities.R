#' Deserialize messages of the \code{SPMIX} package
#'
#' This funciton is a wrapper to deserialize raw vectors using the \code{.proto} files available in this package.
#' We rely on Google Protocol Buffers for the serialization procedure and on \code{\link{RProtoBuf}} package to
#' provide and easy-to-use interface for \code{R} users
#'
#' @param message_type A string containing the name of the Message in which the serialized message will be converted
#' @param raw_vector A vector of type \code{raw}. The message to be unserialized
#' @return An object of class \code{RProtoBuf::Message}, which stores the unserialized message and can be manipulated
#' using the \link{RProtoBuf} package.
#'
#' @export
DeserializeSPMIXProto <- function(message_type, raw_vector) {

  # Check Message Descriptor
  if (message_type == "EigenMatrix") {
    RProtoBuf::readProtoFiles(system.file("/proto/eigen.proto", package = "SPMIX"))
  } else if (message_type == "SamplerParams") {
    RProtoBuf::readProtoFiles(system.file("/proto/sampler_params.proto", package = "SPMIX"))
  } else if (message_type == "UnivariateState") {
    RProtoBuf::readProtoFiles(system.file("/proto/univariate_mixture_state.proto", package = "SPMIX"))
  } else if (message_type == "OptimOptions") {
    RProtoBuf::readProtoFiles(system.file("/proto/optimization_options.proto", package = "SPMIX"))
  } else {
    stop("Input 'message_type' is of uknown type")
  }

  # Read state from proper descriptor
  state <- RProtoBuf::read(get(message_type), raw_vector)
  return(state)
}

#' Compute the chain of the estimated posterior densities over a grid of points
#'
#' \loadmathjax This utility takes as input the deserialized output of the samplers
#' (via \code{\link{DeserializeSPMIXProto}}) and compute the chain of posterior estimates of the densities for each
#' areal location \mjseqn{i=1,\dots,I}.
#'
#' @param deserialized_chains A list of \code{RProtoBuf::Message} which stores the deserialized output of
#' the sampler (either with a fixed or a variable number of components).
#' @param x_grid A numeric vector representing the grid of points on which the estimated densities will be evaluated.
#' @param verbose A bool. If \code{TRUE}, prints the progress of the computation.
#' @return A list of \mjseqn{I} elements, where the \mjseqn{i}-th element is a matrix that represent, by row,
#' the estimated density, the lower and the upper bounds the corresponding credibility interval
#' (if such quantity is required) evaluated over \code{x_grid}.
#'
#' @export
ComputeDensities <- function(deserialized_chains, x_grid, verbose = FALSE) {
  # Check input type for deserialized_chains
  if(!is.list(deserialized_chains)){
    stop("'deserialized_chains' input is not of type 'list'.")
  }
  if(all(!sapply(deserialized_chains,
                 function(x) return(typeof(x)=="S4" &&
                                    class(x)=="Message" &&
                                    x@type == "UnivariateState")))) {
    stop("'deserialized_chains' is a list of wrong type.")
  }
  # Elicit numGroups
  numGroups <- length(deserialized_chains[[1]]$groupParams)
  # Check x_grid
  if(!is.numeric(x_grid)){
    stop("'x_grid' input is not of type 'numeric'.")
  }
  # Extract necessary chains
  H_chain <- sapply(deserialized_chains, function(x) x$num_components)
  means_chain <- lapply(deserialized_chains, function(x) sapply(x$atoms, function(x) x$mean))
  stdev_chain <- lapply(deserialized_chains, function(x) sapply(x$atoms, function(x) x$stdev))
  # Compute Estimated
  estimated_densities <- list()
  for (i in 1:numGroups) {
    weights_chain <- lapply(deserialized_chains, function(x) x$groupParams[[i]]$weights)
    est_dens_mat <- matrix(0,length(deserialized_chains),length(x_grid))
    for (j in 1:length(deserialized_chains)) {
      xgrid_expand <- t(rbind(replicate(H_chain[j], x_grid, simplify = "matrix")))
      est_dens_mat[j,] <- t(as.matrix(weights_chain[[j]])) %*% dnorm(xgrid_expand,means_chain[[j]],stdev_chain[[j]])
    }
    # Build the element of the return list
    estimated_densities[[i]] <- est_dens_mat
    # Print progress (if verbose)
    if(verbose){
      cat(sprintf("\rProcessed Area: %d / %d", i, numGroups))
    }
  }
  # New line
  if(verbose) {
    cat("\n")
  }
  # Return list
  return(estimated_densities)
}


#' Compute the chain of the estimated posterior log-likelihood for each data point
#'
#' \loadmathjax This utility takes as input the deserialized output of the samplers
#' (via \code{\link{DeserializeSPMIXProto}}) and compute the chain of posterior log-likelihood
#' for each data point
#'
#' @param data The data that needs to be fitted by the model. Data are passed as a list of vectors, whose
#' \mjseqn{i}-th element represents the vector of data assigned to the \mjseqn{i}-th location.
#' @param deserialized_chains A list of \code{RProtoBuf::Message} which stores the deserialized output of
#' the sampler (either with a fixed or a variable number of components).
#' @param verbose A bool. If \code{TRUE}, prints the progress of the computation.
#' @return A \mjseqn{T \times N} matrix, \mjseqn{T} being the number of iterations
#' of the MCMC chain and \mjseqn{N} the total number of data points.
#' Element \mjseqn{t,n} of the matrix is the log-likelihood of the \mjseqn{n}-th
#' data point at \mjseqn{t}-th iteration.
#'
#' @export
ComputePosteriorLPDF <- function(data, deserialized_chains, verbose = FALSE) {
  # Check input types
  data <- parseData(data)
  if(!is.list(deserialized_chains)){
    stop("'deserialized_chains' input is not of type 'list'.")
  }
  if(all(!sapply(deserialized_chains,
                 function(x) return(typeof(x)=="S4" &&
                                    class(x)=="Message" &&
                                    x@type == "UnivariateState")))) {
    stop("'deserialized_chains' is a list of wrong type.")
  }
  # Define buffer
  out <- matrix(NA, length(deserialized_chains), 0)
  # Get means and stadard deviations
  means <- lapply(deserialized_chains, function(x){ sapply(x$atoms, function(a){a$mean}) })
  stdevs <- lapply(deserialized_chains, function(x){ sapply(x$atoms, function(a){a$stdev}) })
  for (i in 1:length(data)) {
    # Get cluster allocations in area i
    clus_allocs <- t(sapply(deserialized_chains, function(x){x$groupParams[[i]]$cluster_allocs}))
    # Compute posterior log-likelihood for each datum in area i
    out_area <- matrix(NA, length(deserialized_chains), length(data[[i]]))
    for (j in 1:length(data[[i]])) {
      mean_vector <- sapply(1:length(deserialized_chains), function(l){means[[l]][(1L+clus_allocs[l,j])]})
      stdev_vector <- sapply(1:length(deserialized_chains), function(l){stdevs[[l]][(1L+clus_allocs[l,j])]})
      out_area[,j] <- dnorm(data[[i]][j], mean = mean_vector, sd = stdev_vector, log = T)
    }
    # Log if verbose
    if(verbose){
      cat(sprintf("\rArea %g/%g", i,length(data)))
    }
    # Stack plpdf in general buffer
    out <- cbind(out, out_area)
  }
  # Return output matrix
  return(out)
}

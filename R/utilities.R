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

#' Compute the estimated posterior densities
#'
#' \loadmathjax This utility takes as input the deserialized output of the samplers
#' (via \code{\link{DeserializeSPMIXProto}}) and compute the posterior estimates of the densities for each
#' areal location \mjseqn{i=1,\dots,I}.
#'
#' @param deserialized_chains A list of \code{RProtoBuf::Message} which stores the deserialized output of
#' the sampler (either with a fixed or a variable number of components).
#' @param N An integer representing the number of points of the grid on which the estimated density will be computed.
#' @param ranges A matrix of dimesions \mjseqn{2 \times I}, where \mjseqn{I} is the number of areal locations. This
#' matrix represents the extrema of the grid on which the density estimation will be performed.
#' @param alpha A double in \code{[0,1]} that indicates the level at which compute credibility bands for the
#' estimated densities. The deafult value is set to \code{NULL}, thus meaning that the bounds are not computed.
#' @param names A vector of strings of length \mjseqn{I}, optional input that sets the names of the output elements.
#' Default value is set to \code{NULL}, thus no name will be attributed.
#' @return A list of \mjseqn{I} elements, where the \mjseqn{i}-th element is a matrix that represent, by row,
#' the estimated density, the lower and the upper bounds the corresponding credibility interval
#' (if such quantity is required) over a fixed grid of \code{N} points between \code{ranges[1,i]} and \code{ranges[2,i]}.
#'
#' @export
ComputeDensities <- function(deserialized_chains, N, ranges, alpha = NULL, names = NULL) {

  # Check input type for deserialized_chains
  if(!is.list(deserialized_chains)){
    stop("'deserialized_chains' input is not of type 'list'.")
  }
  if(all(!sapply(deserialized_chains, function(x) return(typeof(x)=="S4" &&
                                                         class(x)=="Message" &&
                                                         x@type == "UnivariateState")))) {
    stop("'deserialized_chains' is a list of wrong type.")
  }

  # Elicit numGroups
  numGroups <- length(deserialized_chains[[1]]$groupParams)

  # Check ranges
  if(!is.matrix(ranges)){
    stop("'ranges' input is not of type 'matrix'.")
  }
  if(all(dim(ranges) != c(2,numGroups))){
    stop("'ranges' dimensions does not match.")
  }

  # Cast to integer N
  N <- as.integer(N)

  # Extract necessary chains
  H_chain <- sapply(deserialized_chains, function(x) x$num_components)
  means_chain <- lapply(deserialized_chains, function(x) sapply(x$atoms, function(x) x$mean))
  stdev_chain <- lapply(deserialized_chains, function(x) sapply(x$atoms, function(x) x$stdev))

  # Compute Estimated
  estimated_densities <- list()
  for (i in 1:numGroups) {
    weights_chain <- lapply(deserialized_chains, function(x) x$groupParams[[i]]$weights)
    x_grid <- seq(ranges[1,i], ranges[2,i], length.out = N)
    est_dens_mat <- matrix(0,length(deserialized_chains),length(x_grid))
    for (j in 1:length(deserialized_chains)) {
      xgrid_expand <- t(rbind(replicate(H_chain[j], x_grid, simplify = "matrix")))
      est_dens_mat[j,] <- t(as.matrix(weights_chain[[j]])) %*% dnorm(xgrid_expand,means_chain[[j]],stdev_chain[[j]])
    }

    # Build the element of the return list
    est_dens <- colMeans(est_dens_mat)
    if (!is.null(alpha)) {
      up_bound <- apply(est_dens_mat, 2, quantile, probs=1-alpha/2)
      low_bound <- apply(est_dens_mat, 2, quantile, probs=alpha/2)
      estimated_densities[[i]] <- rbind("est"=est_dens,"up"=up_bound,"low"=low_bound)
    } else {
      estimated_densities[[i]] <- rbind("est"=est_dens)
    }
  }

  # Set names if passed
  if (!is.null(names)) {
    names(estimated_densities) <- names
  }

  # Return list
  return(estimated_densities)
}

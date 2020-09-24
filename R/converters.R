#' Unserialize raw vectors
#'
#' This funciton is a wrapper to unserialize raw vectors that comes from the serialization of
#' messages using Google Protocol Buffer API using the proto files available in the \code{SPMIX} package.
#'
#' @param message_type A string containing the name of the Message in which the serialized message will be converted
#' @param raw_vector A vector of type \code{raw}. The message to be unserialized
#' @return An object of class \code{RProtoBuf::Message}, which stores the unserialized message and can be manipulated
#' using \link{RProtoBuf}.
#'
#' @export
unserialize_proto <- function(message_type, raw_vector) {

  # Check Message Descriptor
  if (message_type == "EigenMatrix") {
    RProtoBuf::readProtoFiles(system.file("/proto/eigen.proto", package = "SPMIX"))
  } else if (message_type == "SamplerParams") {
    RProtoBuf::readProtoFiles(system.file("/proto/sampler_params.proto", package = "SPMIX"))
  } else if (message_type == "UnivariateState") {
    RProtoBuf::readProtoFiles(system.file("/proto/univariate_mixture_state.proto", package = "SPMIX"))
  } else {
    stop("Input 'message_type' is of uknown type")
  }

  # Read state from proper descriptor
  state <- RProtoBuf::read(get(message_type), raw_vector)
  return(state)
}

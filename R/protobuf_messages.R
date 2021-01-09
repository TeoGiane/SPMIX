#' Parameters of the \code{SPMIX} samplers
#'
#' \loadmathjax This protobuf message provides the data structure to encode the prior hyper-parameters of the
#' \code{SPMIX} samplers available in this package. Using \code{\link{RProtoBuf}} package one can easily
#' retrieve the descriptor and build it from scratch or either parse it from a \code{.asciipb} file.
#'
#' @format A Protobuf::Message with the following fields:
#' \describe{
#'   \item{num_components}{An integer setting the initial number of components.}
#'   \item{p0_params}{A \code{NormalGammaParams} message, that can be simply imagined as a named tuple
#'   storing the prior paramteres for the atom base measure \mjseqn{ P_{0} }, in our case a Normal-Inverse-Gamma
#'   measure defined as \mjseqn{P_{0}( \mu , \sigma ^2) \sim \mathcal{N}( \mu \lvert \sigma ^2; \mu_{0}, \lambda^{-1}
#'   \sigma ^2) \times \operatorname{IG}(\sigma^2; a, b) }. Hence, this message stores four double, i.e.
#'   \mjseqn{(\mu_{0}, a, b, \lambda)}.}
#'   \item{rho_params}{A \code{BetaParams} message that stores the prior parameters for \mjseqn{\rho},
#'   whose prior is a Beta distribution of parameters \mjseqn{(\alpha,\beta)}.}
#'   \item{sigma_params}{a \code{VarianceParams} message that stores the prior parameters for the variance
#'   of the transformed weights. Note that at different models correspond different definitions of such parameter.
#'   In case of fixed dimension, we have a Inverse Wishart distribution for a full matrix \mjseqn{\boldsymbol{\Sigma}},
#'   while in case of varying dimension, we assume \mjseqn{\boldsymbol{\Sigma}=\sigma^2\mathbf{I}} and put
#'   an Inverse Gamma prior on \mjseqn{\sigma^2}. In order to manage both, \code{VarianceParams} is structured
#'   using a \code{oneof} field, thus making it capable of storing both parameters of an Inverse Wishart prior
#'   \mjseqn{\textstyle(\nu,\mathbf{I})} and an Inverse Gamma prior \mjseqn{\textstyle(a_{\sigma^2},b_{\sigma^2})}.}
#' }
#'
#' @examples
#' # Retrieve message descriptor
#' RProtoBuf::readProtoFiles(system.file("proto/sampler_params.proto", package = "SPMIX"))
#' names(SamplerParams)
#'
#' # Build a new message from scratch
#' params <- RProtoBuf::new(SamplerParams)
#'
#' # Read from file (NOT RUN)
#' params <- RProtoBuf::readASCII(SamplerParams, file("path/to/file.asciipb"))
#'
#' @name SamplerParams [Protobuf::Message]
NULL

#' Options for optimization algorithms
#'
#' This protobuf message encodes the options required to the optimization algorithms implemented
#' in the package and used in the reversible jump sampler. Using \code{\link{RProtoBuf}} package one can easily
#' retrieve the descriptor and build it from scratch or either parse it from a \code{.asciipb} file.
#'
#' @format A Protobuf::Message with the following fields:
#' \describe{
#'   \item{max_iter}{An integer storing the maximum number of iterations that the optimization algorithm
#'   have to perform.}
#'   \item{tol}{A double storing the minimum tolerance used to evaluate if the optimization algorithm
#'   has reached convergence.}
#' }
#'
#' @examples
#' # Retrieve message descriptor
#' RProtoBuf::readProtoFiles(system.file("proto/optimization_option.proto", package = "SPMIX"))
#' names(OptimOptions)
#'
#' # Build a new message from scratch
#' params <- RProtoBuf::new(OptimOptions)
#'
#' # Read from file (NOT RUN)
#' params <- RProtoBuf::readASCII(OptimOptions, file("path/to/file.asciipb"))
#'
#' @name OptimOptions [Protobuf::Message]
NULL

#' Draws of the \code{SPMIX} samplers
#'
#' \loadmathjax This protobuf message encodes the structure of the state of the samplers at a given iteration.
#' Using \code{\link{RProtoBuf}} package one can easily retrieve the descriptor and build it from scratch
#' or either parse it from a \code{.asciipb} file.
#'
#' @format A Protobuf::Message with the following fields:
#' \describe{
#'   \item{num_components}{An integer setting the number of components of the mixture.}
#'   \item{atoms}{A repeated field of type \code{GroupParams}, a message that stores, for each areal location
#'   \mjseqn{i=1,\dots,I}, the vector of weights \mjseqn{w_i} and the allocation variables \mjseqn{s_{ij}},
#'   \mjseqn{j=1,\dots,N_{i}}.}
#'   \item{rho}{A double storing the current value of \mjseqn{\rho}.}
#'   \item{sigma}{An \code{EigenMatrix} message storing the matrix \mjseqn{\boldsymbol{\Sigma}}. The latter message
#'   simply stores two integers for the dimensions of the matrix and a repeated field of double to store data
#'   continuousy, as *Eigen* does.}
#'   \item{regression_coefficients}{A repeated packed field that stores the regression coefficients for each datum.
#'   This field is set only in case covariates are passed as input in the sampler, so in case a Bayesian linear
#'   regression is performed.}
#' }
#'
#' @examples
#' # Retrieve message descriptor
#' RProtoBuf::readProtoFiles(system.file("proto/univariate_mixture_state.proto", package = "SPMIX"))
#' names(UnivariateState)
#'
#' # Build a new message from scratch
#' params <- RProtoBuf::new(UnivariateState)
#'
#' # Read from file (NOT RUN)
#' params <- RProtoBuf::readASCII(UnivariateState, file("path/to/file.asciipb"))
#'
#' @name UnivariateState [Protobuf::Message]
NULL

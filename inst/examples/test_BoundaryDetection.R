## Test for boundary detection for the Reversible Jump Sampler ##
# The scenario is described in Sec. 6.1 of Beraha et al. (2020)

# Required libraries
library("SPMIX")
library("ggplot2")

# Helper functions
DataRecast <- function(x, y, labels = row.names(y)) {
  rows <- dim(y)[1]; cols <- dim(y)[2]
  out <- data.frame()
  for (i in 1:rows) {
    tmp <- data.frame(x, y[i,], as.factor(rep(labels[i],cols)));
    names(tmp) <- c("Grid", "Value", "Density"); row.names(tmp) <- NULL
    out <- rbind(out,tmp)
  }
  return(out)
}

###########################################################################
# Data Generation ---------------------------------------------------------

# Generating data
I <- 6; Ns <- rep(100,I); set.seed(230196)
data <- list(); labels <- sprintf("g%d", 1:I)
data[[1]] <- metRology::rt.scaled(n=Ns[1], df=6, mean=-4, sd=1)
data[[2]] <- metRology::rt.scaled(n=Ns[2], df=6, mean=-4, sd=1)
data[[3]] <- sn::rsn(n=Ns[3], xi=4, omega=4, alpha=1); attributes(data[[3]]) <- NULL
data[[4]] <- sn::rsn(n=Ns[4], xi=4, omega=4, alpha=1); attributes(data[[4]]) <- NULL
data[[5]] <- rchisq(n=Ns[5], df=3)
data[[6]] <- rchisq(n=Ns[6], df=3)
names(data) <- labels

# Setting initial W
W <- matrix(1,I,I)

###########################################################################

###########################################################################
# Sampler execution -------------------------------------------------------

# Setting MCMC parameters
burnin = 5000
niter = 5000
thin = 2

# Grab input filenames
params_filename = system.file("input_files/rjsampler_params.asciipb", package = "SPMIX")

# Run Spatial sampler
out <- SPMIXSampler(burnin, niter, thin, data, W, params_filename, type = "rjmcmc")

###########################################################################

# Deserialization
chains <- sapply(out, function(x) DeserializeSPMIXProto("UnivariateState",x))
H_chain <- sapply(chains, function(x) x$num_components)
G_chain <- lapply(chains, function(x) matrix(x$G$data,x$G$rows,x$G$cols))

# reducing visited graphs with how many times they have been visited
Gunique <- unique(G_chain)
freq <- vector()
for (i in 1:length(Gunique)) {
  tmp <- sapply(G_chain, function(x) prod(x == Gunique[[i]]))
  freq[i] <- sum(tmp)
}

# plinks according to mean
plinks <- Reduce('+', G_chain)/length(G_chain)

# Computing estimated densities
data_ranges <- sapply(data, range)
estimated_densities <- ComputeDensities(chains, 500, data_ranges, labels)

# Computing true densities for comparison
true_densities <- list()
true_densities[[1]] <- metRology::dt.scaled(seq(data_ranges[1,1],data_ranges[2,1],length.out = 500), df=6, mean=-4, sd=1)
true_densities[[2]] <- metRology::dt.scaled(seq(data_ranges[1,2],data_ranges[2,2],length.out = 500), df=6, mean=-4, sd=1)
true_densities[[3]] <- sn::dsn(seq(data_ranges[1,3],data_ranges[2,3],length.out = 500), xi=4, omega=4, alpha=1); attributes(true_densities[[3]]) <- NULL
true_densities[[4]] <- sn::dsn(seq(data_ranges[1,4],data_ranges[2,4],length.out = 500), xi=4, omega=4, alpha=1); attributes(true_densities[[4]]) <- NULL
true_densities[[5]] <- dchisq(seq(data_ranges[1,5],data_ranges[2,5],length.out = 500), df=3)
true_densities[[6]] <- dchisq(seq(data_ranges[1,6],data_ranges[2,6],length.out = 500), df=3)
names(true_densities) <- labels


# Comparison plots between estimated and true densities in i-th area
plots_area <- list()
for (i in 1:I) {
  x <- seq(data_ranges[1,i], data_ranges[2,i], length.out = 500)
  y <- rbind(estimated_densities[[i]], true_densities[[i]]); row.names(y) <- c("Estimated", "True")
  tmp <- ggplot(data = DataRecast(x,y), aes(x=Grid, y=Value, group=Density, col=Density)) +
    geom_line(lwd=1) + theme(plot.title = element_text(face="bold", hjust = 0.5)) +
    ggtitle(paste0("Area ", i)) + theme(legend.position = "none")
  plots_area[[i]] <- tmp; rm(list=c('x','y','tmp'))
}
names(plots_area) <- labels

x11(height = 6, width = 8.27); gridExtra::grid.arrange(grobs=plots_area, nrow=2, ncol=3)

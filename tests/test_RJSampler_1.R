## First test for the Reversible Jump Sampler ##

# Required libraries
library("SPMIX")
library("ggplot2")

# Utilites
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
# Data Simulation ---------------------------------------------------------
# Generate data (1 location, mixture of 3 normals)
set.seed(230196)
ngroups <- 1; ncomponents <- 3; N <- 1000
means <- c(-4,1,5); std_devs <- c(1,1,1); weights <- c(1/6,3/6,2/6)
cluster_alloc <- sample(1:ncomponents, prob = weights, size = N, replace = T)
data <- list(); data[[1]] <- rnorm(N, mean = means[cluster_alloc], sd = std_devs[cluster_alloc])

# Generate Matrix W
W <- matrix(0, nrow = 1, ncol = 1, byrow = T)

###########################################################################

###########################################################################
# Sampler Run -------------------------------------------------------------
# Setting MCMC parameters
burninfull = 0;     burnin = 5000
niterfull = 10000;  niter = 5000
thinfull = 1;       thin = 2

# Grab input filenames
params_filename = system.file("input_files/rjsampler_params.asciipb", package = "SPMIX")

# Run Spatial sampler
out <- SPMIXSampler(burninfull, niterfull, thinfull, data, W, params_filename, type = "rjmcmc")

###########################################################################

###########################################################################
# Analyses and Visualization ----------------------------------------------

# Deserialization
chains <- sapply(out, function(x) DeserializeSPMIXProto("UnivariateState",x))
H_chain <- sapply(chains, function(x) x$num_components)

# Traceplot for the whole chain (no burnin, no thinning)
df <- data.frame("Iteration"=1:niterfull, "LowPoints"=H_chain-0.3, "UpPoints"=H_chain+0.3)
plot_traceH <- ggplot(data=df, aes(x=Iteration, y=LowPoints, xend=Iteration, yend=UpPoints)) +
  ylim(range(df[,-1])) + ylab("NumComponents") + geom_segment(lwd=0.1) +
  theme(plot.title = element_text(face="bold", hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  ggtitle("Traceplot of H")
rm(list='df')

x11(height = 3, width = 3); plot_traceH

# Reducing chain according to burnin, niter and thin params
fullout <- out
out <- fullout[burnin+1:length(fullout)]; out<-out[!sapply(out,is.null)]; out <- out[which((1:length(out)) %% thin == 1)]
chains <- sapply(out, function(x) DeserializeSPMIXProto("UnivariateState",x))
H_chain <- sapply(chains, function(x) x$num_components)

# Barplot of the estimated posterior for H
df <- as.data.frame(table(H_chain)/length(H_chain)); names(df) <- c("NumComponents", "Prob.")
plot_postH <- ggplot(data = df, aes(x=NumComponents, y=Prob.)) +
  geom_bar(stat="identity", color="steelblue", fill="lightblue") +#"white") +
  theme(plot.title = element_text(face="bold", hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
  ggtitle("Posterior of H")
rm(list='df')

x11(height = 3, width = 3); plot_postH

# Plotting the average density over iterations and compare with true curve
# Computing estimated density
data_range <- sapply(data, range)
est_dens <- ComputeDensities(chains, 500, data_range)

# Computing true density
x_grid <- seq(data_range[1,1], data_range[2,1], length.out = 500)
xgrid_expand <- t(rbind(replicate(ncomponents, x_grid, simplify = "matrix")))
true_dens <- t(as.matrix(weights)) %*% dnorm(xgrid_expand,means,std_devs)

# Density comparison - visual plot
y <- rbind(est_dens[[1]], true_dens); row.names(y) <- c("Estimated", "True")
plot_densCompare <- ggplot(data = DataRecast(x_grid,y), aes(x=Grid, y=Value, group=Density, col=Density)) +
  geom_line(lwd=1) + theme(plot.title = element_text(face="bold", hjust = 0.5)) +
  ggtitle("Area 1")

x11(width = 5.8, height = 3); plot_densCompare

###########################################################################

###########################################################################
# Sensitivity w.r.t. number of observations -------------------------------

# Clean console and set number of observations
rm(list=ls()); set.seed(230196)
Ns <- c(100,250,500,1000)
W <- matrix(0, nrow = 1, ncol = 1, byrow = T)
plots <- list()

# Sampler run loop
for (i in 1:length(Ns)) {
  # Generate data (1 location, mixture of 3 normals)
  ngroups <- 1; ncomponents <- 3; N <- Ns[i]
  means <- c(-4,1,5); std_devs <- c(1,1,1); weights <- c(1/6,3/6,2/6)
  cluster_alloc <- sample(1:ncomponents, prob = weights, size = N, replace = T)
  data <- list(); data[[1]] <- rnorm(N, mean = means[cluster_alloc], sd = std_devs[cluster_alloc])

  # Grab input filenames
  params_filename = system.file("input_files/rjsampler_params.asciipb", package = "SPMIX")

  # Run Spatial sampler
  out <- SPMIXSampler(0, 10000, 1, data, W, params_filename, type = "rjmcmc")

  # Deserialization
  chains <- sapply(out, function(x) DeserializeSPMIXProto("UnivariateState",x))
  H_chain <- sapply(chains, function(x) x$num_components)

  # Traceplot for the whole chain (no burnin, no thinning)
  df <- data.frame("Iteration"=1:10000, "LowPoints"=H_chain-0.3, "UpPoints"=H_chain+0.3)
  tmp <- ggplot(data=df, aes(x=Iteration, y=LowPoints, xend=Iteration, yend=UpPoints)) +
    ylim(range(df[,-1])) + ylab("NumComponents") + geom_segment(lwd=0.1) +
    theme(plot.title = element_text(face="bold", hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
    ggtitle(paste0("Traceplot of H - ", Ns[i], " obs."))
  rm(list='df')

  # Save plot
  plots[[i]] <- tmp; rm(list='tmp')
}

# Visualization
x11(height = 6, width = 6); gridExtra::grid.arrange(grobs=plots, nrow=2, ncol=2)

###########################################################################

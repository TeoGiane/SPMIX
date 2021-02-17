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
# Data Generation - Scenario I --------------------------------------------

# Generating data
set.seed(230196)
I <- 6; Ns <- rep(100,I); means <- c(rep(-5,2),rep(0,2),rep(5,2))
data <- list(); labels <- sprintf("g%d", 1:I)
for (i in 1:I) {
  data[[i]] <- rnorm(Ns[i], means[i], 1)
}
rm(list = 'i')
names(data) <- labels

# Setting initial W
W <- matrix(1,I,I)

###########################################################################

###########################################################################
# Data Generation - Scenario II -------------------------------------------

# Generating data
set.seed(230196)
I <- 6; Ns <- rep(100,I)
data <- list(); labels <- sprintf("g%d", 1:I)
for (i in 1:I) {
  if (i %in% c(1,2)){
    data[[i]] <- metRology::rt.scaled(n=Ns[i], df=6, mean=-4, sd=1)
  }
  if (i %in% c(3,4)){
    data[[i]] <- sn::rsn(n=Ns[i], xi=4, omega=4, alpha=1); attributes(data[[i]]) <- NULL
  }
  if (i %in% c(5,6)) {
    data[[i]] <- rchisq(n=Ns[i], df=3)
  }
}
rm(list = 'i')
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
out <- Sampler.BoundaryDetection(burnin, niter, thin, data, W, params_filename)

###########################################################################

###########################################################################
# Posterior Analysis - Scenario I -----------------------------------------

# Deserialization
chains <- sapply(out, function(x) DeserializeSPMIXProto("UnivariateState",x))
G_chain <- lapply(chains, function(x) matrix(x$G$data,x$G$rows,x$G$cols))

# Compute plinks and median graph according to mean and estimated graph
plinks <- Reduce('+', G_chain)/length(G_chain)
G_est <- ifelse(plinks > 0.5, 1, 0)

# Computing estimated densities
data_ranges <- sapply(data, range)
estimated_densities <- ComputeDensities(chains, 500, data_ranges, labels)

# Computing true densities for comparison
true_densities <- list()
for (i in 1:I) {
  x <- seq(data_ranges[1,i], data_ranges[2,i], length.out=500)
  true_densities[[i]] <- dnorm(x, means[i], 1)
}
rm(list = c('x','i'))
names(true_densities) <- labels

###########################################################################

###########################################################################
# Posterior Analysis - Scenario II ----------------------------------------

# Deserialization
chains <- sapply(out, function(x) DeserializeSPMIXProto("UnivariateState",x))
G_chain <- lapply(chains, function(x) matrix(x$G$data,x$G$rows,x$G$cols))

# Compute plinks and median graph according to mean and estimated graph
plinks <- Reduce('+', G_chain)/length(G_chain)
G_est <- ifelse(plinks > 0.5, 1, 0)

# Computing estimated densities
data_ranges <- sapply(data, range)
estimated_densities <- ComputeDensities(chains, 500, data_ranges, labels)

# Computing true densities for comparison
true_densities <- list()
for (i in 1:I) {
  x <- seq(data_ranges[1,i], data_ranges[2,i], length.out=500)
  if (i %in% c(1,2)){
    true_densities[[i]] <- metRology::dt.scaled(x, df=6, mean=-4, sd=1)
  }
  if (i %in% c(3,4)){
    true_densities[[i]] <- sn::dsn(x, xi=4, omega=4, alpha=1); attributes(true_densities[[i]]) <- NULL
  }
  if (i %in% c(5,6)) {
    true_densities[[i]] <- dchisq(x, df=3)
  }
}
rm(list = c('x','i'))
names(true_densities) <- labels

###########################################################################

###########################################################################
# Visualization -----------------------------------------------------------

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

# Visualization
x11(height = 6, width = 8.27); gridExtra::grid.arrange(grobs=plots_area, nrow=2, ncol=3)
x11(height = 3, width = 4.5); fields::image.plot(plinks) # Rifare su illustrator

# Save pdf
dev.copy2pdf(device=x11, file="BD_fixed_p05_plinks.pdf"); dev.off()
dev.copy2pdf(device=x11, file="BD_fixed_p05_Densities.pdf"); dev.off()

###########################################################################

# Drafts
AUC <- function(plinks, W_true, thresholds = seq(min(plinks[upper.tri(plinks)]),max(plinks[upper.tri(plinks)]),by = 0.01)) {
  aucs <- vector()
  for (i in 1:length(thresholds)) {
    Wtmp <- ifelse(plinks > thresholds[i], 1, 0)
    aucs[i] <- suppressMessages(as.numeric(pROC::auc(pROC::roc(as.vector(W_true), as.vector(Wtmp)))))
  }
  result <- data.frame("Threshold" = thresholds, "AUC"=aucs)
  return(result)
}

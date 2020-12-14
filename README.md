# SPMIX - Spatial Mixture Models in R
<strong>Author</strong>: Matteo Gianella <br>
<strong>Supervisors</strong>: Prof. Alessandra Guglielmi, Dott. Mario Beraha

This repo collects the developing of the thesis project on Spatial Mixture models.

## Requirements
This package relies, for various purposes, on the following libraries:
<par>
* [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page): it provides advanced and optimized linear algebra utilities, used extensively throughout the code.
* [GSL](https://www.gnu.org/software/gsl/): the GNU Scientific Library is used as linear algebra tool in the <code>/src/polyagamma</code> subfolder, which is a third-party code.
* [Stan Math](https://mc-stan.org/math/): this library provides both probability distributions utilities and the reverse mode automatic differentiation module, required in the samplers. Using this library, also its dependencies are necessairly required.
* [Google Protocol Buffer](https://developers.google.com/protocol-buffers): this library provides tools for quick serialization of structured data.

In order for the code to be available both for Windows and Unix systems, it attempts to provide headers and libraries through already existing R packages. Since this is not possible for all the aforementioned dependencies, some extra effort is required to make them all available on your operating system.

## Required R Packages
The DESCRIPTION file lists all the dependencies of this package. Opening an R terminal, these can be easily installed via
```r
install.packages(c("BH", "Rcpp", "RcppEigen", "RcppParallel", "StanHeaders", "RProtoBuf"))
```
Since this is a developing package, be sure that the development tools for R are installed. If only R is installed, in Linux you can quickly install them via terminal typing:
```shell
sudo apt install r-base-dev
```
## Required Libraries
At this stage, the package requires GSL and Google Protocol Buffers Libraries.

## Installation for Unix Systems
First of all, install the <code>devtools</code> package in R. To install the package, clone this repo, open an R terminal, set the repo directory as working directory via <code>setwd("path/to/the/folder")</code> and then type:
```r
devtools::install()
```
As an alternative, you can download this repo as a tar.gz archive, open a Terminal in the folder where the archive is stored and then type:
```shell
R CMD INSTALL SPMIX-master.tar.gz
```


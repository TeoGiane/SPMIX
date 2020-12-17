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
* [Google Protocol Buffers](https://developers.google.com/protocol-buffers): this library provides tools for quick serialization of structured data.

In order for this package to be available both for Windows and Unix systems, headers and libraries are provided, whenever possible, through already existing R packages. Nevertheless, some aforementioned dependencies require some extra effort to make them available on your operating system (expecially for Windows users).

## Installation - Unix
### R development tools
First of all, <code>SPMIX</code> is an R development package, hence you need to have R installed on your computer as well as all development tools. The [R homepage](https://www.r-project.org/) offers extensive documentation on how to install R on your machine. Moreover, it is also advisable to have an IDE for R installed, in order to simplify future workflow. [RStudio](https://rstudio.com/) is a pretty complete programme and it is highly advisable to have it installed on your system.

Since the development tools package for R has many dependencies that rely on external libraries, it is advisable to install them at this stage in order to avoid several iterations to get <code>devtools</code> installed. Thus, you can easily install those libraries via package manager with
```shell
sudo apt get install libssl-dev libcurl4-openssl-dev libxml2-dev libgit2-dev libnode-dev
```
At the end of this procedure, you will be able to install the <code>devtools</code> package and its dependencies from an R console simply typing
```r
install.packages("devtools")
```
### External Libraries Dependencies
In this package, the <code>gsl</code> and the <code>protobuf</code> libraries are linked to the package as external libraries. As far as the GNU Scientific Library is concerned, it can be easily installed, in most Unix operating systems, through the package manager with commands like
```shell
sudo apt-get install libgsl-dev
```
The protobuf library needs to be installed from source following these [instructions](https://github.com/protocolbuffers/protobuf/blob/master/src/README.md).
Be aware that this library is under continuous development, so you may incurr into some troubles during installation, expecially at checking stage. The following notes may be helpful in order to avoid issues while installing protobuf.

**Note on Compiler**

If you are using GNU compiler version 10.2 (available, for instance, in Ubuntu 20.10), you will have issues at compile time if you download the release of protobuf. So it is advisable to use more stable versions of the compiler (no problems have been reported with g++ version 9.3.0)

**Note on Memory**

The latest release of protobuf, one of the tests requires a huge amount of memory so you may see your <code>make check</code> failed due to the fact that your machine or VM has not enough memory. With 8GB of RAM available, there should not be issues. In case you have less amount of memory, make sure to provide an adeguate amount of swap memory to you Unix machine (the tests were succesfull on a VM with 4GB of RAM and 8GB of swap, for instance). As an alternative, you can pick an older release, in which the aforemontioned test is missing (e.g. protobuf v. 3.13).

### R Packages Dependencies
The DESCRIPTION file lists all the dependencies of this package. Both the <strong>Eigen</strong> and the <strong>Stan Math</strong> libraries are available as R packages, which are, respectively, <code>RcppEigen</code> and <code>StanHeaders</code>, which itself depend on other libraries that will be installed automatically. Moreover, since this package manages compiled code through the <code>[Rcpp](http://www.rcpp.org/)</code> package, this should be installed as well. On the other hand, <code>RProtoBuf</code> is a required package due to the fact that <code>SPMIX</code> relies on Google Protocol Buffers as serialization tool and, hence, an easy-to-use R interface to this API is suggested.

The installation of all these packages is trivial, since you only need a single R command to do it.
```r
install.packages(c("BH", "Rcpp", "RcppEigen", "RcppParallel", "StanHeaders", "RProtoBuf"))
```

### SPMIX
Once all dependencies have been installed, the <code>SPMIX</code> is rather simple to install. To do so, clone this repository via
```shell
git clone https://github.com/TeoGiane/SPMIX
```
Then, open an R terminal and set as current working directory the SPMIX directory and then install via
```r
setwd("path/to/the/folder")
devtools::install()
```
Once installed, you can include the <code>SPMIX</code> package in your workflow as a standard CRAN package with <code>library("SPMIX")</code>.

## Installation - Windows
### R development tools
### External Libraries Dependencies
### R Packages Dependencies
### SPMIX

<!---
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
-->

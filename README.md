# SPMIX - Spatial Mixture Models in R
<strong>Author</strong>: Matteo Gianella <br>
<strong>Supervisors</strong>: Prof. Alessandra Guglielmi, Dott. Mario Beraha

## Table of Contents
1. [Overview](#overview)
2. [Requirements](#requirements)
3. [Installation - Unix](#installation---unix)
4. [Installation - Windows](#installation---windows)

## Overview
SPMIX collects sampling schemes that performs density estimation for spatially dependent areal data, in a Bayesian Non-Parametric setting. Data on each area are modelled as a finite mixture of Gaussian kernels and the weights provides the spatial dependence among neighbours via the logistic multivariate CAR prior. The number of components of the mixture can be either fixed or variable: in the second case, a prior on such quantity is added and a reversible jump scheme is adopted to estimate the posterior distribution of the model.

## Requirements
This package relies, for various purposes, on the following libraries:

* [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page): it provides advanced and optimized linear algebra utilities, used extensively throughout the code.
* [GSL](https://www.gnu.org/software/gsl/): the GNU Scientific Library is used as linear algebra tool in the <code>/src/polyagamma</code> subfolder, which is a third-party code used in the samplers.
* [Stan Math](https://mc-stan.org/math/): this library provides both probability distributions utilities and the reverse mode automatic differentiation module, required in the samplers. Using this library, also its dependencies are necessairly required.
* [Google Protocol Buffers](https://developers.google.com/protocol-buffers): this library provides tools for quick serialization of structured data. This protocol is used as serialization tool for samplers inputs and outputs.

In order for this package to be available both for Windows and Unix systems, headers and libraries are provided, whenever possible, through already existing R packages. Nevertheless, some aforementioned dependencies require some extra effort to make them available on your operating system (expecially for Windows users).

## Installation - Unix
### R development tools
<code>SPMIX</code> is an R development package, hence you need to have R installed on your computer as well as all development tools. The [R homepage](https://www.r-project.org/) offers extensive documentation on how to install R on your machine. Moreover, it is also advisable to have an IDE for R installed, in order to simplify future workflow. As a suggestion, [RStudio](https://rstudio.com/) is a pretty complete programme and it is highly advisable to have it installed on your system.

Since the development tools package for R has many dependencies that rely on external libraries, it is advisable to install them at this stage in order to avoid several iterations to install the <code>devtools</code> package. To install those libraries, type on terminal
```shell
sudo apt-get install libssl-dev libcurl4-openssl-dev libxml2-dev libgit2-dev libnode-dev
```
**Remark**: The command above is valid for Ubuntu/Debian users. On other platforms, please use the corresponding package managing tool to install them.
At the end of this procedure, you will be able to install the <code>devtools</code> package and its dependencies from an R console simply typing
```r
install.packages("devtools")
```
### External Libraries Dependencies
In this package, the <code>gsl</code> and the <code>protobuf</code> libraries are linked to the package as external libraries. As far as the GNU Scientific Library is concerned, it can be easily installed, in most Unix operating systems, through the package manager with commands like
```shell
sudo apt-get install libgsl-dev
```
The protobuf library needs to be installed from source following the [protobuf C++ Installation Instructions](https://github.com/protocolbuffers/protobuf/blob/master/src/README.md).
Be aware that this library is under continuous development, so you may incurr into some troubles during installation, expecially at checking stage. The following notes may be helpful in order to avoid issues while installing protobuf.

**Note on Compiler**: If you are using GNU compiler version 10.2 (available, for instance, in Ubuntu 20.10), you will have issues at compile time if you download the release of protobuf. So it is advisable to use more stable versions of the compiler (no problems have been reported with g++ version 9.3.0)

**Note on Memory**: The latest release of protobuf, one of the tests requires a huge amount of memory, so you may see your <code>make check</code> failed due to the fact that your machine or VM has not enough memory. With 8GB of RAM available, there should not be issues. In case you have a lower amount of memory, make sure to provide an adeguate amount of swap memory to you Unix machine (the tests were succesfull on a VM with 4GB of RAM and 8GB of swap, for instance). As an alternative, you can pick an older release, in which the aforemontioned test is missing (e.g. protobuf v. 3.13).

### R Packages Dependencies
The DESCRIPTION file lists all the dependencies of this package. Both the <strong>Eigen</strong> and the <strong>Stan Math</strong> headers and libraries are available as R packages, which are, respectively, <code>RcppEigen</code> and <code>StanHeaders</code>, which themselves depend on other libraries (e.g. <code>RcppParallel</code> or <code>rstan</code>) that will be installed automatically. <code>RcppProgress</code>, instead, offers display classes to print samplers' progresses during execution. LaTeX support for documentation in the R helper is offered by <code>mathjaxr</code>. Finally, since this package manages compiled code through the <code>Rcpp</code> package, this should be installed as well. On the other hand, <code>RProtoBuf</code> is a required package due to the fact that <code>SPMIX</code> relies on Google Protocol Buffers as serialization tool and, hence, an easy-to-use R interface to this API is suggested.

The installation of all these packages is trivial, since you only need a single R command to do it.
```r
install.packages(c("BH", "Rcpp", "RcppEigen", "RcppParallel", "RcppProgress",
                   "RProtoBuf", "StanHeaders", "mathjaxr", "rstan"))
```

### SPMIX
Once all dependencies have been installed, the <code>SPMIX</code> is extremely simple to install. To do so, open R and simply type
```r
devtools::install_github("TeoGiane/SPMIX")
```
This command will automatically download, build and install <code>SPMIX</code> in your package library. Once installed, you can import it in your workflow as a standard package with <code>library("SPMIX")</code>.

## Installation - Windows
### R development tools
<code>SPMIX</code> is an R development package, hence you need to have R installed on your computer as well as all development tools. The [R homepage](https://www.r-project.org/) offers extensive documentation on how to install R on your machine. Moreover, it is also advisable to have an IDE for R installed, in order to simplify future workflow. As a suggestion, [RStudio](https://rstudio.com/) is a pretty complete programme and it is highly advisable to have it installed on your system. Finally, in order to install the R development tools package in Windows systems, all the required compilers needs to be available on your machine. In windows, these are provided in a toolchain bundle called [Rtools](https://cran.r-project.org/bin/windows/Rtools/). At the previous link you can find all the instructions to set up Rtools on your system.

<!---
In order to add the Rtools binary path to PATH variable and set compatibility flags for the compiler to avoid unnecessary warnings at compile time, type in the Rtools Bash (available when Rtools is installed) the following commands:
```shell
echo 'PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"' >> Documents/.Renviron
mkdir -p Documents/.R
echo 'CXX17FLAGS = -O2 -mtune=native -mmmx -msse -msse2 -msse3 -mssse3 -msse4.1 -msse4.2' >> Documents/.R/Makevars.win
```
-->

At the end of this procedure, you will be able to install the <code>devtools</code> package and its dependencies from an R console simply typing
```r
install.packages("devtools")
```

### External Libraries Dependencies
In this package, the <code>gsl</code> and the <code>protobuf</code> libraries are linked to the package as external libraries. While Unix systems are organized in repositories and manages libraries via package managing programmes, on Windows things are more complicated, in general. The general approach adopted here is to exploit the R toolchain provided with Rtools to compile and install the external libraries dependencies emulating a Unix environment. Main reference for this procedure is the [Rtools Packages Repository](https://github.com/r-windows/rtools-packages), which offers an extremely automatized procedure to build and install the required libraries for both 32 and 64 bits architectures. In the following we report the necessary instructions, with the aim of furtherly simplifying the procedure.

Open the Rtools Bash terminal (it has been installed with Rtools, you can find it simply searching for "Rtools Bash" on Windows search bar) and simply follow these instructions:

* **Enable all upstream MSYS2 build tools**: this will allow the installation of further packages required by <code>protobuf</code>. Type:
```shell
sed -i '78,79 s/^#//' $RTOOLS40_HOME\\etc\\pacman.conf
echo 'SigLevel = Never' >> $RTOOLS40_HOME\\etc\\pacman.conf
```
* **Clone the Rtools Packages repo**: after git installation (you can skip the command on the first line in case you have [Git](https://gitforwindows.org/) already installed on your OS), clone the <code>rtools-packages</code> repository. You can choose whichever location you like to clone the repo. Use the following commands:
```shell
pacman -S git
git clone https://github.com/r-windows/rtools-packages.git
cd rtools-packages
```
* **Build and Install the gsl library**: the following commands also clean the folder from files created during the build. Type:
```shell
cd mingw-w64-gsl
makepkg-mingw --syncdeps --noconfirm
pacman -U mingw-w64-i686-gsl-2.6-1-any.pkg.tar.xz
pacman -U mingw-w64-x86_64-gsl-2.6-1-any.pkg.tar.xz
rm -f -r pkg src *.xz *.gz
cd ../
```
* **Build and Install the protobuf library**: the following commands also clean the folder from files created during the build. Type:
```shell
cd mingw-w64-protobuf
makepkg-mingw --syncdeps --noconfirm
pacman -U mingw-w64-i686-protobuf-3.12.4-1-any.pkg.tar.xz
pacman -U mingw-w64-x86_64-protobuf-3.12.4-1-any.pkg.tar.xz
rm -f -r pkg src *.xz *.gz
cd ../
```

### R Packages Dependencies
The DESCRIPTION file lists all the dependencies of this package. Both the <strong>Eigen</strong> and the <strong>Stan Math</strong> headers and libraries are available as R packages, which are, respectively, <code>RcppEigen</code> and <code>StanHeaders</code>, which themselves depend on other libraries (e.g. <code>RcppParallel</code> or <code>rstan</code>) that will be installed automatically. <code>RcppProgress</code>, instead, offers display classes to print samplers' progresses during execution. LaTeX support for documentation in the R helper is offered by <code>mathjaxr</code>. Finally, since this package manages compiled code through the <code>Rcpp</code> package, this should be installed as well. On the other hand, <code>RProtoBuf</code> is a required package due to the fact that <code>SPMIX</code> relies on Google Protocol Buffers as serialization tool and, hence, an easy-to-use R interface to this API is suggested.

The installation of all these packages is trivial, since you only need a single R command to do it.
```r
install.packages(c("BH", "Rcpp", "RcppEigen", "RcppParallel", "RcppProgress",
                   "RProtoBuf", "StanHeaders", "mathjaxr", "rstan"))
```

### SPMIX
Once all dependencies have been installed, the <code>SPMIX</code> is extremely simple to install. To do so, open R and simply type
```r
devtools::install_github("TeoGiane/SPMIX")
```
This command will automatically download, build and install <code>SPMIX</code> in your package library. Once installed, you can import it in your workflow as a standard package with <code>library("SPMIX")</code>.

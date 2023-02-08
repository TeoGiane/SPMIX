# SPMIX - Spatial Mixture Models in R

<code>SPMIX</code> collects sampling schemes that performs density estimation for spatially dependent areal data, in a Bayesian Nonparametric setting. Data on each area are modelled as a finite mixture of Gaussian kernels and the weights provides the spatial dependence among neighbours via the logistic multivariate CAR prior. The number of components of the mixture can be either fixed or variable: in the second case, a prior on such quantity is added and a reversible jump scheme is adopted to estimate the posterior distribution of the model.

## Dependencies
<code>SPMIX</code> depends on:

- [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page): A C++ template library for linear algebra.
- [GNU Scientific Library](https://www.gnu.org/software/gsl/): A numerical library for C and C++ programmers. It is required in <code>src/polyagamma/</code> subfolder.
- [Stan Math Library](https://mc-stan.org/math/): A BSD-3 licensed C++, reverse-mode automatic differentiation library designed to facilitate the construction and utilization of algorithms that utilize derivatives.
- [Google Protocol Buffers](https://developers.google.com/protocol-buffers): A language-neutral, platform-neutral extensible mechanism for serializing structured data.

## Prerequisites
<code>SPMIX</code> is an <code>R</code> package, so make sure to have <code>R</code> installed on your machine (see the [project homepage](https://www.r-project.org/) otherwise) as well as the 
<code>devtools</code> package. Once <code>R</code> is installed, it can be insalled with the following command called from an <code>R</code> terminal.
```R
> install.packages("devtools")
```

**Note for Linux users**: <code>devtools</code> and all <code>R</code> packages in Linux are built from source. Thus, you need to install some external libraries via package manager. Usually, the following command on Ubuntu/Debian leads to a succesfull installation of <code>devtools</code>.

```shell
# apt install libssl-dev libcurl4-openssl-dev libxml2-dev libgit2-dev libnode-dev
```

<br>

# Installation

This package links <code>gsl</code> and <code>protobuf</code> as external libraries, while <code>Eigen</code> and <code>stan/math</code> are provided as <code>R</code> packages. All <code>R</code> dependencies are automatically installed via <code>devtools</code>.

External libraries are handled differently according to your operating system. Follow the instruction for [Linux](#linux-systems) or [Windows](#windows-systems) machines accordingly.

<br>

## Linux Systems

1. Install <code>gsl</code> via package manager.

```shell
# apt install libgsl-dev
```

2. Install <code>protobuf</code> from source using [CMake](https://cmake.org/). The latter can be installed via package manager.

```shell
# apt install cmake
```

3. Download and extract the latest Protocol Buffers Release from [here](https://github.com/protocolbuffers/protobuf/releases/latest/). For a minimal installation, download the <code>protobuf-cpp-<VERSION\>.\*</code> asset.

4. Once inside the folder containing <code>protobuf</code> source files, build and install it using the following command.

```shell
$ cmake . -Dprotobuf_BUILD_SHARED_LIBS=ON
$ cmake --build . [--parallel <NUM_THREADS>]
$ ctest --verbose
# cmake --install
# ldconfig
```
5. Install <code>SPMIX</code> and all its package dependencies.

```R
> devtools::install_github("TeoGiane/SPMIX")
```

<br>

## Windows Systems

On Windows machines, external libraries will be available during installation automatically thanks to the [rwinlib](https://github.com/rwinlib) project. Thus, you only need to install <code>SPMIX</code> via:

```R
> devtools::install("TeoGiane/SPMIX")
```

<br>

# License
This package is licensed under the [MIT License](https://github.com/TeoGiane/SPMIX/blob/master/LICENSE.md).

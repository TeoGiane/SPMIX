# Load RProtoBuf automatically when this package is loaded
.onLoad <- function (libname, pkgname) {
  require(RProtoBuf)
}


# Unloading the dynamic library when SPMIX is unloaded
.onUnload <- function (libpath) {
  library.dynam.unload("SPMIX", libpath)
}

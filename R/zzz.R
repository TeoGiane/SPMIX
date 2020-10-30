# Unloading the dynamic library when SPMIX is unloaded
.onUnload <- function (libpath) {
  library.dynam.unload("SPMIX", libpath)
}

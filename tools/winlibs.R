# gsl library to build against it with RTools
GSL_VERSION <- commandArgs(TRUE)
if (!file.exists(sprintf("../windows/gsl-%s/include/gsl/gsl_version.h", GSL_VERSION))) {
  if (getRversion() < "3.3.0") setInternet2()
  download.file(sprintf("https://github.com/rwinlib/gsl/archive/v%s.zip", GSL_VERSION), "lib.zip", quiet = TRUE)
  dir.create("../windows", showWarnings = FALSE)
  unzip("lib.zip", exdir = "../windows")
  unlink("lib.zip")
}

# protobuf library to build against it with RTools
PROTOBUF_VERSION <- commandArgs(TRUE)
if (!file.exists(sprintf("../windows/protobuf-%s/include/google/protobuf/descriptor.h", PROTOBUF_VERSION))) {
  if (getRversion() < "3.3.0") setInternet2()
  download.file(sprintf("https://github.com/rwinlib/protobuf/archive/v%s.zip", PROTOBUF_VERSION), "lib.zip", quiet = TRUE)
  unzip("lib.zip", exdir = "../windows")
  unlink("lib.zip")
}

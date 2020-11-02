// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// alr
Eigen::VectorXd alr(Eigen::VectorXd x, bool pad_zero);
RcppExport SEXP _SPMIX_alr(SEXP xSEXP, SEXP pad_zeroSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type x(xSEXP);
    Rcpp::traits::input_parameter< bool >::type pad_zero(pad_zeroSEXP);
    rcpp_result_gen = Rcpp::wrap(alr(x, pad_zero));
    return rcpp_result_gen;
END_RCPP
}
// inv_alr
Eigen::VectorXd inv_alr(Eigen::VectorXd x, bool padded_zero);
RcppExport SEXP _SPMIX_inv_alr(SEXP xSEXP, SEXP padded_zeroSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type x(xSEXP);
    Rcpp::traits::input_parameter< bool >::type padded_zero(padded_zeroSEXP);
    rcpp_result_gen = Rcpp::wrap(inv_alr(x, padded_zero));
    return rcpp_result_gen;
END_RCPP
}
// runSpatialSampler
std::vector<Rcpp::RawVector> runSpatialSampler(int burnin, int niter, int thin, const std::vector<std::vector<double>>& data, const Eigen::MatrixXd& W, Rcpp::S4 params, const std::vector<Eigen::MatrixXd>& covariates, bool display_progress);
RcppExport SEXP _SPMIX_runSpatialSampler(SEXP burninSEXP, SEXP niterSEXP, SEXP thinSEXP, SEXP dataSEXP, SEXP WSEXP, SEXP paramsSEXP, SEXP covariatesSEXP, SEXP display_progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< const std::vector<std::vector<double>>& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type W(WSEXP);
    Rcpp::traits::input_parameter< Rcpp::S4 >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< const std::vector<Eigen::MatrixXd>& >::type covariates(covariatesSEXP);
    Rcpp::traits::input_parameter< bool >::type display_progress(display_progressSEXP);
    rcpp_result_gen = Rcpp::wrap(runSpatialSampler(burnin, niter, thin, data, W, params, covariates, display_progress));
    return rcpp_result_gen;
END_RCPP
}
// readMatrixFromCSV
Eigen::MatrixXd readMatrixFromCSV(std::string filename);
RcppExport SEXP _SPMIX_readMatrixFromCSV(SEXP filenameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type filename(filenameSEXP);
    rcpp_result_gen = Rcpp::wrap(readMatrixFromCSV(filename));
    return rcpp_result_gen;
END_RCPP
}
// readDataFromCSV
std::vector<std::vector<double>> readDataFromCSV(std::string filename);
RcppExport SEXP _SPMIX_readDataFromCSV(SEXP filenameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type filename(filenameSEXP);
    rcpp_result_gen = Rcpp::wrap(readDataFromCSV(filename));
    return rcpp_result_gen;
END_RCPP
}
// spmixLogLikelihood
double spmixLogLikelihood(const Rcpp::S4& state, const std::vector<std::vector<double>>& data, const Eigen::MatrixXd& W, const Rcpp::S4& params);
RcppExport SEXP _SPMIX_spmixLogLikelihood(SEXP stateSEXP, SEXP dataSEXP, SEXP WSEXP, SEXP paramsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type state(stateSEXP);
    Rcpp::traits::input_parameter< const std::vector<std::vector<double>>& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type W(WSEXP);
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type params(paramsSEXP);
    rcpp_result_gen = Rcpp::wrap(spmixLogLikelihood(state, data, W, params));
    return rcpp_result_gen;
END_RCPP
}
// stan_HelloWorld
void stan_HelloWorld();
RcppExport SEXP _SPMIX_stan_HelloWorld() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    stan_HelloWorld();
    return R_NilValue;
END_RCPP
}
// fromProto_tostring
std::vector<Rcpp::RawVector> fromProto_tostring();
RcppExport SEXP _SPMIX_fromProto_tostring() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(fromProto_tostring());
    return rcpp_result_gen;
END_RCPP
}
// readingStates
void readingStates(std::vector<Rcpp::RawVector> raw_vect);
RcppExport SEXP _SPMIX_readingStates(SEXP raw_vectSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<Rcpp::RawVector> >::type raw_vect(raw_vectSEXP);
    readingStates(raw_vect);
    return R_NilValue;
END_RCPP
}
// messageFromR
void messageFromR(Rcpp::S4 params);
RcppExport SEXP _SPMIX_messageFromR(SEXP paramsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::S4 >::type params(paramsSEXP);
    messageFromR(params);
    return R_NilValue;
END_RCPP
}
// newton_opt_test
void newton_opt_test(const Rcpp::S4& state, const std::vector<std::vector<double>>& data, const Eigen::MatrixXd& W, const Rcpp::S4& params, const Rcpp::S4& options);
RcppExport SEXP _SPMIX_newton_opt_test(SEXP stateSEXP, SEXP dataSEXP, SEXP WSEXP, SEXP paramsSEXP, SEXP optionsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type state(stateSEXP);
    Rcpp::traits::input_parameter< const std::vector<std::vector<double>>& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type W(WSEXP);
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type options(optionsSEXP);
    newton_opt_test(state, data, W, params, options);
    return R_NilValue;
END_RCPP
}
// grad_ascent_test
void grad_ascent_test(const Rcpp::S4& options);
RcppExport SEXP _SPMIX_grad_ascent_test(SEXP optionsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type options(optionsSEXP);
    grad_ascent_test(options);
    return R_NilValue;
END_RCPP
}
// RJsampler_test
void RJsampler_test(const std::vector<std::vector<double>>& data, const Eigen::MatrixXd& W, const Rcpp::S4& params, const Rcpp::S4& options);
RcppExport SEXP _SPMIX_RJsampler_test(SEXP dataSEXP, SEXP WSEXP, SEXP paramsSEXP, SEXP optionsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<std::vector<double>>& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type W(WSEXP);
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::S4& >::type options(optionsSEXP);
    RJsampler_test(data, W, params, options);
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SPMIX_alr", (DL_FUNC) &_SPMIX_alr, 2},
    {"_SPMIX_inv_alr", (DL_FUNC) &_SPMIX_inv_alr, 2},
    {"_SPMIX_runSpatialSampler", (DL_FUNC) &_SPMIX_runSpatialSampler, 8},
    {"_SPMIX_readMatrixFromCSV", (DL_FUNC) &_SPMIX_readMatrixFromCSV, 1},
    {"_SPMIX_readDataFromCSV", (DL_FUNC) &_SPMIX_readDataFromCSV, 1},
    {"_SPMIX_spmixLogLikelihood", (DL_FUNC) &_SPMIX_spmixLogLikelihood, 4},
    {"_SPMIX_stan_HelloWorld", (DL_FUNC) &_SPMIX_stan_HelloWorld, 0},
    {"_SPMIX_fromProto_tostring", (DL_FUNC) &_SPMIX_fromProto_tostring, 0},
    {"_SPMIX_readingStates", (DL_FUNC) &_SPMIX_readingStates, 1},
    {"_SPMIX_messageFromR", (DL_FUNC) &_SPMIX_messageFromR, 1},
    {"_SPMIX_newton_opt_test", (DL_FUNC) &_SPMIX_newton_opt_test, 5},
    {"_SPMIX_grad_ascent_test", (DL_FUNC) &_SPMIX_grad_ascent_test, 1},
    {"_SPMIX_RJsampler_test", (DL_FUNC) &_SPMIX_RJsampler_test, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_SPMIX(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

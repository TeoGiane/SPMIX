/* This source file contains all utilities that are exported to R through Rcpp. Any other function 
 * should be put here in order to preserve coherence in the code
 */

#ifndef SPMIX_EXPORTS
#define SPMIX_EXPORTS


#include <deque>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include <Rcpp.h>

#include "utils.h"
/* UNCOMMENT TO SE THE COMPILE ERROR */
/*#include "sampler.h"
#include "univariate_mixture_state.pb.h"*/

// using return_t = std::tuple<std::vector<std::string>, double>;
using return_t = std::vector<std::string>; // Simplified to check if Rcpp is able to manage the output.


/*TODO: finire documentazione*/

//' Additive log ratio
//'
//' This utility computes the additive log ratio transorm of a given vector. This transformation is
//' defined as: (METTI DEF)
//' @param x Vector of double in \ifelse{html}{\out{R<sup>H</sup>}}{\eqn{R^{H}}}.
//' @param pad_zero Boolean, default to FALSE.
//' @return Vector of double in \ifelse{html}{\out{R<sup>H-1</sup>}}{\eqn{R^{H-1}}} i.e. alr(x).
//' @export
// [[Rcpp::export]]
Eigen::VectorXd alr(Eigen::VectorXd x, bool pad_zero = false) {
	return utils::Alr(x, pad_zero);
}

/*TODO: finire documentazione*/

//' Inverse additive log ratio
//'
//' This utility computes the inverse additive log ratio transorm of a given vector. This transofrmation is
//' defined as: (METTI DEF)
//' @param x Vector of double in \ifelse{html}{\out{R<sup>H-1</sup>}}{\eqn{R^{H-1}}}.
//' @param padded_zero Boolean, default to FALSE.
//' @return Vector of double in \ifelse{html}{\out{R<sup>H</sup>}}{\eqn{R^{H}}} i.e.
//' \ifelse{html}{\out{alr <sup>-1</sup>(x)}}{\eqn{alr^{-1}(x)}}.
//' @export
// [[Rcpp::export]]
Eigen::VectorXd inv_alr(Eigen::VectorXd x, bool padded_zero = false) {
  return utils::InvAlr(x, padded_zero);
}

/* FIRST ATTEMPT OF SAMPLER (no covariates, read from file) - UNCOMMENT TO SE THE COMPILE ERROR */
/*return_t runSpatialSampler(
    int burnin, int niter, int thin,
    const std::vector<std::vector<double>> &data,
    const Eigen::MatrixXd &W, const SamplerParams &params) {

    SpatialMixtureSampler spSampler(params, data, W);
    spSampler.init();

    std::vector<std::string> out;
    int log_every = 200;

    auto start = std::chrono::high_resolution_clock::now();
    for (int i=0; i < burnin; i++) {
        spSampler.sample();
        if ((i + 1) % log_every == 0)
            Rcpp::Rcout << "Burn-in, iter #" << i+1 << " / " << burnin << std::endl;
    }

    for (int i=0; i < niter; i++) {
        spSampler.sample();
        if ((i +1) % thin == 0) {
            std::string s;
            spSampler.getStateAsProto().SerializeToString(&s);
            out.push_back((py::bytes) s);
        }
        if ((i + 1) % log_every == 0)
            Rcpp::Rcout << "Running, iter #" << i+1 << " / " << niter << std::endl;
    }
    auto end = std::chrono::high_resolution_clock::now();
    double duration = std::chrono::duration<double>(end - start).count();
    
    Rcpp::Rcout << "Duration: " << duration << std::endl;
    return out;
    //return std::make_tuple(out, duration);
}



//' First try on the spatial sampler, taking input from files
//' @export
// [[Rcpp::export]]
std::vector<std::string> SPMIX_SamplerFromFiles(int burnin, int niter, int thin,
    							std::string data_file, std::string w_file, std::string params_file) {

    std::vector<std::vector<double>> data = utils::readDataFromCSV(data_file);
    Eigen::MatrixXd W = utils::readMatrixFromCSV(w_file);
    SamplerParams params = loadTextProto<SamplerParams>(params_file);
    return runSpatialSampler(burnin, niter, thin, data, W, params);
}
*/

#endif // SPMIX_EXPORTS
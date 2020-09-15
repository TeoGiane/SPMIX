#define STRICT_R_HEADERS
#include <Rcpp.h>
#include <stan/math.hpp>
#include <Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct>
#include <string>

#include "PolyaGamma.h"
#include "sampler_params.pb.h"
#include "univariate_mixture_state.pb.h"
#include "recordio.h"
#include "utils.h"

//' Simple test with stan/math C++ library
//'
//' Simply computes logN(1|2,3)
//' @export
// [[Rcpp::export]]
void stan_HelloWorld() {
	Rcpp::Rcout << "log normal(1 | 2, 3) = " << stan::math::normal_log(1, 2, 3) << std::endl;
	PolyaGamma pg;
	Rcpp::Rcout << "A PolyaGamma object has been instanciated!" << std::endl;
	Rcpp::Rcout << "Testing for KroneckerProduct..." << std::endl;
	Eigen::MatrixXd A(3,3); A << 1, 2, 3, 4, 5, 6, 7, 8, 9;
	Eigen::MatrixXd B(2,2); B << 4, 3, 2, 1;
	Rcpp::Rcout << "A:\n" << A << std::endl << std::endl;
	Rcpp::Rcout << "B:\n" << B << std::endl << std::endl;
	Eigen::MatrixXd C = Eigen::kroneckerProduct(A,B).eval();
	Rcpp::Rcout << "A x B:\n" << C << std::endl << std::endl;
  return;
}


//' Serialization testing
//' @export
// [[Rcpp::export]]
std::vector<Rcpp::RawVector> fromProto_tostring() {
	
	Rcpp::Rcout << "Reading params..." << std::endl;
	SamplerParams params;
    params = loadTextProto<SamplerParams>("/home/m_gianella/Documents/R-Packages/SPMIX_input/sampler_params.asciipb");
    Rcpp::Rcout << params.DebugString() << std::endl;

    Rcpp::Rcout << "Creating states..." << std::endl;
    std::vector<UnivariateMixtureState> states;
    UnivariateMixtureAtom* atom;
    for (int i = 0; i < 5; i++) {
        UnivariateMixtureState curr;
        curr.set_num_components(i + 3);
        for (int j = 0; j < i+3; j++) {
            atom = curr.add_atoms();
            atom->set_mean(-1.0 * i);
            atom->set_stdev(0.5 * i);
        }
        states.push_back(curr);
    }

    for (int i = 0; i < states.size(); ++i) {
    	Rcpp::Rcout << "Printing state " << i+1 << std::endl;
    	Rcpp::Rcout << states[i].DebugString() << std::endl << std::endl;
    }

    Rcpp::Rcout << "Serializing messages... ";
    std::vector<std::string> out;
    for (int i = 0; i < states.size(); ++i) {
    	std::string s;
		states[i].SerializeToString(&s);
		out.push_back(s);
    }

    Rcpp::Rcout << "done!" << std::endl;


    return utils::str2raw(out);
/*    Rcpp::Rcout << "Restoring original messages from serialized strings:" << std::endl;
	for (int i = 0; i < out.size(); ++i) {
		Rcpp::Rcout << "Message " << i+1 << " is:" << std::endl;
		UnivariateMixtureState curr;
		//std::string s(out(i));
		curr.ParseFromString(out[i]);
		Rcpp::Rcout << curr.DebugString() << std::endl << std::endl;
	}

	// Make everything RAW
	std::vector<Rcpp::RawVector> test;
	for (int i = 0; i < out.size(); ++i) {
		Rcpp::RawVector tmp(out[i].size());
		std::copy(out[i].begin(), out[i].end(), tmp.begin());
		test.push_back(tmp);
	}

    return test;*/
}


//' Translation from serialization testing
//' @export
// [[Rcpp::export]]
void readingStates(std::vector<Rcpp::RawVector> raw_vect){

	Rcpp::Rcout << "Restoring original messages from raw vectors:" << std::endl;
	std::vector<std::string> str_vect(utils::raw2str(raw_vect));
	for (int i = 0; i < str_vect.size(); ++i) {
		Rcpp::Rcout << "Message " << i+1 << " is:" << std::endl;
		UnivariateMixtureState curr;
		curr.ParseFromString(str_vect[i]);
		Rcpp::Rcout << curr.DebugString() << std::endl << std::endl;
	}
	return;
}
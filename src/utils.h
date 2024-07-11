#ifndef SRC_UTILS_HPP
#define SRC_UTILS_HPP

#include <map>
#include <random>
#include <vector>
#include <utility>
#include <sstream>
#include <string>
#include <fstream>

#include <stan/math.hpp>
#include <Eigen/Dense>
#define STRICT_R_HEADERS
#include <Rcpp.h>


namespace utils {

double trunc_normal_rng( double mu, double sigma, double lower, double upper, std::mt19937_64& rng);

double trunc_normal_lpdf(double x, double mu, double sigma, double lower, double upper);

Eigen::VectorXd Alr(Eigen::VectorXd x, bool pad_zero = false);

template<typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> InvAlr(Eigen::Matrix<T, Eigen::Dynamic, 1> x, bool padded_zero = false) {
    int D;

    if (padded_zero)
        D = x.size();
    else
        D = x.size() + 1;

    Eigen::Matrix<T, Eigen::Dynamic, 1> out(D);
    out.head(D - 1) = x.head(D-1);
    out(D - 1) = 0;
    T norm = stan::math::log_sum_exp(out);
    out = (out.array() - norm).exp();
    return out;
}

//Eigen::VectorXd InvAlr(Eigen::VectorXd x, bool padded_zero = false);

std::vector<std::vector<double>> readDataFromCSV(std::string filename);

Eigen::MatrixXd readMatrixFromCSV(std::string filename);

Eigen::VectorXd removeElem(Eigen::VectorXd vec, unsigned int toRemove);

Eigen::MatrixXd removeRow(Eigen::MatrixXd matrix, unsigned int rowToRemove);

Eigen::MatrixXd removeColumn(Eigen::MatrixXd matrix, unsigned int colToRemove);

Eigen::MatrixXd removeRowColumn(Eigen::MatrixXd matrix, unsigned int toRemove);

std::vector<int> findConnectedComponents(const Eigen::MatrixXd& adjacency);

void _dephtFirstSearch(const Eigen::MatrixXd &adjacency, int curr_node,
                       std::vector<bool> *visited,
                       std::vector<int> *node2comp,
                       int curr_comp);

double matrix_normal_prec_lpdf(Eigen::MatrixXd x, Eigen::MatrixXd m, Eigen::MatrixXd A, Eigen::MatrixXd B);

Rcpp::RawVector str2raw(const std::string & str);

std::string raw2str(const Rcpp::RawVector & raw_vect);

std::pair<double,double> range(const std::vector<std::vector<double>> & vectvect);

int min(const std::vector<double> & vect);

} // namespace utils

#endif // SRC_UTILS_HPP

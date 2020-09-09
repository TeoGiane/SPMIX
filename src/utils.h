#ifndef SRC_UTILS_HPP
#define SRC_UTILS_HPP

#include <map>
#include <random>
#include <vector>
#include <sstream>
#include <string>
#include <fstream>
#include <Eigen/Dense>
#include <stan/math/prim/mat.hpp>


namespace utils {

double trunc_normal_rng(
    double mu, double sigma, double lower, double upper,
    std::mt19937_64& rng);

double trunc_normal_lpdf(double x, double mu, double sigma, double lower, double upper);

Eigen::VectorXd Alr(Eigen::VectorXd x, bool pad_zero = false);

Eigen::VectorXd InvAlr(Eigen::VectorXd x, bool padded_zero = false);

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

double matrix_normal_prec_lpdf(
    Eigen::MatrixXd x, Eigen::MatrixXd m, Eigen::MatrixXd A,
    Eigen::MatrixXd B);

} // namespace utils

#endif // SRC_UTILS_HPP

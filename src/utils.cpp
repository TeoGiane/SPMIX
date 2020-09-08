#include "utils.hpp"
#include <iostream>

namespace utils {


// generate from truncated normal by rejection sampling
// !! might not be the best idea
double trunc_normal_rng(
        double mu, double sigma, double lower, double upper,
        std::mt19937_64& rng) {
    while (true) {
        double val = stan::math::normal_rng(mu, sigma, rng);
        if (val <= upper && val >= lower)
            return val;

    }
}

double trunc_normal_lpdf(double x, double mu, double sigma, double lower, double upper) {
    double out = stan::math::normal_lpdf(x, mu, sigma);
    out -= stan::math::log_diff_exp(
        stan::math::normal_lcdf(upper, mu, sigma),
        stan::math::normal_lcdf(lower, mu, sigma));

    return out;
}

Eigen::VectorXd Alr(Eigen::VectorXd x, bool pad_zero) {
    int D = x.size();
    Eigen::VectorXd out = x.head(D-1);
    out /= x(D-1);
    out = out.array().log();
    if (pad_zero) {
        Eigen::VectorXd out2 = Eigen::VectorXd::Zero(D);
        out2.head(D -1) = out;
        return out2;
    } else {
        return out;
    }
}

Eigen::VectorXd InvAlr(Eigen::VectorXd x, bool padded_zero) {
    int D;

    if (padded_zero)
        D = x.size();
    else
        D = x.size() + 1;

    Eigen::VectorXd out(D);
    out.head(D - 1) = x.head(D-1);
    out(D - 1) = 0;
    double norm = stan::math::log_sum_exp(out);
    out = (out.array() - norm).exp();
    return out;
}

std::vector<std::vector<double>> readDataFromCSV(std::string filename) {
    std::ifstream infile(filename);

    std::map<int, std::vector<double>> out;

    int group;
    double datum;
    std::string line;
    char delim;

    // skip header
    std::getline(infile, line);

    while (std::getline(infile, line)) {
      std::istringstream iss(line);
      if (!(iss >> group >> delim >> datum) ) { break; }
      out[group].push_back(datum);
    }

    // maps are sorted
    int ngroups = out.rbegin()->first;
    bool startsFromZero = out.begin()->first == 0;
    if (startsFromZero)
        ngroups += 1;

    std::vector<std::vector<double>> data(ngroups);
    for (int g=0; g < ngroups; g++) {
        if (startsFromZero)
            data[g] = out[g];
        else
            data[g] = out[g + 1];
    }

    return data;
}

Eigen::MatrixXd readMatrixFromCSV(std::string filename) {
    std::cout << "readMatrixFromCSV" << std::endl;
    int MAXBUFSIZE = ((int) 1e6);
    int cols = 0, rows = 0;
    double buff[MAXBUFSIZE];

    // Read numbers from file into buffer.
    std::ifstream infile;
    infile.open(filename);
    double d;
    char delim;
    while (! infile.eof())
        {
        std::string line;
        getline(infile, line);
        int temp_cols = 0;
        std::istringstream stream(line);
        while(! stream.eof()) {
            stream >> d;
            stream >> delim;
            buff[cols*rows+temp_cols++] = d;
        }

        if (temp_cols == 0)
            continue;

        if (cols == 0)
            cols = temp_cols;

        rows++;
        }

    infile.close();

    rows--;

    // Populate matrix with numbers.
    Eigen::MatrixXd result(rows,cols);
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            result(i,j) = buff[ cols*i+j ];

    return result;

}

Eigen::VectorXd removeElem(Eigen::VectorXd vec, unsigned int toRemove)
{
    unsigned int size = vec.size()-1;
    if(toRemove < size)
        vec.segment(toRemove, size - toRemove) = vec.segment(toRemove+1, size-toRemove);

    vec.conservativeResize(size);
    return vec;
}

Eigen::MatrixXd removeRow(Eigen::MatrixXd matrix, unsigned int rowToRemove)
{
    unsigned int numRows = matrix.rows()-1;
    unsigned int numCols = matrix.cols();

    if( rowToRemove < numRows )
        matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) =
        matrix.block(rowToRemove+1,0,numRows-rowToRemove,numCols);

    matrix.conservativeResize(numRows,numCols);
    return matrix;
}

Eigen::MatrixXd removeColumn(Eigen::MatrixXd matrix, unsigned int colToRemove)
{
    unsigned int numRows = matrix.rows();
    unsigned int numCols = matrix.cols()-1;

    if( colToRemove < numCols )
        matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.block(0,colToRemove+1,numRows,numCols-colToRemove);

    matrix.conservativeResize(numRows,numCols);
    return matrix;
}

Eigen::MatrixXd removeRowColumn(Eigen::MatrixXd matrix, unsigned int toRemove)
{
    matrix = removeRow(matrix, toRemove);
    matrix = removeColumn(matrix, toRemove);
    return matrix;
}

std::vector<int> findConnectedComponents(const Eigen::MatrixXd& adjacency)
{
    std::vector<bool> visited(adjacency.rows(), false);
    std::vector<int> out(adjacency.rows());
    int curr_comp = 0;
    int curr_node = 0;
    while (std::find(visited.begin(), visited.end(), false) != visited.end()) {
        _dephtFirstSearch(adjacency, curr_node, &visited, &out, curr_comp);
        curr_comp += 1;
        auto it = std::find(visited.begin(), visited.end(), false);
        if (it == visited.end())
            break;
        else
            curr_node = std::distance(visited.begin(), it);
    }
    return out;
}

void _dephtFirstSearch(const Eigen::MatrixXd &adjacency, int curr_node,
                       std::vector<bool> *visited,
                       std::vector<int> *node2comp,
                       int curr_comp)
{
    visited->at(curr_node) = true;
    node2comp->at(curr_node) = curr_comp;
    for (int k = 0; k < adjacency.cols(); k++)
    {
        if ((adjacency(curr_node, k) > 0) && !(visited->at(k)))
            _dephtFirstSearch(adjacency, k, visited, node2comp, curr_comp);
    }
    return;
}

double matrix_normal_prec_lpdf(
    Eigen::MatrixXd x, Eigen::MatrixXd m, Eigen::MatrixXd A,
    Eigen::MatrixXd B)
{
    double out = 0.5 * A.rows() * stan::math::log_determinant_spd(B);
    out += 0.5 * B.rows() * stan::math::log_determinant_spd(A);
    double exp = (B * (x - m).transpose() * A * (x - m)).trace();
    out -= 0.5 * exp;
    return out;
}
}

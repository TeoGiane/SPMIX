#include "utils.h"

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

/*Eigen::VectorXd InvAlr(Eigen::VectorXd x, bool padded_zero) {
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
}*/

std::vector<std::vector<double>> readDataFromCSV(std::string filename) {

    // Check if file exists
    if (!std::filesystem::exists(filename))
    	throw std::runtime_error("Input file does not exists.");

    // Start reading file
    std::cout << "readDataFromCSV ... ";
    std::ifstream infile(filename);
    std::map<int, std::vector<double>> out;

    int group;
    double datum;
    std::string line;
    char delim;

    // Skip header
    std::getline(infile, line);

    while (std::getline(infile, line)) {
      std::istringstream iss(line);
      if (!(iss >> group >> delim >> datum) ) { break; }
      out[group].push_back(datum);
    }

    // Maps are sorted
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

    std::cout << "done!" << std::endl;
    return data;
}

Eigen::MatrixXd readMatrixFromCSV(std::string filename) {

    // Check if file exists
    if (!std::filesystem::exists(filename))
    	throw std::runtime_error("Input file does not exists.");

    std::cout << "readMatrixFromCSV ... ";
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

    std::cout << "done!" << std::endl;
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

Rcpp::RawVector str2raw(const std::string & str){

    Rcpp::RawVector out(str.size());
    std::copy(str.begin(), str.end(), out.begin());
    return out;
}

std::string raw2str(const Rcpp::RawVector & raw_vect){

    std::string out;
    for (auto elem : raw_vect)
        out.push_back(elem);
    return out;
}

std::pair<double,double> range(const std::vector<std::vector<double>> & vectvect) {

    int N = vectvect.size();
    std::vector<double> mins(N), maxs(N);
    for (int i = 0; i < N; ++i){
        auto [minit, maxit] = std::minmax_element(vectvect[i].begin(), vectvect[i].end());
        mins[i] = *minit; maxs[i] = *maxit;
    }
    return std::make_pair(*std::min_element(mins.begin(), mins.end()),
                          *std::max_element(maxs.begin(), maxs.end()));
}

int min(const std::vector<double> & vect) {
  int min_idx = 0;
  for(int j = 1; j < vect.size(); j++) {
    if(vect[j] <= vect[min_idx]) {
      min_idx = j;
    }
  }
  return min_idx;
}

std::vector<Eigen::VectorXd> post_lpdf_from_state(const spmix::UnivariateState & current_state, const std::vector<std::vector<double>> & data) {
    // Get dimensions
    int numGroups = current_state.groupparams_size();
    int numComponents = current_state.num_components();
    // Pre-cache component parameters
    std::vector<double> means(numComponents), stdevs(numComponents);
    for (int h = 0; h < numComponents; ++h) {
        means[h] = current_state.atoms(h).mean();
        stdevs[h] = current_state.atoms(h).stdev();
    }
    // Compute posterior lpdfs for each group
    std::vector<Eigen::VectorXd> post_lpdfs(numGroups);
    for (int i = 0; i < numGroups; ++i) {
        // Get number of data points
        int numDataPoints = data[i].size();
        post_lpdfs[i] = Eigen::VectorXd(numDataPoints);        
        // Pre-compute log weights for this group
        std::vector<double> log_weights(numComponents);
        for (int h = 0; h < numComponents; ++h) {
            log_weights[h] = std::log(current_state.groupparams(i).weights(h));
        }
        // Compute posterior lpdf for each data point
        for (int j = 0; j < numDataPoints; ++j) {
            Eigen::VectorXd log_probs(numComponents);
            for (int h = 0; h < numComponents; ++h) {
                log_probs(h) = log_weights[h] + stan::math::normal_lpdf(data[i][j], means[h], stdevs[h]);
            }
            post_lpdfs[i](j) = stan::math::log_sum_exp(log_probs);
        }
    }
    return post_lpdfs;
}

std::vector<Eigen::VectorXd> pred_lpdf_from_state(const spmix::UnivariateState & current_state, const Eigen::VectorXd & grid) {
    // Get dimensions
    int numGroups = current_state.groupparams_size();
    int numComponents = current_state.num_components();
    int gridSize = grid.size();
    
    // Pre-cache component parameters
    std::vector<double> means(numComponents), stdevs(numComponents); //, log_weights_component(numComponents);
    for (int h = 0; h < numComponents; ++h) {
        means[h] = current_state.atoms(h).mean();
        stdevs[h] = current_state.atoms(h).stdev();
    }
    
    // Pre-compute component lpdf matrix: [gridSize x numComponents]
    Eigen::MatrixXd component_lpdfs(gridSize, numComponents);
    for (int j = 0; j < gridSize; ++j) {
        for (int h = 0; h < numComponents; ++h) {
            component_lpdfs(j, h) = stan::math::normal_lpdf(grid(j), means[h], stdevs[h]);
        }
    }
    
    // Compute predictive lpdfs for each group
    std::vector<Eigen::VectorXd> pred_lpdfs(numGroups);
    for (int i = 0; i < numGroups; ++i) {
        // Pre-compute log weights for this group
        std::vector<double> log_weights(numComponents);
        for (int h = 0; h < numComponents; ++h) {
            log_weights[h] = std::log(current_state.groupparams(i).weights(h));
        }
        
        // Vectorize: add log weights to each column and compute log_sum_exp per row
        Eigen::MatrixXd weighted_lpdfs = component_lpdfs;
        for (int h = 0; h < numComponents; ++h) {
            weighted_lpdfs.col(h).array() += log_weights[h];
        }
        
        pred_lpdfs[i] = Eigen::VectorXd(gridSize);
        for (int j = 0; j < gridSize; ++j) {
            pred_lpdfs[i](j) = stan::math::log_sum_exp(weighted_lpdfs.row(j));
        }
    }
    
    return pred_lpdfs;
}

// std::vector<std::vector<double>> subsample_data(const std::vector<std::vector<double>>& data, int num_samples, std::mt19937_64& rng) {

//     if (num_samples < 0) {
//         throw std::invalid_argument("Number of samples cannot be negative.");
//     }

//     std::vector<std::vector<double>> subsampled_data;
//     subsampled_data.reserve(data.size());

//     for (const auto & group_data : data) {
//         int original_size = group_data.size();
//         if (num_samples >= original_size) {
//             subsampled_data.push_back(group_data);
//         } else {
//             std::vector<double> temp_data = group_data;
//             std::shuffle(temp_data.begin(), temp_data.end(), rng);
//             subsampled_data.emplace_back(temp_data.begin(), temp_data.begin() + num_samples);
//         }
//     }

//     return subsampled_data;
// }

} // namespace utils

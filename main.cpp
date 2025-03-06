#include <iostream>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <ctime>
#include <algorithm>
#include <stdexcept>
//#include <core/lattice/lat-hal.h>
//#include <pke/openfhe.h>

#include <getopt.h>
#include <execinfo.h>
#include <signal.h>
#include <unistd.h>
#include <math/dftransform.h>
#include "PSA-cryptocontext.h"
using namespace std;
// For convenience, define a type alias for a matrix of doubles.
using Mat = std::vector<std::vector<double>>;

// -------------------------
// Logistic Regression Helpers
// -------------------------

// Sigmoid activation function
double sigmoid(double z) {
    return 1.0 / (1.0 + std::exp(-z));
}

// Compute logistic (cross-entropy) loss for a dataset.
double computeLoss(const Mat &X, const std::vector<double> &y, const std::vector<double> &weights) {
    size_t m = X.size();
    double loss = 0.0;
    const double eps = 1e-15;
    for (size_t i = 0; i < m; i++) {
        double z = 0.0;
        for (size_t j = 0; j < weights.size(); j++) {
            z += X[i][j] * weights[j];
        }
        double pred = sigmoid(z);
        loss += - (y[i] * std::log(pred + eps) + (1 - y[i]) * std::log(1 - pred + eps));
    }
    return loss / m;
}

// Train logistic regression using gradient descent with validation-based checkpointing.
// Returns the best weights (lowest validation loss).
std::vector<double> trainLogisticRegression(const Mat &X_train, const std::vector<double> &y_train,
                                              const Mat &X_val, const std::vector<double> &y_val,
                                              double learningRate = 0.1, int epochs = 1000) {
    size_t m = X_train.size();
    if (m == 0) return {};
    size_t n = X_train[0].size();

    // Initialize weights to small random values.
    std::vector<double> weights(n, 0.0);
    std::srand(static_cast<unsigned>(std::time(0)));
    for (size_t j = 0; j < n; j++) {
        weights[j] = ((double) std::rand() / RAND_MAX - 0.5) * 0.1;
    }
    
    std::vector<double> bestWeights = weights;
    double bestValLoss = std::numeric_limits<double>::max();

    // Gradient descent training loop.
    for (int epoch = 0; epoch < epochs; epoch++) {
        std::vector<double> gradients(n, 0.0);
        double trainLoss = 0.0;
        const double eps = 1e-15;
        // Loop over each training sample.
        for (size_t i = 0; i < m; i++) {
            double z = 0.0;
            for (size_t j = 0; j < n; j++) {
                z += X_train[i][j] * weights[j];
            }
            double pred = sigmoid(z);
            double error = pred - y_train[i];
            for (size_t j = 0; j < n; j++) {
                gradients[j] += error * X_train[i][j];
            }
            trainLoss += - (y_train[i] * std::log(pred + eps) + (1 - y_train[i]) * std::log(1 - pred + eps));
        }
        trainLoss /= m;

        // Update weights.
        for (size_t j = 0; j < n; j++) {
            weights[j] -= learningRate * gradients[j] / m;
        }

        // Compute validation loss.
        double valLoss = computeLoss(X_val, y_val, weights);

        // Save best weights if validation loss improves.
        if (valLoss < bestValLoss) {
            bestValLoss = valLoss;
            bestWeights = weights;
        }

        if (epoch % 100 == 0) {
            std::cout << "Epoch " << epoch 
                      << ", Training Loss: " << trainLoss 
                      << ", Validation Loss: " << valLoss << std::endl;
        }
    }
    return bestWeights;
}

// Compute accuracy on a given dataset.
double computeAccuracy(const Mat &X, const std::vector<double> &y, const std::vector<double> &weights) {
    size_t correct = 0;
    size_t m = X.size();
    size_t n = weights.size();
    
    for (size_t i = 0; i < m; i++) {
        double z = 0.0;
        for (size_t j = 0; j < n; j++) {
            z += X[i][j] * weights[j];
        }
        int predictedLabel = sigmoid(z) >= 0.5 ? 1 : 0;
        if (predictedLabel == static_cast<int>(y[i])) {
            correct++;
        }
    }
    return static_cast<double>(correct) / m;
}

// Function to shuffle data (features and labels) in unison.
void shuffleData(Mat &X, std::vector<double> &y) {
    size_t m = X.size();
    std::vector<size_t> indices(m);
    for (size_t i = 0; i < m; i++) {
        indices[i] = i;
    }
    std::random_shuffle(indices.begin(), indices.end());
    Mat X_shuffled;
    std::vector<double> y_shuffled;
    for (size_t i = 0; i < m; i++) {
        X_shuffled.push_back(X[indices[i]]);
        y_shuffled.push_back(y[indices[i]]);
    }
    X = X_shuffled;
    y = y_shuffled;
}

// Compute the number of unique monomials for N variables:
// num_monomial = (N+2 choose 3) = (N+2)*(N+1)*N / 6
int computeNumMonomials(int N) {
    return (N + 2) * (N + 1) * N / 6;
}

// Constructs the 3-row x-dependent matrix.
// Each row has length = prefix_size + num_monomial, where prefix_size = N + 1.
vector<vector<double>> constructXDependentMatrix(const vector<double>& x) {
    int N = x.size();
    int prefix_size = N + 1;
    int num_monomial = computeNumMonomials(N);
    int total_size = prefix_size + num_monomial;

    // Build monomial vectors (v1, v2, v3) from all unique triples (i,j,k) with 0 <= i <= j <= k < N.
    vector<double> v1, v2, v3;
    for (int i = 0; i < N; i++) {
        for (int j = i; j < N; j++) {
            for (int k = j; k < N; k++) {
                v1.push_back(x[i]);  // First factor.
                v2.push_back(x[j]);  // Second factor.
                v3.push_back(x[k]);  // Third factor.
            }
        }
    }

    // Create a 3 x total_size matrix.
    vector<vector<double>> xMatrix(3, vector<double>(total_size, 0.0));

    // Row 0: first prefix_size entries are ones, then v1.
    for (int i = 0; i < prefix_size; i++) {
        xMatrix[0][i] = 1.0;
    }
    for (int i = 0; i < num_monomial; i++) {
        xMatrix[0][prefix_size + i] = v1[i];
    }

    // Row 1: first prefix_size entries are ones, then v2.
    for (int i = 0; i < prefix_size; i++) {
        xMatrix[1][i] = 1.0;
    }
    for (int i = 0; i < num_monomial; i++) {
        xMatrix[1][prefix_size + i] = v2[i];
    }

    // Row 2: first entry is 1, next N entries are the original x's, then v3.
    xMatrix[2][0] = 1.0;
    for (int i = 0; i < N; i++) {
        xMatrix[2][1 + i] = x[i];
    }
    for (int i = 0; i < num_monomial; i++) {
        xMatrix[2][prefix_size + i] = v3[i];
    }

    return xMatrix;
}

// Constructs the coefficient vector using the weights w.
// The length is prefix_size + num_monomial.
vector<double> constructCoefficientVector(const vector<double>& w) {
    int N = w.size();
    int prefix_size = N + 1;
    int num_monomial = computeNumMonomials(N);
    int total_size = prefix_size + num_monomial;

    vector<double> coeff(total_size, 0.0);

    // First prefix_size entries: first entry 0.5, then w[i]/4.
    coeff[0] = 0.5;
    for (int i = 0; i < N; i++) {
        coeff[1 + i] = w[i] / 4.0;
    }

    // For the monomial coefficients:
    int index = prefix_size;
    for (int i = 0; i < N; i++) {
        for (int j = i; j < N; j++) {
            for (int k = j; k < N; k++) {
                int m = 0; // Multiplicity factor.
                if (i == j && j == k) {
                    m = 1;
                } else if (i == j || j == k || i == k) {
                    m = 3;
                } else {
                    m = 6;
                }
                coeff[index++] = m * (w[i] * w[j] * w[k]);
            }
        }
    }

    return coeff;
}

void handler(int sig) {
    void *array[10];
    size_t size;

    // get void*'s for all entries on the stack
    size = backtrace(array, 10);

    // print out all the frames to stderr
    fprintf(stderr, "Error: signal %d:\n", sig);
    backtrace_symbols_fd(array, size, STDERR_FILENO);
    exit(1);
}

double slap(int argc, char **argv, Mat& inputmatrix) {
    // Set up signal handling.
    signal(SIGSEGV, handler);
    std::cout << "Hello, World!" << std::endl;
    
    // Default parameters.
    unsigned int plain_bits = 40; // log t
    // Set number of users from the inputmatrix: all rows except the last one.
    unsigned int num_users = inputmatrix.size() - 1;
    unsigned int iters = 1;
    unsigned int k_prime = 1;
    Scheme scheme1 = NS;
    unsigned int N = 1;
    
    // Process command-line arguments.
    int c;
    while ((c = getopt(argc, argv, "t:n:i:k:N:")) != -1) {
        switch(c) {
            case 't':
                plain_bits = atoi(optarg);
                break;
            case 'n':
                // Ignore user-specified n; we use inputmatrix.size()-1.
                break;
            case 'i':
                iters = atoi(optarg);
                break;
            case 'k':
                k_prime = atoi(optarg);
                break;
            case 'N':
                N = atoi(optarg);
                break;
            default:
                std::cout << "Invalid argument: " << c;
                if (optarg != nullptr) {
                    std::cout << ' ' << optarg;
                }
                std::cout << std::endl;
                return 1;
        }
    }
    
    if (!plain_bits) {
        throw std::runtime_error("Must have nonempty plaintext space");
    }
    if (num_users == 0) {
        throw std::runtime_error("Must have at least some users (x-dependent rows)");
    }
    if (!iters) {
        throw std::runtime_error("Must have at least some iterations");
    }
    
    unsigned int MAX_CTEXTS_DEFAULT = 20;
    PSACryptocontext pp(plain_bits, num_users, iters, scheme1);
    
    std::vector<double> poly_noise_times;
    std::vector<double> poly_enc_times;
    
    // Set up the polynomial encryption environment.
    pp.PolynomialEnvSetup(poly_noise_times, poly_enc_times);
    
    // Create an exponent vector of ones with length equal to the number of columns.
    std::vector<double> expvec(inputmatrix[0].size(), 1);
    
    // Encrypt each x-dependent row (all rows except the last one).
    for (int i = 0; i < num_users; i++) {
        pp.PolynomialEncryption(inputmatrix[i], expvec, i, poly_noise_times, poly_enc_times);
    }
    
    std::vector<double> decrypt_times;
    
    // Use the last row of inputmatrix as the coefficient (constants) vector.
    std::vector<double> constants = inputmatrix[inputmatrix.size() - 1];
    
    // Decrypt using the coefficient vector.
    std::vector<double> outputvec = pp.PolynomialDecryption(constants, 1, decrypt_times);
    
    std::cout << "Final output: " << outputvec << std::endl;
    
    // Sum the output vector and return the result.
    double sum = std::accumulate(outputvec.begin(), outputvec.end(), 0.0);
    return sum;
}

double sigmoidapprox(int argc, char **argv, std::vector<double>& x, std::vector<double>& w) {
    Mat xMatrix = constructXDependentMatrix(x);
    std::vector<double> coeff = constructCoefficientVector(w);
    vector<vector<double>> inputMatrix = xMatrix;
    inputMatrix.push_back(coeff);
    return slap(argc, argv, inputMatrix);
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <data_file.csv>" << std::endl;
        return 1;
    }

    std::string dataFile = argv[1];

    Mat featuresData;
    std::vector<double> labelsData;
    std::vector<std::string> featureNames;
    std::string labelName;

    std::ifstream file(dataFile);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << dataFile << std::endl;
        return 1;
    }
    std::cout << "Opened file: " << dataFile << std::endl;

    // Read header line.
    std::string headerLine;
    if (!std::getline(file, headerLine)) {
        std::cerr << "Error: File is empty." << std::endl;
        return 1;
    }
    std::istringstream headerStream(headerLine);
    std::vector<std::string> headers;
    std::string colName;
    while (std::getline(headerStream, colName, ',')) {
        headers.push_back(colName);
    }
    if (headers.size() < 2) {
        std::cerr << "Error: File must contain at least one feature and one label column." << std::endl;
        return 1;
    }
    labelName = headers.back();
    headers.pop_back();
    featureNames = headers;
    size_t numFeatures = featureNames.size();
    std::cout << "Found " << numFeatures << " feature columns." << std::endl;

    // Read data rows.
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::istringstream ss(line);
        std::vector<double> row;
        std::string value;
        while (std::getline(ss, value, ',')) {
            try {
                row.push_back(std::stod(value));
            } catch (...) {
                row.clear();
                break;
            }
        }
        if (row.size() != numFeatures + 1) {
            //std::cerr << "Warning: Skipping inconsistent row." << std::endl;
            continue;
        }
        std::vector<double> features(row.begin(), row.end() - 1);
        double label = row.back();
        featuresData.push_back(features);
        labelsData.push_back(label);
    }
    file.close();
    if (featuresData.empty()) {
        std::cerr << "Error: No data loaded." << std::endl;
        return 1;
    }
    std::cout << "Loaded " << featuresData.size() << " samples." << std::endl;

    // Shuffle and split data: 60% train, 20% validation, 20% test.
    shuffleData(featuresData, labelsData);
    size_t total = featuresData.size();
    size_t trainSize = total * 0.6;
    size_t valSize = total * 0.2;
    size_t testSize = total - trainSize - valSize;
    Mat X_train(featuresData.begin(), featuresData.begin() + trainSize);
    std::vector<double> y_train(labelsData.begin(), labelsData.begin() + trainSize);
    Mat X_val(featuresData.begin() + trainSize, featuresData.begin() + trainSize + valSize);
    std::vector<double> y_val(labelsData.begin() + trainSize, labelsData.begin() + trainSize + valSize);
    Mat X_test(featuresData.begin() + trainSize + valSize, featuresData.end());
    std::vector<double> y_test(labelsData.begin() + trainSize + valSize, labelsData.end());

    std::cout << "Training samples: " << X_train.size() 
              << ", Validation samples: " << X_val.size() 
              << ", Test samples: " << X_test.size() << std::endl;

    // 1. Train logistic regression to get best weights.
    double learningRate = 0.1;
    int epochs = 1000;
    std::vector<double> bestWeights = trainLogisticRegression(X_train, y_train, X_val, y_val, learningRate, epochs);

    // 2. Compute testing accuracy using the standard logistic regression model.
    double standardAccuracy = computeAccuracy(X_test, y_test, bestWeights);
    std::cout << "\nStandard Logistic Regression Test Accuracy: " 
              << standardAccuracy * 100.0 << "%" << std::endl;

    // 3. For each test sample, use sigmoidapprox() to get an approximate prediction.
    // Here we assume that sigmoidapprox() returns an approximation of sigmoid(w'x).
    // We threshold the approximation at 0.5 to assign a label.
    int correctApprox = 0;
    for (size_t i = 0; i < X_test.size(); i++) {
        // For each test sample, use its feature vector.
        std::vector<double> testSample = X_test[i];
        // Here, we pass dummy command-line parameters.
        int dummy_argc = 1;
        char* dummy_argv[] = { (char*)"sigmoidapprox", nullptr };
        double approxOutput = sigmoidapprox(dummy_argc, dummy_argv, testSample, bestWeights);
        int predictedLabel = (approxOutput >= 0.5) ? 1 : 0;
        if (predictedLabel == static_cast<int>(y_test[i])) {
            correctApprox++;
        }
    }
    double approxAccuracy = static_cast<double>(correctApprox) / X_test.size();
    std::cout << "Sigmoid Approximation Test Accuracy (via slap): " 
              << approxAccuracy * 100.0 << "%" << std::endl;

    // 4. Compare accuracies.
    std::cout << "\nComparison:" << std::endl;
    std::cout << "Standard Logistic Regression Accuracy: " << standardAccuracy * 100.0 << "%" << std::endl;
    std::cout << "Sigmoid Approximation Accuracy: " << approxAccuracy * 100.0 << "%" << std::endl;

    return 0;
}
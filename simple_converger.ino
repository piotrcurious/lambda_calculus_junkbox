#include <vector>
#include <algorithm>
#include <numeric>

// Function to calculate the polynomial coefficients
std::vector<float> polynomialFit(const std::vector<float>& x, const std::vector<float>& y, int order) {
    size_t n = x.size();
    int parameters = order + 1;
    std::vector<std::vector<float>> B(n, std::vector<float>(parameters));
    std::vector<float> S(parameters * 2 - 1);
    std::vector<float> beta(parameters);

    // Fill matrix B and vector S
    for (size_t i = 0; i < n; i++) {
        for (int j = 0; j < parameters; j++) {
            B[i][j] = pow(x[i], j);
        }
    }

    for (int i = 0; i < parameters * 2 - 1; i++) {
        for (size_t j = 0; j < n; j++) {
            S[i] += pow(x[j], i);
        }
    }

    // Solve the normal equations using lambda calculus
    auto lambda = & {
        return std::inner_product(B.begin(), B.end(), B.begin(), 0.0f,
                                  std::plus<float>(),
                                  i, j {
                                      return lhs[i] * rhs[j];
                                  });
    };

    std::vector<std::vector<float>> A(parameters, std::vector<float>(parameters));
    for (int i = 0; i < parameters; i++) {
        for (int j = 0; j < parameters; j++) {
            A[i][j] = lambda(i, j);
        }
    }

    // Calculate vector beta (coefficients)
    for (int i = 0; i < parameters; i++) {
        for (size_t j = 0; j < n; j++) {
            beta[i] += B[j][i] * y[j];
        }
    }

    // Solve the linear system A * coefficients = beta
    std::vector<float> coefficients = solveLinearSystem(A, beta);
    return coefficients;
}


/*
// Function to solve the linear system Ax=b using Gaussian elimination
std::vector<float> solveLinearSystem(std::vector<std::vector<float>>& A, std::vector<float>& b) {
    size_t n = A.size();
    for (size_t i = 0; i < n; i++) {
        // Search for maximum in this column
        float maxEl = abs(A[i][i]);
        size_t maxRow = i;
        for (size_t k = i + 1; k < n; k++) {
            if (abs(A[k][i]) > maxEl) {
                maxEl = abs(A[k][i]);
                maxRow = k;
            }
        }

        // Swap maximum row with current row (column by column)
        for (size_t k = i; k < n; k++) {
            std::swap(A[maxRow][k], A[i][k]);
        }
        std::swap(b[maxRow], b[i]);

        // Make all rows below this one 0 in current column
        for (size_t k = i + 1; k < n; k++) {
            float c = -A[k][i] / A[i][i];
            for (size_t j = i; j < n; j++) {
                if (i == j) {
                    A[k][j] = 0;
                } else {
                    A[k][j] += c * A[i][j];
                }
            }
            b[k] += c * b[i];
        }
    }

    // Solve equation Ax=b for an upper triangular matrix A
    std::vector<float> x(n);
    for (int i = n - 1; i >= 0; i--) {
        x[i] = b[i] / A[i][i];
        for (int k = i - 1; k >= 0; k--) {
            b[k] -= A[k][i] * x[i];
        }
    }
    return x;
}
*/
#include <vector>

// Function to perform LU Factorization on matrix A
void LUFactorization(std::vector<std::vector<float>>& A, std::vector<std::vector<float>>& L, std::vector<std::vector<float>>& U) {
    size_t n = A.size();
    L = std::vector<std::vector<float>>(n, std::vector<float>(n, 0));
    U = std::vector<std::vector<float>>(n, std::vector<float>(n, 0));

    for (size_t i = 0; i < n; i++) {
        // Upper Triangular
        for (size_t k = i; k < n; k++) {
            float sum = 0;
            for (size_t j = 0; j < i; j++)
                sum += (L[i][j] * U[j][k]);
            U[i][k] = A[i][k] - sum;
        }

        // Lower Triangular
        for (size_t k = i; k < n; k++) {
            if (i == k)
                L[i][i] = 1; // Diagonal as 1
            else {
                float sum = 0;
                for (size_t j = 0; j < i; j++)
                    sum += (L[k][j] * U[j][i]);
                L[k][i] = (A[k][i] - sum) / U[i][i];
            }
        }
    }
}

// Function to solve the linear system Ax=b using LU Factorization
std::vector<float> solveLinearSystem(std::vector<std::vector<float>>& A, std::vector<float>& b) {
    size_t n = A.size();
    std::vector<std::vector<float>> L, U;

    // Perform LU Factorization
    LUFactorization(A, L, U);

    // Solve Ly=b for y using forward substitution
    std::vector<float> y(n, 0);
    for (size_t i = 0; i < n; i++) {
        y[i] = b[i];
        for (size_t j = 0; j < i; j++) {
            y[i] -= L[i][j] * y[j];
        }
        y[i] /= L[i][i];
    }

    // Solve Ux=y for x using back substitution
    std::vector<float> x(n, 0);
    for (int i = n - 1; i >= 0; i--) {
        x[i] = y[i];
        for (size_t j = i + 1; j < n; j++) {
            x[i] -= U[i][j] * x[j];
        }
        x[i] /= U[i][i];
    }

    return x;
}



// Example usage
int main() {
    // Assuming we have a buffer of data points (x, y)
    std::vector<float> x; // x-values
    std::vector<float> y; // corresponding y-values

    // Fill x and y with your data
    // ...

    // Fit the data with a polynomial of order 3
    std::vector<float> coefficients = polynomialFit(x, y, 3);

    // Output the coefficients
    for (float coeff : coefficients) {
        Serial.println(coeff);
    }

    return 0;
}

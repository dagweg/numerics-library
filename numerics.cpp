#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

// Matrix Processor class
template <class T>
class MP
{
public:
    // finds minor
    // input ci - the row to cross out
    // input cj - the column to cross out
    static vector<vector<T>> minor(vector<vector<T>> matrix, int ci, int cj)
    {
        // invalid case return the original matrix
        if (ci > matrix.size() || cj > matrix.size())
        {
            return matrix;
        }

        vector<vector<T>> newMatrix;

        for (int i = 0; i < matrix.size(); i++)
        {
            vector<T> row;

            // skip this row if its to be crossed out
            if (i == ci)
                continue;

            for (int j = 0; j < matrix.size(); j++)
            {
                // skip this column if its to be crossed out
                if (j == cj)
                    continue;

                // the rest add to the innerMatrix
                row.push_back(matrix[i][j]);
            }
            newMatrix.push_back(row);
        }

        return newMatrix;
    };

    // determinant for a square matrix of order N
    static T detNxN(vector<vector<T>> matrix)
    {

        T determinant = 0;

        switch (matrix.size())
        {
        case 0:
            return 0;
        case 1:
            return matrix[0][0];
        case 2:
            return det2x2(matrix);
        case 3:
            return det3x3(matrix);
        default:
            return _detNxN(matrix);
            break;
        }

        return determinant;
    };

    // determinant for a square matrix of order 3
    static T det3x3(vector<vector<T>> matrix)
    {
        T determinant = 0;

        for (int i = 0; i < matrix.size(); i++)
        {
            determinant += matrix[0][i] * pow(-1, i) * det2x2(minor(matrix, 0, i));
        }

        return determinant;
    }

    // determinant for a square matrix of order 2
    static T det2x2(vector<vector<T>> matrix)
    {
        return matrix[0][0] * matrix[1][1] - matrix[1][0] * matrix[0][1];
    }

    // finds cofactor of a matrix
    static vector<vector<T>> cofactor(vector<vector<T>> matrix)
    {
        vector<vector<T>> cofactorMatrix(vector<vector<T>>(matrix.size(), vector<T>(matrix.size(), 0.0)));

        for (int i = 0; i < matrix.size(); i++)
        {
            for (int j = 0; j < matrix.size(); j++)
            {
                cofactorMatrix[i][j] = pow(-1, i + j) * detNxN(minor(matrix, i, j));
            }
        }
        return cofactorMatrix;
    };

    // finds transpose of a matrix
    static vector<vector<T>> transpose(vector<vector<T>> matrix)
    {
        vector<vector<T>> transposedMatrix = matrix;

        for (int i = 0; i < matrix.size(); i++)
        {
            for (int j = i + 1; j < matrix.size(); j++)
            {
                swap(transposedMatrix[i][j], transposedMatrix[j][i]);
            }
        }

        return transposedMatrix;
    };

    // scales each entry of the matrix with some scalar
    static vector<vector<T>> scaleMatrix(vector<vector<T>> matrix, T scalar)
    {
        vector<vector<T>> scaledMatrix = matrix;
        for (int i = 0; i < matrix.size(); i++)
        {
            for (int j = 0; j < matrix.size(); j++)
            {
                scaledMatrix[i][j] *= scalar;
            }
        }
        return scaledMatrix;
    }

    // calcuates inverse
    static vector<vector<T>> inverse(vector<vector<T>> matrix)
    {
        vector<vector<T>> adjugate = transpose(cofactor(matrix));
        return scaleMatrix(adjugate, 1 / detNxN(matrix));
    };

    // Ax = B
    // A - coeff_matrix
    // Xo - init_guess for x
    // B - rhs_vector
    // error_tolerance - epsilon
    // returns the approximated x-values (x0,x1,x2,....xn)
    static vector<T> gaussJacobiSolver(vector<vector<T>> coeff_matrix, vector<T> init_guess, vector<T> rhs_vector, double epsilon)
    {
        // STEP 1: reorder diagonals such that they aren't zero
        for (int i = 0; i < coeff_matrix.size(); i++)
        {
            if (coeff_matrix[i][i] != 0)
                continue;

            for (int j = i + 1; j < coeff_matrix.size(); j++)
            {
                if (coeff_matrix[j][j] == 0)
                    continue;

                // we found the one, so swap the rows
                vector<T> temp_coeff = coeff_matrix[i];
                coeff_matrix[i] = coeff_matrix[j];
                coeff_matrix[j] = temp_coeff;

                // Swap also the corresponding rhs_vector i.e B
                T temp_rhs = rhs_vector[i];
                rhs_vector[i] = rhs_vector[j];
                rhs_vector[j] = temp_rhs;
            }
        }

        vector<T> ret; // this is what we return

        // STEP 2: approximate till tolerance is reached
        double cur_tolerance = INFINITY;

        while (cur_tolerance > epsilon)
        {
            // initialize an empty array with length of the guess
            vector<T> next_guess(init_guess.size(), 0);

            for (int i = 0; i < next_guess.size(); i++)
            {

                double summation = 0;
                for (int j = 0; j < next_guess.size(); j++)
                {
                    // skip where i equals j
                    if (i == j)
                        continue;

                    summation += coeff_matrix[i][j] * init_guess[j];
                }

                // using gauss-jacobi formula for next-guess
                next_guess[i] = 1 / coeff_matrix[i][i] * (rhs_vector[i] - summation);
            }

            vector<T> difference;
            for (int i = 0; i < init_guess.size(); i++)
            {
                difference.emplace_back(next_guess[i] - init_guess[i]);
            }

            cur_tolerance = euclideanNorm(difference); // to compare to the epsilon value (error tolerance)
            ret = init_guess = next_guess;
        }

        return ret;
    }

    static double euclideanNorm(const vector<double> &vectr)
    {
        double sum = 0.0;
        for (int i = 0; i < vectr.size(); i++)
        {
            sum += pow(vectr[i], 2);
        }
        return sqrt(sum);
    }

private:
    // recursively calculates the determinant of any sized NxN matrix
    static T _detNxN(vector<vector<T>> matrix)
    {
        if (matrix.size() == 3)
            return det3x3(matrix);

        T determinant = 0;
        for (int j = 0; j < matrix.size(); j++)
        {
            determinant += matrix[0][j] * pow(-1, j) * _detNxN(minor(matrix, 0, j));
        }
        return determinant;
    }
};

template <typename T>
void printMatrix(const vector<vector<T>> &matrix)
{
    for (const vector<T> &v : matrix)
    {
        for (const T &t : v)
        {
            cout << t << "  ";
        }
        cout << endl;
    }
}

template <typename T>
void printVector(const vector<T> &vectr)
{
    for (const T &t : vectr)
    {
        cout << t << "\n";
    }
}

int main(void)
{
    MP<double> m;

    // Matrix Invertion example
    vector<vector<double>> testMatrix = {{1, 2, 3}, {3, 2, 1}, {2, 1, 3}};

    vector<vector<double>> invertedMatrix = m.inverse(testMatrix);

    cout << "Inverse of the matrix:\n";
    printMatrix<double>(testMatrix);

    cout << "\nis:\n";
    printMatrix<double>(invertedMatrix);

    // Gauss Jacobi SLE approximation example
    vector<vector<double>> coeff_matrix = {
        {10, 3, 1},
        {2, -10, 3},
        {1, 3, 10},
    };

    vector<double> B = {14, -5, 14};
    vector<double> xo = {0, 0, 0}; // initial guess
    double epsilon = 1e-2;
    vector<double> soln = m.gaussJacobiSolver(coeff_matrix, xo, B, epsilon);

    cout << "\nThe gauss jacobi approximation for the system \n";
    printMatrix<double>(coeff_matrix);
    cout << "\nWith an initial guess of:\n";
    printVector<double>(xo);
    cout << "\nand error tolerance: " << epsilon << endl;
    cout << "\nis:\n";
    printVector<double>(soln);

    return 0;
}
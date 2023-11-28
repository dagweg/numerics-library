#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

// Matrix Processor class
template <typename T>
class MP
{
private:
    static T _determinant;

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

int main(void)
{
    MP<double> m;

    // example
    vector<vector<double>> testMatrix = {{1, 2, 3}, {3, 2, 1}, {2, 1, 3}};

    vector<vector<double>> invertedMatrix = m.inverse(testMatrix);

    for (vector<double> &v : invertedMatrix)
    {
        for (double &d : v)
            cout << d << " ";
        cout << endl;
    }

    return 0;
}
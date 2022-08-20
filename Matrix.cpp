#include "Matrix.hpp"
#include <fstream>
#include <stdexcept>
#include <ctime>
#include <cassert>
#include <cmath>
#include <memory>

Matrix::Matrix(int m, int n, bool random) {
    matrix = create_matrix(m, n);
    for (int i = 0; i < _row; ++i) {
        for (int j = 0; j < _col; ++j) {
            matrix[i * _col + j] = 0.;
        }
    }
    if (random) {
        std::srand(std::time(nullptr));
        for (int i = 0; i < _row; ++i) {
            for (int j = 0; j < _col; ++j) {
                matrix[i * _col + j] = std::rand();
            }
        }
    }
}

Matrix::~Matrix() {
    delete[] matrix;
    delete[] _U_matrix;
    delete[] _L_matrix;
//    delete[] inversematrix;
}


Matrix::Matrix(std::string file_name) {
    std::ifstream in;
    int m;
    int n;
    std::string line;
    in.open(file_name);
    if (!in.is_open()) {
        throw std::runtime_error("file_error");
    }
    in >> m;
    in >> n;
    matrix = create_matrix(m, n);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            if (in.eof()) {
                throw std::runtime_error("wrong input matrix");
            }
            in >> matrix[i * n + j];
        }
    }
}

Matrix::Matrix(int m, int n, double *array) {
    matrix = create_matrix(m, n);
    for (int i = 0; i < _row; ++i) {
        for (int j = 0; j < _col; ++j) {
            matrix[i * _col + j] = array[i * _row + j];
        }
    }

}

Matrix::Matrix(Matrix &tmp) {
    matrix = create_matrix(tmp._row, tmp._col);
    for (int i = 0; i < _row; ++i) {
        for (int j = 0; j < _col; ++j) {
            matrix[i * _col + j] = tmp.matrix[i * _row + j];
        }
    }
}

double *Matrix::create_matrix(int m, int n) {
    _row = m;
    _col = n;
    return new double[m * n];
}

void Matrix::print_matrix() {
//    std::cout << std::scientific;
    for (int i = 0; i < _row; ++i) {
        for (int j = 0; j < _col; ++j) {
            std::cout << matrix[i * _col + j] << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}

void Matrix::matrix_to_file(std::string filename) {
    std::ofstream fout;
    if (!filename.empty()) {
        fout.open(filename + ".txt");

    } else {
        fout.open("matrix_output.txt");
    }
    for (int i = 0; i < _row; ++i) {
        for (int j = 0; j < _col; ++j) {
            fout << matrix[i * _col + j] << " ";
        }
        fout << "\n";
    }
    fout << "\n";
}

Matrix Matrix::operator*(const Matrix &right) {
    if (_col != right._row) {
//        std::cout << _col << " " <<  right._row << std::endl;
        throw std::runtime_error("wrong matrix for multiplication");
    }
    Matrix result(_row, right._col);
    for (int i = 0; i < result._row; ++i) {
        for (int j = 0; j < result._col; ++j) {
            for (int k = 0; k < _col; ++k) {
                double num = matrix[i * _col + k] * right.matrix[k * right._col + j];
                if (fabs(num) < 1e-13) {
                    num = 0.;
                }
                result.matrix[i * result._col + j] += num;
            }
        }
    }
    return result;
}

void Matrix::print_array(double *array, int row, int col) const {
    std::cout << std::scientific;
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            std::cout << array[i * col + j] << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}


void Matrix::LU_decomposition() {
    if (_row != _col) {
        throw std::runtime_error("matrix is not n*n");
    }
    if (_L_matrix == nullptr) {
        _L_matrix = create_matrix(_row, _row);
        _U_matrix = create_matrix(_row, _row);
    }

    for (int i = 0; i < _row; ++i) {
        for (int j = 0; j < _col; ++j) {
            _U_matrix[i * _col + j] = matrix[i * _col + j];
        }
    }

    for (int i = 0; i < _row; ++i) {
        for (int j = 0; j < _col; ++j) {
            if (j < i) {
                _L_matrix[j * _col + i] = 0;
            } else {
                _L_matrix[j * _col + i] = _U_matrix[j * _col + i] / _U_matrix[i * _col + i];
            }
        }
    }

    for (int k = 1; k < _col; ++k) {
        for (int i = k - 1; i < _col; ++i) {
            for (int j = i; j < _col; ++j) {
                _L_matrix[j * _col + i] = _U_matrix[j * _col + i] / _U_matrix[i * _col + i];
            }
        }
        for (int i = k; i < _col; ++i) {
            for (int j = k - 1; j < _col; ++j) {
                _U_matrix[i * _col + j] =
                        _U_matrix[i * _col + j] - _L_matrix[i * _col + (k - 1)] * _U_matrix[(k - 1) * _col + j];
            }
        }
    }


}

void Matrix::print_LU_decomposition() {
    LU_decomposition();
    std::cout << std::scientific;
    for (int i = 0; i < _row; ++i) {
        for (int j = 0; j < _col; ++j) {
            std::cout << _L_matrix[i * _col + j] << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";

    for (int i = 0; i < _row; ++i) {
        for (int j = 0; j < _col; ++j) {
            std::cout << _U_matrix[i * _col + j] << " ";
        }
        std::cout << "\n";
    }
}

bool Matrix::operator==(const Matrix &right) {
    if (_row != right._row || _col != right._col) {
        return false;
    }
    for (int i = 0; i < _row; ++i) {
        for (int j = 0; j < _col; ++j) {
            if (matrix[i * _col + j] - right.matrix[i * _col + j] > 1e-13)
                return false;
        }
    }
    return true;
}

Matrix &Matrix::operator=(const Matrix &right) {
    _row = right._row;
    _col = right._col;
    delete[] matrix;
    matrix = create_matrix(_row, _col);
    for (int i = 0; i < _row; ++i) {
        for (int j = 0; j < _col; ++j) {
            matrix[i * _col + j] = right.matrix[i * _col + j];
        }
    }
    return *this;
}

int Matrix::get_row_size() {
    return _row;
}

int Matrix::get_col_size() {
    return _col;
}

double *Matrix::get_L() {
    return _L_matrix;
}

double *Matrix::get_U() {
    return _U_matrix;
}

void Matrix::Gauss(Matrix &matrix_A, Matrix &matrix_b) {
    Matrix matrix3;
    if (matrix_A._L_matrix == nullptr) {
        matrix_A.LU_decomposition();
    }

    Matrix matrix_L(matrix_A.get_row_size(), matrix_A.get_col_size(), matrix_A.get_L());
    Matrix matrix_U(matrix_A.get_row_size(), matrix_A.get_col_size(), matrix_A.get_U());

    assert(matrix_A == matrix_L * matrix_U);

    Matrix matrix_y(_row, 1);

    for (int i = 0; i < matrix_L._col; ++i) {
        double sum = 0;
        for (int j = 0; j < i; ++j) {
            sum += matrix_L.matrix[i * matrix_L._col + j] * matrix_y.matrix[j];
        }
        matrix_y.matrix[i] = matrix_b.matrix[i] - sum;

    }
    assert(matrix_L * matrix_y == matrix_b);

    for (int i = matrix_U._col - 1; i >= 0; --i) {
        double sum = 0;
        for (int j = i + 1; j < matrix_U._col; ++j) {
            sum += matrix_U.matrix[i * matrix_U._col + j] * matrix[j];
        }
        matrix[i] = (matrix_y.matrix[i] - sum) / matrix_U.matrix[i * matrix_U._col + i];
    }
    assert(matrix_A * (*this) == matrix_b);
}

Matrix Matrix::transpose() const {
    Matrix A_T(_row, _col);
    for (int i = 0; i < _row; ++i) {
        for (int j = 0; j < _col; ++j) {
            A_T.matrix[i * _col + j] = matrix[j * _col + i];
        }
    }
    return A_T;
}

double Matrix::max_eigenvalue() {
    int k = 0;

    Matrix x_0(_col, 1);
    Matrix x_1(_col, 1);

    for (int i = 0; i < x_0._row; ++i) {
        x_0.matrix[i] = 1;
    }
    double sum1 = 0, sum2 = 0;
    for (int i = 0; i < x_0._row; ++i) {
        sum1 += x_0.matrix[i] * x_0.matrix[i];
    }
    sum1 = sqrt(sum1);
    for (int i = 0; i < x_1._row; ++i) {
        x_0.matrix[i] = x_0.matrix[i] / sum1;
    }
    x_1 = *this * x_0;
    while (true) {
        for (int i = 0; i < x_0._row; ++i) {
            sum2 += x_1.matrix[i] * x_1.matrix[i];
        }
        sum2 = sqrt(sum2);
        for (int i = 0; i < x_0._row; ++i) {
            x_1.matrix[i] = x_1.matrix[i] / sum2;
            x_0.matrix[i] = x_1.matrix[i];
        }
        if (fabs(sum2 - sum1) / fabs(sum2) < 1e-7) {
            break;
        }
        x_1 = *this * x_0;
        sum1 = sum2;
        k++;
    }
    return sum2;
}

double Matrix::norm() const {
    return sqrt(fabs((transpose() * *this).max_eigenvalue()));
}

Matrix Matrix::inverse() {
    Matrix inversematrix(_row, _col);
    if (_row != _col) {
        return inversematrix;
    }
    if (fabs(determinant()) < 1e-7) {
        return inversematrix;
    }
    Matrix E_i(_row, 1);
    Matrix result(_row, 1);
    for (int i = 0; i < _row; ++i) {
        for (int j = 0; j < _col; ++j) {
            if (j == i) {
                E_i.matrix[j] = 1;
            } else {
                E_i.matrix[j] = 0;
            }
        }
        result.Gauss(*this, E_i);
        for (int j = 0; j < _col; ++j) {
            inversematrix.matrix[i * _col + j] = result.matrix[j];
        }
    }

    for (int i = 0; i < _row; ++i) {
        for (int j = 0; j < i; ++j) {
            double c = inversematrix.matrix[i * _col + j];
            inversematrix.matrix[i * _col + j] = inversematrix.matrix[j * _col + i];
            inversematrix.matrix[j * _col + i] = c;
        }
    }
    Matrix E(_row, _col);
    for (int i = 0; i < _row; ++i) {
        for (int j = 0; j < _col; ++j) {
            if (i == j) {
                E.matrix[i * _col + j] = 1;
            }
        }
    }
    assert(*this * inversematrix == E);

    return inversematrix;
}

double Matrix::determinant() {
    LU_decomposition();
    double det = 1;
    for (int i = 0; i < _row; ++i) {
        det *= _U_matrix[i * _col + i];
    }
    return det;
}

double Matrix::condition_number() {
    Matrix inv_m;
    inv_m = inverse();
    double cond = norm() * inv_m.norm();
    return cond;
}

/*bool Matrix::diagonal() {
    double sum;
    int k = 1;
    for (int i = 0; i < _row; ++i) {
        sum = 0;
        for (int j = 0; j < _col; ++j) {
            sum += fabs(matrix[i * _col + j]);
        }
        sum -= fabs(matrix[i * _col + i]);
        if (sum > matrix[i * _col + i]) {
          k = 0;
        }
    }
    return (k == 1);
}*/

bool Matrix::convergence(Matrix &x1, Matrix &x2) const {
    double norm = 0;
    for (int i = 0; i < _row; ++i) {
        norm += (x1.matrix[i] - x2.matrix[i]) * (x1.matrix[i] - x2.matrix[i]);
    }
    return (sqrt(norm) < 1e-7);
}

void Matrix::Seidel(Matrix &A, Matrix &matrix_b) {
    Matrix p(_row, 1);

    for (int i = 0; i < _row; ++i) {
        matrix[i] = 0;
    }
    int m = 0;
    do {
        for (int i = 0; i < _row; ++i) {
            p.matrix[i] = matrix[i];
        }
        for (int i = 0; i < A._row; ++i) {
            double temp = 0;
            for (int j = 0; j < A._col; ++j) {
                if (j != i) {
                    temp += A.matrix[i * A._col + j] * matrix[j];
                }
            }
            matrix[i] = (matrix_b.matrix[i] - temp) / A.matrix[i * A._col + i];
        }
        ++m;
    } while (!convergence(*this, p));
}

double diff(Matrix &A, Matrix &B) {
    if (A.get_row_size() != B.get_row_size())
        return -1;
    double d = 0;
    for (int i = 0; i < A.get_row_size(); ++i) {
        d += fabs(A.matrix[i] - B.matrix[i]);
    }
    return d;
}
#pragma once

#include <iostream>
#include <string>
#include <vector>

#define EPS 1e-13;

class Matrix {
public:
    Matrix(int m = 1, int n = 1, bool random = false);

    Matrix(std::string file_name);

    Matrix(Matrix &tmp);

    Matrix(int m, int n, double *array);

    ~Matrix();

    void matrix_to_file(std::string filename = "");

    int get_row_size();

    int get_col_size();

    Matrix operator*(const Matrix &right);

    Matrix &operator=(const Matrix &right);

    bool operator==(const Matrix &right);

    void Gauss(Matrix &A, Matrix &matrix_b);

    void Seidel(Matrix &A, Matrix &matrix_b);

    double condition_number();

    double *matrix;

private:
    int _row;
    int _col;
    double *_L_matrix = nullptr;
    double *_U_matrix = nullptr;

    double *create_matrix(int m = 1, int n = 1);

    void _get_determinant();

    void print_array(double *array, int row, int col) const;

    bool convergence(Matrix &x1, Matrix &x2) const;

    double max_eigenvalue();

    double norm() const;

    Matrix inverse();

    double determinant();

    void print_matrix();

    Matrix transpose() const;

    void print_LU_decomposition();

    double *get_L();

    double *get_U();

    void LU_decomposition();


};

double diff(Matrix &A, Matrix &B);
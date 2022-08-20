#include <iostream>
#include "Matrix.hpp"

int main() {

    {
        Matrix matrix_A("matrix_A.txt");
        Matrix matrix_b("matrix_b.txt");
        Matrix matrix_x(100, 1);
        std::cout << "Gauss method solution" << std::endl;
        matrix_x.Gauss(matrix_A, matrix_b);
        std::cout << "Output to file x.txt" << std::endl;
        matrix_x.matrix_to_file("x");

        std::cout << "Condition number: " << matrix_A.condition_number() << std::endl;

        Matrix matrix_z(100, 1);
        std::cout << "Seidel method solution" << std::endl;
        matrix_z.Seidel(matrix_A, matrix_b);
        matrix_z.matrix_to_file("z");
        std::cout << "Output to file z.txt" << std::endl;

        std::cout << "diff: " << diff(matrix_x, matrix_z) << std::endl;

    }

    return 0;
}

#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <ctime>

using namespace std;

typedef vector<vector<int>> Matrix;

Matrix generate_matrix(int rows, int cols, mt19937& rng) {
    uniform_int_distribution<int> dist(0, 10);
    Matrix mat(rows, vector<int>(cols));
    for (auto& row : mat)
        for (auto& val : row)
            val = dist(rng);
    return mat;
}

void write_matrix(ofstream& out, const Matrix& mat) {
    out << mat.size() << " " << mat[0].size() << "\n";
    for (const auto& row : mat) {
        for (int val : row)
            out << val << " ";
        out << "\n";
    }
}

int main() {
    int n = 1;

    mt19937 rng(time(nullptr));
    uniform_int_distribution<int> size_dist(5000, 5000);

    ofstream matrix_file("pair_rand_5000.txt");
    matrix_file << n << "\n";

    for (int i = 0; i < n; ++i) {
        int rows_A = size_dist(rng);
        int common_dim = size_dist(rng);
        int cols_B = size_dist(rng);

        Matrix A = generate_matrix(rows_A, common_dim, rng);
        Matrix B = generate_matrix(common_dim, cols_B, rng);

        write_matrix(matrix_file, A);
        write_matrix(matrix_file, B);
    }

    matrix_file.close();
    cout << "Matrix pairs saved to file\n";
    return 0;
}

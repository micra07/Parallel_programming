#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>

using namespace std;
using namespace chrono;

typedef vector<vector<int>> Matrix;

Matrix read_matrix(ifstream& in) {
    int rows, cols;
    in >> rows >> cols;
    Matrix mat(rows, vector<int>(cols));
    for (auto& row : mat)
        for (int& val : row)
            in >> val;
    return mat;
}



Matrix multiply_matrices(const Matrix& A, const Matrix& B) {
    int n = A.size();
    int m = B[0].size();
    int common = A[0].size();
    Matrix result(n, vector<int>(m, 0));
#pragma omp parallel for shared(A, B, result)
    for (int i = 0; i < n; ++i) {
        for (int k = 0; k < common; ++k) {
            for (int j = 0; j < m; ++j) {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return result;
}

void write_matrix(ofstream& out, const Matrix& mat) {
    out << mat.size() << " " << mat[0].size() << "\n";
    for (const auto& row : mat) {
        for (int val : row)
            out << val << " ";
        out << "\n";
    }
}

int main(int argc, char* argv[]) {
    if (argc < 3 || argc % 2 != 1) {
        cerr << "Usage: " << argv[0] << " 10_rand_pairs.txt 10_rand_pairs_result1.txt [100_rand_pairs.txt 100_rand_pairs_result.txt \
    1000_rand_pairs.txt 1000_rand_pairs_result.txt 10000_rand_pairs.txt 10000_rand_pairs_result.txt pair_rand_100.txt pair_rand_100_result.txt \
    pair_rand_500.txt pair_rand_500_result.txt pair_rand_1000.txt pair_rand_1000_result.txt pair_rand_2000.txt pair_rand_2000_result.txt \
    10_five_100_pairs.txt 10_five_100_pairs_result.txt 100_five_100_pairs.txt 100_five_100_pairs_result.txt 1000_five_100_pairs.txt 1000_five_100_pairs_result.txt]\n";
        return 1;
    }

    ofstream log_file("log.txt");
    if (!log_file) {
        cerr << "Could not open log.txt for writing.\n";
        return 1;
    }

    long long total_time_us = 0;
    int file_pairs = (argc - 1) / 2;

    for (int i = 1; i < argc; i += 2) {
        string matrix_filename = argv[i];
        string result_filename = argv[i + 1];

        ifstream matrix_file(matrix_filename);
        if (!matrix_file) {
            cerr << "Cannot open input file: " << matrix_filename << endl;
            continue;
        }

        ofstream result_file(result_filename);
        if (!result_file) {
            cerr << "Cannot open output file: " << result_filename << endl;
            continue;
        }

        int n;
        matrix_file >> n;

        cout << "Processing: " << matrix_filename << " (" << n << " pairs)\n";

        auto start = high_resolution_clock::now();
        for (int pair_id = 0; pair_id < n; ++pair_id) {
            Matrix A = read_matrix(matrix_file);
            Matrix B = read_matrix(matrix_file);
            Matrix C = multiply_matrices(A, B);
            write_matrix(result_file, C);
        }
        auto end = high_resolution_clock::now();

        long long duration_ms = duration_cast<milliseconds>(end - start).count();
        //long long duration_us = duration_cast<microseconds>(end - start).count();
        total_time_us += duration_ms;
        //total_time_us += duration_us;

        log_file << matrix_filename << " -> " << result_filename << ": "
            << duration_ms << " ms\n";
        //<< duration_us << " 탎\n";

        cout << "Finished: " << matrix_filename << " => Time: " << duration_ms << " ms\n";
        //cout << "Finished: " << matrix_filename << " => Time: " << duration_us << " 탎\n";

        matrix_file.close();
        result_file.close();
    }

    log_file << "Total: " << total_time_us << " ms\n";
    cout << "All files processed. Total time: " << total_time_us << " ms\n";
    //log_file << "Total: " << total_time_us << " 탎\n";
    //cout << "All files processed. Total time: " << total_time_us << " 탎\n";

    log_file.close();
    return 0;
}

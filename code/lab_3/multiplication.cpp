#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <mpi.h>

using namespace std;
using namespace chrono;

typedef vector<int> FlatMatrix;

void read_matrix_flat(ifstream& in, FlatMatrix& mat, int& rows, int& cols) {
    in >> rows >> cols;
    mat.resize(rows * cols);
    for (int i = 0; i < rows * cols; ++i)
        in >> mat[i];
}

void write_matrix_flat(ofstream& out, const FlatMatrix& mat, int rows, int cols) {
    out << rows << " " << cols << "\n";
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j)
            out << mat[i * cols + j] << " ";
        out << "\n";
    }
}

void transpose_matrix(FlatMatrix& mat, int rows, int cols) {
    FlatMatrix temp(cols * rows);
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            temp[j * rows + i] = mat[i * cols + j];
    mat = temp;
}

void multiply_local(const FlatMatrix& A_part, const FlatMatrix& B_T, FlatMatrix& C_part,
    int local_rows, int A_cols, int B_cols) {
    for (int i = 0; i < local_rows; ++i) {
        for (int j = 0; j < B_cols; ++j) {
            int sum = 0;
            for (int k = 0; k < A_cols; ++k) {
                sum += A_part[i * A_cols + k] * B_T[j * A_cols + k];
            }
            C_part[i * B_cols + j] = sum;
        }
    }
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    int rank, num_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    if ((argc < 3) || (argc % 2 != 1)) {
        if (rank == 0) {
            cerr << "Usage: " << argv[0] << " 10_rand_pairs.txt 10_rand_pairs_result1.txt [100_rand_pairs.txt 100_rand_pairs_result.txt "
                << "1000_rand_pairs.txt 1000_rand_pairs_result.txt 10000_rand_pairs.txt 10000_rand_pairs_result.txt pair_rand_100.txt pair_rand_100_result.txt "
                << "pair_rand_500.txt pair_rand_500_result.txt pair_rand_1000.txt pair_rand_1000_result.txt pair_rand_2000.txt pair_rand_2000_result.txt "
                << "10_five_100_pairs.txt 10_five_100_pairs_result.txt 100_five_100_pairs.txt 100_five_100_pairs_result.txt 1000_five_100_pairs.txt 1000_five_100_pairs_result.txt]\n";
        }
        MPI_Finalize();
        return 1;
    }

    ofstream log_file;
    long long total_time_ms = 0;

    if (rank == 0) {
        log_file.open("log.txt");
        if (!log_file) {
            cerr << "Could not open log.txt for writing.\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    for (int arg = 1; arg < argc; arg += 2) {
        string input_file = argv[arg];
        string output_file = argv[arg + 1];

        ifstream fin;
        ofstream fout;
        int pair_count = 0;
        long long file_time_ms = 0;

        if (rank == 0) {
            fin.open(input_file);
            if (!fin) {
                cerr << "Cannot open input file: " << input_file << endl;
                continue;
            }
            fout.open(output_file);
            if (!fout) {
                cerr << "Cannot open output file: " << output_file << endl;
                continue;
            }
            fin >> pair_count;
            cout << "Processing: " << input_file << " (" << pair_count << " pairs)\n";
        }

        auto file_start = high_resolution_clock::now();
        MPI_Bcast(&pair_count, 1, MPI_INT, 0, MPI_COMM_WORLD);

        for (int p = 0; p < pair_count; ++p) {
            int A_rows = 0, A_cols = 0, B_rows = 0, B_cols = 0;
            FlatMatrix A, B;

            if (rank == 0) {
                read_matrix_flat(fin, A, A_rows, A_cols);
                read_matrix_flat(fin, B, B_rows, B_cols);

                if (A_cols != B_rows) {
                    cerr << "Matrix dimensions mismatch: A cols (" << A_cols
                        << ") != B rows (" << B_rows << ")\n";
                    MPI_Abort(MPI_COMM_WORLD, 1);
                }
            }

            int dims[4] = { A_rows, A_cols, B_rows, B_cols };
            MPI_Bcast(dims, 4, MPI_INT, 0, MPI_COMM_WORLD);
            A_rows = dims[0]; A_cols = dims[1];
            B_rows = dims[2]; B_cols = dims[3];

            int local_rows = A_rows / num_procs;
            int remainder = A_rows % num_procs;
            if (rank < remainder) local_rows++;

            int offset = 0;
            for (int i = 0; i < rank; ++i) {
                offset += (A_rows / num_procs) + (i < remainder ? 1 : 0);
            }

            FlatMatrix A_part(local_rows * A_cols);
            FlatMatrix B_T(B_cols * B_rows);
            FlatMatrix C_part(local_rows * B_cols);

            if (rank == 0) {
                for (int i = 0; i < local_rows; ++i)
                    for (int j = 0; j < A_cols; ++j)
                        A_part[i * A_cols + j] = A[(offset + i) * A_cols + j];

                for (int dest = 1; dest < num_procs; ++dest) {
                    int dest_rows = A_rows / num_procs;
                    if (dest < remainder) dest_rows++;
                    int dest_offset = 0;
                    for (int i = 0; i < dest; ++i) {
                        dest_offset += (A_rows / num_procs) + (i < remainder ? 1 : 0);
                    }
                    MPI_Send(A.data() + dest_offset * A_cols, dest_rows * A_cols,
                        MPI_INT, dest, 0, MPI_COMM_WORLD);
                }

                transpose_matrix(B, B_rows, B_cols);
                B_T = B;
            }
            else {
                MPI_Recv(A_part.data(), local_rows * A_cols,
                    MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }

            MPI_Bcast(B_T.data(), B_rows * B_cols, MPI_INT, 0, MPI_COMM_WORLD);

            auto start = high_resolution_clock::now();
            multiply_local(A_part, B_T, C_part, local_rows, A_cols, B_cols);
            auto end = high_resolution_clock::now();

            file_time_ms += duration_cast<milliseconds>(end - start).count();

            FlatMatrix C;
            if (rank == 0)
                C.resize(A_rows * B_cols);

            if (rank == 0) {
                for (int i = 0; i < local_rows; ++i)
                    for (int j = 0; j < B_cols; ++j)
                        C[(offset + i) * B_cols + j] = C_part[i * B_cols + j];

                for (int src = 1; src < num_procs; ++src) {
                    int src_rows = A_rows / num_procs;
                    if (src < remainder) src_rows++;
                    int src_offset = 0;
                    for (int i = 0; i < src; ++i) {
                        src_offset += (A_rows / num_procs) + (i < remainder ? 1 : 0);
                    }
                    MPI_Recv(C.data() + src_offset * B_cols, src_rows * B_cols,
                        MPI_INT, src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
            else {
                MPI_Send(C_part.data(), local_rows * B_cols,
                    MPI_INT, 0, 0, MPI_COMM_WORLD);
            }

            if (rank == 0) {
                write_matrix_flat(fout, C, A_rows, B_cols);
            }
        }

        if (rank == 0) {
            auto file_end = high_resolution_clock::now();
            long long total_file_time = duration_cast<milliseconds>(file_end - file_start).count();
            total_time_ms += total_file_time;

            log_file << input_file << " -> " << output_file << ": " << total_file_time << " ms\n";
            cout << "Finished: " << input_file << " => Time: " << total_file_time << " ms\n";

            fin.close();
            fout.close();
        }
    }

    if (rank == 0) {
        log_file << "Total: " << total_time_ms << " ms\n";
        cout << "All files processed. Total time: " << total_time_ms << " ms\n";
        log_file.close();
    }

    MPI_Finalize();
    return 0;
}

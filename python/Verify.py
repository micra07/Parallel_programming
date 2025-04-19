import sys
print("Аргументы:", sys.argv)

def read_matrix(f):
    rows, cols = map(int, f.readline().split())
    matrix = [list(map(int, f.readline().split())) for _ in range(rows)]
    return matrix

def multiply(A, B):
    n, m, p = len(A), len(B[0]), len(B)
    result = [[0] * m for _ in range(n)]
    for i in range(n):
        for k in range(p):
            for j in range(m):
                result[i][j] += A[i][k] * B[k][j]
    return result

def matrices_equal(M1, M2):
    if len(M1) != len(M2) or len(M1[0]) != len(M2[0]):
        return False
    for row1, row2 in zip(M1, M2):
        if row1 != row2:
            return False
    return True

def check_file_pair(input_path, result_path):
    with open(input_path, 'r') as fin, open(result_path, 'r') as fres:
        n = int(fin.readline())
        for pair_idx in range(n):
            A = read_matrix(fin)
            B = read_matrix(fin)
            C_expected = multiply(A, B)
            C_actual = read_matrix(fres)
            if not matrices_equal(C_expected, C_actual):
                print(f"[❌] Ошибка в файле '{input_path}' / '{result_path}', пара #{pair_idx + 1}")
                return False
    print(f"[✅] Проверка пройдена: {input_path} / {result_path}")
    return True

def main():
    if len(sys.argv) < 3 or len(sys.argv) % 2 != 1:
        print("Usage: python verify_matrix_results.py 10_rand_pairs.txt 10_rand_pairs_result1.txt [100_rand_pairs.txt 100_rand_pairs_result.txt \
    1000_rand_pairs.txt 1000_rand_pairs_result.txt 10000_rand_pairs.txt 10000_rand_pairs_result.txt pair_rand_100.txt pair_rand_100_result.txt \
    pair_rand_500.txt pair_rand_500_result.txt pair_rand_1000.txt pair_rand_1000_result.txt pair_rand_2000.txt pair_rand_2000_result.txt \
    10_five_100_pairs.txt 10_five_100_pairs_result.txt 100_five_100_pairs.txt 100_five_100_pairs_result.txt 1000_five_100_pairs.txt 1000_five_100_pairs_result.txt]")
        sys.exit(1)

    all_correct = True
    for i in range(1, len(sys.argv), 2):
        input_file = sys.argv[i]
        result_file = sys.argv[i+1]
        correct = check_file_pair(input_file, result_file)
        if not correct:
            all_correct = False

    if all_correct:
        print("Все файлы проверены успешно ✅")
    else:
        print("Обнаружены ошибки ❌")

if __name__ == "__main__":
    main()

/**
 * Вариант задания 4.
 * Нахождение обратной матрицы для квадратной матрицы размера nxn.
 * Выполнила студентка группы БПИ 199 Вахитова Диана.
 */
#include <iostream>
#include <fstream>
#include <thread>
#include <omp.h>
#include <sstream>
#include <vector>

using namespace std;

// Вывод матрицы.
void PrintMatrix(vector<vector<float>>& ar, int n, int m)
{
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) cout << ar[i][j] << "  ";
        printf("\n");
    }
    return;
}

// Вывод обратной матрицы.
void PrintInverse(vector<vector<float>>& ar, int n, int m)
{
    for (int i = 0; i < n; ++i) {
        for (int j = n; j < m; ++j) printf("%.3f  ", ar[i][j]);
        printf("\n");
    }
    return;
}

// Нахождение обратной матрицы.
void InverseOfMatrix(vector<vector<float>>& matrix, int order) {
    float temp_;
    // Печать первоначальной матрицы.
    printf("=== Matrix ===\n");
    PrintMatrix(matrix, order, order);

    for (int i = 0; i < order; ++i) matrix[i].resize(2 * order);
    for (int i = 0; i < order; ++i) {

        for (int j = 0; j < 2 * order; ++j) {
            if (j == (i + order))
                matrix[i][j] = 1;
        }
    }
    for (int i = order - 1; i > 0; --i) {
        if (matrix[i - 1][0] < matrix[i][0]) {
            vector<float> temp = matrix[i];
            matrix[i] = matrix[i - 1];
            matrix[i - 1] = temp;
        }
    }

    printf("\n=== Augmented Matrix ===\n");
    PrintMatrix(matrix, order, order * 2);
    // Тут в происходит обратный ход Гаусса, я бы сказала это неудобно параллелить.
    for (int i = 0; i < order; i++) {
        for (int j = 0; j < order; j++) {
            if (j != i) {
                temp_ = matrix[j][i] / matrix[i][i];
                for (int k = 0; k < 2 * order; k++) {
                    matrix[j][k] -= matrix[i][k] * temp_;
                }
            }
        }
    }
    return;
}

// Метод делает из диагональной маттрицы единичную. Невероятная польза, это происходит в t раз быстрее. 
void threadFunc(size_t i, vector<vector<float>>* A) {
    vector<vector<float>>& matrix = *A;
    float temp_ = matrix[i][i];
    for (size_t j = 0; j < matrix.size() * 2; ++j) {
        matrix[i][j] = matrix[i][j] / temp_;
    }
    printf("\n thread %d calculated row number: %d", omp_get_thread_num(), i);
}

int main(int argsNumber, char** args) {
    // Обработка аргументов командной строки. 
    if (argsNumber != 3) {
        std::cout << "Wrong Console Data" << std::endl;
        std::exit(EXIT_FAILURE);
    }
    // Для оптимизации используется работа с файлами. 
    int N, t;
    std::ifstream in{ args[1] };
    in >> N >> t;
    std::vector<std::vector<float>> matrix(N);
    std::vector<float> temp(N);
    // Заполняем матрицу. 
    for (auto i = 0; i < N; ++i) {
        for (auto j = 0; j < N; ++j) {
            in >> temp[j];
        }
        matrix[i] = temp;
    }

    InverseOfMatrix(matrix, N);

    #pragma omp parallel for
    for (int i = 0; i < matrix.size(); ++i) {
        threadFunc((size_t)i, &matrix);
    }

    // Печатаем результат. 
    printf("\n=== Inversed Matrix ===\n");
    PrintInverse(matrix, N, N * 2);
}
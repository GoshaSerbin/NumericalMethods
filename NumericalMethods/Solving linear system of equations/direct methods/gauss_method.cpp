#include"header.h"

//Находит максимальный по модулю коэффициент в столбце k матрицы A
size_t max_coef_in_col(T** const A, size_t n, size_t k) {
	size_t i_max = k;
	T max = abs(A[k][k]);

	for (size_t i = k + 1; i < n; i++) {
		if (abs(A[i][k]) > max) {
			max = abs(A[i][k]);
			i_max = i;
		}
	}
	return i_max;
}

//Обратный ход метода Гаусса
void gauss_reverse(const size_t& n, T** const A, T* const b, T* x) {
	for (int i = n - 1; i >= 0; i--) {
		T sum = 0;
		for (size_t j = i + 1; j < n; j++) {
			sum += A[i][j] * x[j];
		}
		x[i] = (b[i] - sum) / A[i][i];
	}
}

//Метод Гаусса решения СЛАУ Ax=b, записывает решение в вектор x
int gauss_method(const size_t& n, T** const A, T* const b, T* x) {

	T** _A = new T * [n];
	for (size_t i = 0; i < n; i++) {
		_A[i] = new T[n];
		for (size_t j = 0; j < n; j++) {
			_A[i][j] = A[i][j];
		}
	}
	T* _b = new T[n];
	for (size_t i = 0; i < n; i++) {
		_b[i] = b[i];
	}
	//Прямой ход
	for (size_t k = 0; k < n - 1; k++) {
		size_t i_max = max_coef_in_col(A, n, k);
		if (k != i_max) {
			std::swap(_A[k], _A[i_max]);
			std::swap(_b[k], _b[i_max]);
		}

		if (abs(_A[k][k]) < EPS) {
			return -1;
		}

		for (size_t i = k + 1; i < n; i++) {
			if (abs(_A[i][k]) > EPS) {
				T c = _A[i][k] / _A[k][k];
				_A[i][k] = 0;
				for (size_t j = k + 1; j < n; j++) {
					_A[i][j] = _A[i][j] - c * _A[k][j];
				}
				_b[i] = _b[i] - c * _b[k];
			}

		}

	}

	if (abs(_A[n - 1][n - 1]) < EPS) {
		return -1;
	}
	//Обратный ход
	gauss_reverse(n, _A, _b, x);



	for (size_t i = 0; i < n; i++)	delete[] _A[i];
	delete[] _A;

	delete[] _b;

	return 0;
}



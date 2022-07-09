#include"header.h"


//Показать систему уравнений
void showLinearSystem(const size_t& n, T** const A, T* const b) {
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < n; j++) {
			std::cout << std::setw(15) << A[i][j] << "\t";
		}
		std::cout << "  |  " << b[i];
		std::cout << "\n                                                                  |  \n";
	}
}

//Вывести на экран матрицу A
void showMatrix(const size_t& n, T** const A) {
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < n; j++) {
			std::cout << std::setw(15) << A[i][j] << "\t";
		}
		std::cout << "\n\n";
	}
}

//Выводит произведение матриц
void showMatrixMultiplication(size_t n, T** const A, T** const B) {

	T** C = new T * [n];
	for (size_t i = 0; i < n; i++)	C[i] = new T[n];


	for (size_t i = 0; i < n; i++)
		for (size_t j = 0; j < n; j++) {
			T sum = 0;
			for (size_t k = 0; k < n; k++) {
				sum += A[i][k] * B[k][j];
			}
			C[i][j] = sum;
		}

	showMatrix(n, C);

	for (size_t i = 0; i < n; i++)	delete[] C[i];
	delete[] C;
}

//Норма вектора
T norm1(size_t n, T* x) {
	T sum = 0;
	for (size_t i = 0; i < n; i++) {
		sum += abs(x[i]);
	}

	return sum;
}

T norm1(size_t n, T** A) {

	T sum_max = 0;
	for (size_t j = 0; j < n; j++) {
		T sum = 0;
		for (size_t i = 0; i < n; i++) {
			sum += abs(A[i][j]);
		}
		if (sum > sum_max) sum_max = sum;
	}
	return sum_max;
}

T normInf(size_t n, T* x) {
	T x_max = 0;
	for (size_t i = 0; i < n; i++) {
		if (x_max > abs(x[i])) {
			x_max = abs(x[i]);
		}
	}

	return x_max;
}

T normInf(size_t n, T** A) {

	T sum_max = 0;
	for (size_t i = 0; i < n; i++) {
		T sum = 0;
		for (size_t j = 0; j < n; j++) {
			sum += abs(A[i][j]);
		}
		if (sum > sum_max) sum_max = sum;
	}
	return sum_max;
}


//Возвращает обратную матрицу к A
T** inverse(size_t n, T** const A) {

	T** invA = new T * [n];
	for (size_t i = 0; i < n; i++)	invA[i] = new T[n];

	T* x = new T[n];

	T** _A = new T * [n];
	for (size_t i = 0; i < n; i++) {
		_A[i] = new T[n];
		for (size_t j = 0; j < n; j++) {
			_A[i][j] = A[i][j];
		}
	}

	T** E = new T * [n];
	for (size_t i = 0; i < n; i++) {
		E[i] = new T[n];
		for (size_t j = 0; j < n; j++) {
			if (i == j) {
				E[i][j] = 1;
			}
			else E[i][j] = 0;
		}
	}

	T** Q = new T * [n];
	for (size_t i = 0; i < n; i++)	Q[i] = new T[n];

	T** R = new T * [n];
	for (size_t i = 0; i < n; i++)	R[i] = new T[n];


	QR_decomposition(n, A, Q, R);
	for (size_t j = 0; j < n; j++) {
		QR_solve(n, Q, R, E[j], x);

		for (size_t i = 0; i < n; i++) {
			invA[i][j] = x[i];
		}
	}


	for (size_t i = 0; i < n; i++)	delete[] E[i];
	delete[] E;

	for (size_t i = 0; i < n; i++)	delete[] _A[i];
	delete[] _A;

	for (size_t i = 0; i < n; i++)	delete[] Q[i];
	delete[] Q;

	for (size_t i = 0; i < n; i++)	delete[] R[i];
	delete[] R;

	delete[] x;

	return invA;
}

//Транспонирует матрицу
void transpose(size_t n, T** const A) {

	for (size_t i = 0; i < n; i++)
		for (size_t j = 0; j < i; j++) {
			T A_ij = A[i][j];
			A[i][j] = A[j][i];
			A[j][i] = A_ij;
		}

}

//внесение возмущений в вектор b
void disturbance(size_t n, T* b) {
	for (size_t i = 0; i < n; i++) {
		b[i] = b[i] * 1 + 2.0*disturbancePercent * (0.01 * (rand() % 101) - 0.5);
	}
}

//Норма вектора невязки Ax*-b
T residual(const size_t& n, T** const A, T* const b, T* x) {
	T* _A = new T[n];

	for (size_t i = 0; i < n; i++) {
		T sum = 0;
		for (size_t j = 0; j < n; j++) {
			sum += A[i][j] * x[j];
		}
		_A[i] = sum - b[i];
	}
	T r = norm1(n, _A);

	delete[] _A;
	return r;
}


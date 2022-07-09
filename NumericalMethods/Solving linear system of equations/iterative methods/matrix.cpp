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
void MatrixMultiplication(size_t n, T** const A, T** const B, T** C) {


	for (size_t i = 0; i < n; i++)
		for (size_t j = 0; j < n; j++) {
			T sum = 0;
			for (size_t k = 0; k < n; k++) {
				sum += A[i][k] * B[k][j];
			}
			C[i][j] = sum;
		}

}

//Норма вектора
//Возможно сделать апгрейд и добавить параметр "Inf", чтобы уменьшить количество функций
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
		if (x_max < abs(x[i])) {
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


//Разложение матрицы A на Q*R
void QR_decomposition(size_t n, T** A, T** Q, T** R) {

	for (size_t i = 0; i < n; i++)
		for (size_t j = 0; j < n; j++)
			if (i != j) Q[i][j] = 0;
			else	Q[i][j] = 1;

	for (size_t i = 0; i < n; i++)
		for (size_t j = 0; j < n; j++)
			R[i][j] = A[i][j];

	T c = 0, s = 0;
	for (size_t k = 0; k < n - 1; k++)
		for (size_t i = k + 1; i < n; i++)
			if (abs(R[i][k]) > EPS) {

				T sq = sqrt(pow(R[k][k], 2) + pow(R[i][k], 2));
				c = R[k][k] / sq;
				s = R[i][k] / sq;
				for (size_t j = k; j < n; j++) {
					T R_kj = R[k][j];
					R[k][j] = c * R[k][j] + s * R[i][j];
					R[i][j] = -s * R_kj + c * R[i][j];
					T Q_kj = Q[k][j];
					Q[k][j] = c * Q[k][j] + s * Q[i][j];
					Q[i][j] = -s * Q_kj + c * Q[i][j];

				}
				R[i][k] = 0;
			}

	transpose(n, Q);
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

//Решение СЛАУ Ax=b с помощью Q и R
void QR_solve(size_t n, T** Q, T** R, T* b, T* x) {


	transpose(n, Q);

	//Находим b*, транспонируя Q, тк она ортогональна
	T* _b = new T[n];
	for (size_t i = 0; i < n; i++)	_b[i] = b[i];

	for (size_t i = 0; i < n; i++) {
		T sum = 0;
		for (size_t k = 0; k < n; k++) {
			sum += b[k] * Q[i][k];
		}
		_b[i] = sum;
	}


	gauss_reverse(n, R, _b, x);

	transpose(n, Q);

	delete[] _b;

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

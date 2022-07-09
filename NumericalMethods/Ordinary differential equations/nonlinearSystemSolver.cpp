#include"header.h"

T norm1(T* x, size_t n) {
	T sum = 0;
	for (size_t i = 0; i < n; i++) {
		sum += abs(x[i]);
	}

	return sum;
}

T norm1(T** A, size_t n) {

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

T normInf(T* x, size_t n) {
	T x_max = 0;
	for (size_t i = 0; i < n; i++) {
		if (x_max < abs(x[i])) {
			x_max = abs(x[i]);
		}
	}

	return x_max;
}

T normInf(T** A, size_t n) {

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

//Транспонирует матрицу
void transpose(size_t n, T** const A) {

	for (size_t i = 0; i < n; i++)
		for (size_t j = 0; j < i; j++) {
			T A_ij = A[i][j];
			A[i][j] = A[j][i];
			A[j][i] = A_ij;
		}

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

//Решение СЛАУ Ax=b методом QR-разложения, решение записывается в вектор x
int QR_method(const size_t& n, T** const A, T* const b, T* x) {

	T** Q = new T * [n];
	for (size_t i = 0; i < n; i++)	Q[i] = new T[n];

	T** R = new T * [n];
	for (size_t i = 0; i < n; i++) 	R[i] = new T[n];

	QR_decomposition(n, A, Q, R);

	T det = 1;
	for (size_t i = 0; i < n; i++) {
		det = det * R[i][i];
	}
	if (abs(det) < EPS) {
		return -1;
	}

	QR_solve(n, Q, R, b, x);

	for (size_t i = 0; i < n; i++)		delete[] Q[i];
	delete[] Q;

	for (size_t i = 0; i < n; i++) 		delete[] R[i];
	delete[] R;

	return 0;
}

//нейлинейная система относительно y (y_{n+1})
T* systemForImplicitEuler(T* (*f)(T, T*, size_t, T*), T t, T* y_n, T tau, T* y, size_t N, T* vec) {
	T* vec2 = new T[N];

	for (size_t i = 0; i < N; i++) {
		vec[i] = (y[i] - y_n[i]) / tau - f(t, y, N, vec2)[i];
	}

	delete[] vec2;
	return vec;
}

T* systemForSymmetric(T* (*f)(T, T*, size_t, T*), T t_n, T* y_n, T tau, T* y, size_t N, T* vec) {
	T* vec2 = new T[N];
	T* vec3 = new T[N];
	for (size_t i = 0; i < N; i++) {
		vec[i] = (y[i] - y_n[i]) / tau - (f(t_n, y_n, N, vec2)[i] + f(t_n + tau, y, N, vec3)[i]) / 2;
	}

	delete[] vec2;
	delete[] vec3;
	return vec;
}

bool inRegion(T* Y, T* d1, T* d2, size_t N) {
	bool inReg = true;
	for (size_t i = 0; i < N; i++) {
		if (Y[i] > d2[i] || Y[i] < d1[i]) {
			inReg = false;
			break;
		}
	}
	return inReg;
}

//Решение нелинейной системы из предыдущей лабораторной
size_t NSolveNewton(T* (*system)(T* (*f)(T, T*, size_t, T*), T, T*, T, T*, size_t, T*), T* (*f)(T, T*, size_t, T*), T t, T* y_n, T tau, T* d1, T* d2, T* initY, size_t N, T* YNew) {


	(T * (*)(T, T*, size_t, T*))* system;
	T** JacobiMatrix = new T * [N];	for (size_t i = 0; i < N; i++) { JacobiMatrix[i] = new T[N]; }
	T* YLast = new T[N];
	for (size_t i = 0; i < N; i++) {
		YNew[i] = initY[i];
	}

	//Начальное приближение = XNew
	T* dY = new T[N];
	T eps = difference;//или dx
	T count = 1;


	T* vec = new T[N];
	T* vec2 = new T[N];
	if (norm1(system(f, t, y_n, tau, YNew, N, vec), N) >= EPSZERO) {
		do {
			count++;
			if (count > 30) {
				std::cout << "\nНе фортануло\n";
				return -1;
			}
			for (size_t i = 0; i < N; i++) {
				YLast[i] = YNew[i];
			}
			for (size_t i = 0; i < N; i++) {
				for (size_t j = 0; j < N; j++) {
					T* Yk = new T[N];
					for (size_t k = 0; k < N; k++) {
						Yk[k] = YLast[k];
						if (k == j) {
							Yk[k] += eps;
						}
					}
					JacobiMatrix[i][j] = (system(f, t, y_n, tau, Yk, N, vec)[i] - system(f, t, y_n, tau, YLast, N, vec2)[i]) / eps;//Вместо функции возвращающей вектор логичней бы было сделать вектор из функций
					delete[] Yk;
				}
			}
			T* b = new T[N];
			for (size_t i = 0; i < N; i++) {
				b[i] = -system(f, t, y_n, tau, YLast, N, vec)[i];
			}
			QR_method(N, JacobiMatrix, b, YNew);

			for (size_t i = 0; i < N; i++) {
				YNew[i] = YNew[i] + YLast[i];
			}
			delete[] b;


			//	std::cout << "cash2";
			for (size_t i = 0; i < 2; i++) {
				dY[i] = YNew[i] - YLast[i];
			}

			T counter = 0;
			while (!inRegion(YNew, d1, d2, N)) {
				counter++;
				if (counter > 20) break;
				for (size_t i = 0; i < N; i++) {
					YNew[i] = (YNew[i] + YLast[i]) / 2.;
				}
			}
			//	std::cout << "cash3";
		} while (norm1(dY, N) > EPSZERO);
	}
	else {
		for (size_t i = 0; i < N; i++) {
			YNew[i] = initY[i];
		}
	}

	//std::cout << "\nКоличество итераций: " << norm1(system(f, t, y_n, tau, YNew, N, vec), N)<<"\n";


	for (size_t i = 0; i < 2; i++) {
		delete[] JacobiMatrix[i];
	}
	delete[] JacobiMatrix;

	delete[] dY;
	delete[] YLast;

	delete[] vec;
	delete[] vec2;

	if (count > 30) return -1;
	return 0;
}

#include"header.h"


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
				c = R[k][k] /sq;
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
void QR_solve(size_t n,T** Q,T** R,T* b,T* x) {


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

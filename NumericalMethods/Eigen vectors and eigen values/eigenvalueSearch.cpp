#include"header.h"

extern size_t iterationCounter;

extern size_t multiOperCounter;


//рекурсивная функция, находящая по цепочке собственные числа
void eigenvalueN(size_t n, size_t size, T** const A, bool withShift, bool forHessenburg, T* eigenvalues) {

	if (size == 1) {
		eigenvalues[n - 1] = A[0][0];
	}
	else {


		T** Ak = new T * [size];
		for (size_t i = 0; i < size; i++)	Ak[i] = new T[size];

		for (size_t i = 0; i < size; i++) {
			for (size_t j = 0; j < size; j++) {
				Ak[i][j] = A[i][j];
			}
		}

		T** Q = new T * [size];
		for (size_t i = 0; i < size; i++)	Q[i] = new T[size];

		T** R = new T * [size];
		for (size_t i = 0; i < size; i++) 	R[i] = new T[size];

		T count = 0;
		while (norm1(size-1,Ak[size-1])>EPS && count<15) {
			count++;
			iterationCounter++;
			T shift = Ak[size - 1][size - 1];
			if (withShift && size < n + 1) {
				for (size_t i = 0; i < size; i++) {
					Ak[i][i] -= shift;
				}
			//	std::cout << shift << "\n";
			}

			iterationCounter++;

			
			if (forHessenburg) {
				QR_decForHessenburg(size, Ak, Q, R);

			}
			else {
				QR_decomposition(size, Ak, Q, R);
				//std::cout << "Q^\n";
				//showMatrix(size, Q);
				//std::cout << "R^\n";
				//showMatrix(size, R);
			}


			for (size_t i = 0; i < size; i++) {
				for (size_t j = 0; j < size; j++) {
					if (forHessenburg && i > j + 1) {
						Ak[i][j] = 0;
					}
					else {
						T sum1 = 0;
						for (size_t k = 0; k < size; k++) {
							sum1 += R[i][k] * Q[k][j];
							multiOperCounter += 1;
						}
						Ak[i][j] = sum1;
					}

				}
			}

			if (withShift && size < n+1) {
				for (size_t i = 0; i < size; i++) {
					Ak[i][i] += shift;
				}
			}
			//std::cout << "Матрица Ak:\n";
			//showMatrix(size, Ak);
			//std::cout << "\n";


		}
		//std::cout << "Матрица Ak:\n";
		//showMatrix(size, Ak);
		//std::cout << "\n";
		for (size_t i = 0; i < size; i++)		delete[] Q[i];
		delete[] Q;

		for (size_t i = 0; i < size; i++) 		delete[] R[i];
		delete[] R;

		eigenvalues[n - size] = Ak[size-1][size-1];
		eigenvalueN(n, size-1, Ak, withShift,forHessenburg, eigenvalues);




		for (size_t i = 0; i < size; i++) 		delete[] Ak[i];
		delete[] Ak;
	}



}

void eigenvaluesSearch(size_t n, T** const A, bool withShift, bool forHessenburg, T* eigenvalues) {//QR method

	if (forHessenburg) {
		T** H = new T * [n];
		for (size_t i = 0; i < n; i++)	H[i] = new T[n];

		toHessenbergMatrix(n, A, H);


		eigenvalueN(n, n, H, withShift, forHessenburg, eigenvalues);

		for (size_t i = 0; i < n; i++) 		delete[] H[i];
		delete[] H;
	}
	else {
		eigenvalueN(n, n, A, withShift, forHessenburg, eigenvalues);
	}



}

void eigenvectorsSearch(size_t n, T** const A, T* eigenvalues,T** const eigenvectors) {//reverse iteration method

	T** Q = new T * [n];
	for (size_t i = 0; i < n; i++)	Q[i] = new T[n];

	T** R = new T * [n];
	for (size_t i = 0; i < n; i++) 	R[i] = new T[n];

	T** B = new T * [n];
	for (size_t i = 0; i < n; i++) 	B[i] = new T[n];


	for (size_t s = 0; s < n; s++) {

		for (size_t i = 0; i < n; i++) {
			for (size_t j = 0; j < n; j++) {
				if (i == j) B[i][j] = A[i][j] - eigenvalues[s];
				else	B[i][j] = A[i][j];
			}
		}

		

		T* x = new T[n];//начальное приближение
		for (size_t i = 0; i < n; i++) x[i] = 1/sqrt(n);
	//	x[0] = -0.1269; x[1] = 0.0428; x[2] = -0.03045; x[3] = 0.9905;
		T* y = new T[n];
		for (size_t i = 0; i < n; i++) y[i] = 0;

		QR_decomposition(n, B,Q,R);

		T count = 0;
		while (count<10 ) {
			count++;
			QR_solve(n, Q, R, x, y);
			T normY = normF(n, y);
			for (size_t i = 0; i < n; i++) x[i] = y[i]/normY;
			//for (size_t i = 0; i < n; i++) {
			//	std::cout << std::setw(15) << y[i] << '\t';
			//}
			//std::cout << "\n";
		}

		for (size_t i = 0; i < n; i++) {
			eigenvectors[s][i] = x[i];
		}

		delete[] x;
		delete[] y;
	}


	for (size_t i = 0; i < n; i++)		delete[] Q[i];
	delete[] Q;

	for (size_t i = 0; i < n; i++) 		delete[] R[i];
	delete[] R;

	for (size_t i = 0; i < n; i++) 		delete[] B[i];
	delete[] B;


}

void toHessenbergMatrix(size_t n, T** const A, T** const H) {


	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < n; j++) {
			H[i][j] = A[i][j];
		}
	}


	T c = 0, s = 0;

	for (size_t t = 0; t < n - 2; t++)//t=k-1
		for (size_t l = t + 2; l < n; l++)
			if (abs(H[l][t]) > EPSZERO) {
				size_t k = t + 1;
				T sq = sqrt(pow(H[k][k-1], 2) + pow(H[l][k-1], 2));
				c = H[k][k-1] / sq;
				s = H[l][k-1] / sq;

			//	multiOperCounter += 5;

				for (size_t j = k-1; j < n; j++) {
					T H_kj = H[k][j];
					H[k][j] = c * H[k][j] + s * H[l][j];
					H[l][j] = -s * H_kj + c * H[l][j];
				//	multiOperCounter += 4;
				}
				for (size_t j = k-1; j < n; j++) {
					T H_jk = H[j][k];
					H[j][k] = c * H[j][k] + s * H[j][l];
					H[j][l] = -s * H_jk + c * H[j][l];
				//	multiOperCounter += 4;
				}


				H[l][t] = 0;

			}
			else {
				H[l][t] = 0;
			}

	//std::cout << "H:\n";
	//showMatrix(n, H);

}

void eigenSearch(size_t n, T** const A, T* eigenvalues, T** eigenvectors) {


	T** Q = new T * [n];
	for (size_t i = 0; i < n; i++)	Q[i] = new T[n];

	T** R = new T * [n];
	for (size_t i = 0; i < n; i++) 	R[i] = new T[n];

	T** B = new T * [n];
	for (size_t i = 0; i < n; i++) 	B[i] = new T[n];


	for (size_t s = 0; s < n; s++) {

		
		T* x = new T[n];//начальное приближение
		for (size_t i = 0; i < n; i++) x[i] = eigenvectors[s][i];
	//	x[0] = -0.1269; x[1] = 0.0428; x[2] = -0.03045; x[3] = 0.9905;
		T normX = normF(n, x);
		for (size_t i = 0; i < n; i++) x[i] /= normX;
		T* y = new T[n];
		for (size_t i = 0; i < n; i++) y[i] = 0;

		T count = 0;
		T aproxLambda;
		while (count < 1) {
			count++;
			iterationCounter = 0;
			T sum;
			T* Ax = new T[n];
			for (size_t i = 0; i < n; i++) {	
				sum = 0;
				for (size_t j = 0; j < n; j++) {
					sum += A[i][j] * x[j];
				}
				Ax[i] = sum;
			}
			//for (size_t i = 0; i < n; i++) std::cout << Ax[i] << '\t';
			sum = 0;
			for (size_t i = 0; i < n; i++) {
				sum += Ax[i] * x[i];
			}
			aproxLambda = sum;
			//std::cout << "sum=" << sum << "\n";

			for (size_t i = 0; i < n; i++) {
				for (size_t j = 0; j < n; j++) {
					if (i == j) B[i][j] = A[i][j] - aproxLambda;
					else	B[i][j] = A[i][j];
				}
			}



			QR_decomposition(n, B, Q, R);

			QR_solve(n, Q, R, x, y);
			T normY = normF(n, y);
			for (size_t i = 0; i < n; i++) x[i] = y[i] / normY;
			//for (size_t i = 0; i < n; i++) {
			//	std::cout << std::setw(15) << y[i] << '\t';
			//}
			//std::cout << "\n";
		}

		for (size_t i = 0; i < n; i++) {
			eigenvectors[s][i] = x[i];
		}
		eigenvalues[s] = aproxLambda;


		delete[] x;
		delete[] y;
	}


	for (size_t i = 0; i < n; i++)		delete[] Q[i];
	delete[] Q;

	for (size_t i = 0; i < n; i++) 		delete[] R[i];
	delete[] R;

	for (size_t i = 0; i < n; i++) 		delete[] B[i];
	delete[] B;

}
#include"Header.h"

void initialApproximation(const size_t& n, T* x, T* const b) {
	for (size_t i = 0; i < n; i++) {
		x[i] = b[i];
	}


}

int simpleIterationMethod(const size_t& n, T** const A, T* const b, T* x, T tau, STOP_CRITERION sc, T* trueSol) {
	
	initialApproximation(n, x,b);

	T** C = new T * [n];//C=-(A-E)
	for (size_t i = 0; i < n; i++)	C[i] = new T[n];
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < n; j++) {
			if (i != j) {
				C[i][j] = - tau*A[i][j];
			}
			else {
				C[i][j] = -(tau*A[i][j] - 1);
			}
		}
	}

	T normC = normInf(n, C);
	std::cout << "         Метод простой итерации:";
	std::cout << "\nКубическая норма матрицы С: ||C||inf = " << normInf(n, C);
	std::cout << "\nОктаэдрическая норма матрицы С: ||C||1 = " << norm1(n, C);

	//if (normC >= 1) {
	//	std::cout << "\nATTENTION";
	//	return -1;
	//}



	T* x_new = new T[n];
	for (size_t i = 0; i < n; i++) {
		x_new[i] = x[i];
	}


	size_t k = 0;
	T* dx = new T[n];
	bool STOP = false;
	do{
		
		for (size_t i = 0; i < n; i++) {
			x[i] = x_new[i];
		}

		for (size_t i = 0; i < n; i++) {
			x_new[i] = b[i];
			for (size_t j = 0; j < n; j++) {
					x_new[i] -= A[i][j] * x[j];			
			}
			x_new[i] *= tau;
			x_new[i]+=x[i];
		}
		k++;
		for (size_t i = 0; i < n; i++) {
			dx[i] = x_new[i] - x[i];
		}

		switch (sc) {


		case CRITERION_NORM:
			if (normInf(n, dx) < (1 - normC) / normC * EPS) STOP = true;

			if (normC < 1) break;


		case CRITERION_EPS1:
			if (normInf(n, dx) < EPS) STOP = true;
			break;


		case CRITERION_TRUESOLVE:
			for (size_t i = 0; i < n; i++) {
				dx[i] = x_new[i] - trueSol[i];
			}
			if (normInf(n, dx) < EPS) STOP = true;
			break;

		case CRITERION_RESIDUAL:
			if (residual(n,A,b,x_new) < EPS) STOP = true;
			break;

		}

	} while (!STOP);


	std::cout << "\nКоличество итераций: " << k;




	for (size_t i = 0; i < n; i++) {
		x[i] = x_new[i];
	}
	delete[] x_new;

	for (size_t i = 0; i < n; i++) {
		delete[] C[i];
	}
	delete[] C;


	return 0;
}

int JacobiMethod(const size_t& n, T** const A, T* const b, T* x, STOP_CRITERION sc, T* trueSol) {

	initialApproximation(n, x,b);

	T** C = new T * [n];//C=-D^-1 * (L+U), где A=D+L+U
	for (size_t i = 0; i < n; i++)	C[i] = new T[n];
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < n; j++) {
			if (i != j) {
				C[i][j] = -A[i][j]/A[i][i];
			}
			else {
				C[i][j] = 0;
			}
		}
	}
	T normC = normInf(n, C);


	std::cout << "\n\n              Метод Якоби:";
	std::cout << "\nКубическая норма матрицы С: ||C||inf = " << normInf(n, C);
	std::cout << "\nОктаэдрическая норма матрицы С: ||C||1 = " << norm1(n, C);

	//if (normC >= 1) {
	//	std::cout << "\nATTENTION";
	//		return -1;
	//}

	T* x_new = new T[n];
	for (size_t i = 0; i < n; i++) {
		x_new[i] = x[i];
	}


	bool STOP = false;
	T* dx = new T[n];
	size_t k = 0;
	do{

		for (size_t i = 0; i < n; i++) {
			x[i] = x_new[i];
		}

		for (size_t i = 0; i < n; i++) {
			T Ax_i = 0;
			for (size_t j = 0; j < n; j++) {
				if(j!=i)	Ax_i += A[i][j] * x[j];
			}
			x_new[i] = (b[i] - Ax_i)/A[i][i];
		}
		k++;

		for (size_t i = 0; i < n; i++) {
			dx[i] = x_new[i] - x[i];
		}


		switch (sc) {

		case CRITERION_NORM:
			if (normInf(n, dx) < (1 - normC) / normC * EPS) STOP = true;
			if(normC<1) break;
		case CRITERION_EPS1:
			if (norm1(n, dx) < EPS) STOP = true;
			break;

		case CRITERION_TRUESOLVE:
			for (size_t i = 0; i < n; i++) {
				dx[i] = x_new[i] - trueSol[i];
			}
			if (norm1(n, dx) < EPS) STOP = true;
			break;

		case CRITERION_RESIDUAL:
			if (residual(n, A, b, x_new) < EPS) STOP = true;
			break;

		}

	} while (!STOP);
	std::cout << "\nКоличество итераций: " << k;

	for (size_t i = 0; i < n; i++) {
		delete[] C[i];
	}
	delete[] C;

	for (size_t i = 0; i < n; i++) {
		x[i] = x_new[i];
	}
	delete[] x_new;
	delete[] dx;


	return 0;
}


int SeidelsMethod(const size_t& n, T** const A, T* const b, T* x, STOP_CRITERION sc, T* trueSol) {

	initialApproximation(n, x,b);




	T** B1 = new T * [n];// A=D+L+U
	//C=(E+wD^-1 L)^-1 (-D^-1 U)
	for (size_t i = 0; i < n; i++)	B1[i] = new T[n];
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < n; j++) {
			if (i > j) {
				B1[i][j] = -A[i][j]/A[i][i];
			}
			else {
				B1[i][j] = 0;			
			}
		}
	}


	T** B2 = new T * [n];
	for (size_t i = 0; i < n; i++)	B2[i] = new T[n];
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < n; j++) {
			if (i < j) {
				B2[i][j] = -A[i][j] / A[i][i];
			}
			else {
				B2[i][j] = 0;
			}
		}
	}

	T** B = new T * [n];
	for (size_t i = 0; i < n; i++)	B[i] = new T[n];
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < n; j++) {

				B[i][j] = B1[i][j] +B2[i][j];
		}
	}

	T** M1 = new T * [n];// A=D+L+U
	//C=(E+wD^-1 L)^-1 (-D^-1 U)
	for (size_t i = 0; i < n; i++)	M1[i] = new T[n];
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < n; j++) {
			if (i >= j) {
				M1[i][j] = A[i][j];
			}
			else {
				M1[i][j] = 0;
			}
		}
	}


	T** M2 = new T * [n];
	for (size_t i = 0; i < n; i++)	M2[i] = new T[n];

	M2 = inverse(n, M1);

	T** C = new T * [n];
	for (size_t i = 0; i < n; i++)	C[i] = new T[n];
	MatrixMultiplication(n, M2, A, C);


	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < n; j++) {

			if (i != j) {
				C[i][j] = -C[i][j];
			}
			else {
				C[i][j] = -C[i][j] + 1;
			}
		}
	}

	for (size_t i = 0; i < n; i++) {
		delete[] M1[i];
	}
	delete[] M1;

	for (size_t i = 0; i < n; i++) {
		delete[] M2[i];
	}
	delete[] M2;

	T normB1 = normInf(n, B1);
	T normB2 = normInf(n, B2);


	for (size_t i = 0; i < n; i++) {
		delete[] B1[i];
	}
	delete[] B1;

	for (size_t i = 0; i < n; i++) {
		delete[] B2[i];
	}
	delete[] B2;


	if (normB1 + normB2 >= 1) {
		std::cout << "\nATTENTION";
		std::cout << "\n";
	//		return -1;
	}
	std::cout << "\n\n              Метод Зейделя:";
	std::cout << "\nКубическая норма матрицы С: ||C||inf = " << normInf(n, C);
	std::cout << "\nОктаэдрическая норма матрицы С: ||C||1 = " << norm1(n, C);
	std::cout << "\n ||B1||1+||B2||1 = " << normB1 + normB2;
	std::cout <<"\n"<< normB1 << "\t"<< (1 - normInf(n, B)) <<"\n";
	T* x_new = new T[n];
	for (size_t i = 0; i < n; i++) {
		x_new[i] = x[i];
	}

	bool STOP = false;
	T* dx = new T[n];
	size_t k = 0;
	do {

		for (size_t i = 0; i < n; i++) {
			x[i] = x_new[i];
		}

		T Ax_i = 0;//x[0]
		for (size_t j = 1; j < n; j++) {
			Ax_i += A[0][j] * x[j] / A[0][0];
		}
		x_new[0] = -Ax_i + b[0] / A[0][0];

		Ax_i = 0;//x[1]
		for (size_t j = 2; j < n; j++) {
			Ax_i += A[1][j] * x[j] / A[1][1];
		}
		x_new[1] = -A[1][0]*x_new[0]/A[1][1] - Ax_i + b[1] / A[1][1];

		for (size_t i = 2; i < n-1; i++) {

			T sum1 = 0;
			T sum2 = 0;
			for (size_t j = 0; j <= i-1; j++) {
				sum1 += A[i][j] * x_new[j] / A[i][i];
			}

			for (size_t j = i+1; j < n; j++) {
				sum2 += A[i][j] * x[j] / A[i][i];
			}

			x_new[i] = -sum1 - sum2 + b[i] / A[i][i];
	//		if (i == 1) { std::cout << "hey" << x_new[i] << "hey"; }
		}

		Ax_i = 0;//x[n-1]
		for (size_t j = 0; j < n-1; j++) {
			Ax_i += A[n-1][j] * x_new[j] / A[n-1][n-1];
		}
		x_new[n-1] = - Ax_i + b[n-1] / A[n-1][n-1];

		k++;

		for (size_t i = 0; i < n; i++) {
			dx[i] = x_new[i] - x[i];
		}
		switch (sc) {
		case CRITERION_NORM:
			//if (norm1(n, dx) < (1 - normC) / normC * EPS) STOP = true;			
			if (normInf(n, dx) < (1 - normInf(n,B)) / normB2*EPS) STOP = true;
			break;
		case CRITERION_EPS1:
			if (norm1(n, dx) < EPS) STOP = true;
			break;

		case CRITERION_TRUESOLVE:
			for (size_t i = 0; i < n; i++) {
				dx[i] = x_new[i] - trueSol[i];
			}
			if (norm1(n, dx) < EPS) STOP = true;
			break;


		case CRITERION_RESIDUAL:
			if (residual(n, A, b, x_new) < EPS) STOP = true;
			break;

		}


	} while (!STOP);

	std::cout << "\nКоличество итераций: " << k;


	for (size_t i = 0; i < n; i++) {
		x[i] = x_new[i];
	}



	for (size_t i = 0; i < n; i++) {
		delete[] C[i];
	}
	delete[] C;


	for (size_t i = 0; i < n; i++) {
		delete[] B[i];
	}
	delete[] B;

	delete[] dx;
	delete[] x_new;

	return 0;
}


int relaxationMethod(const size_t& n, T** const A, T* const b, T* x, T omega, STOP_CRITERION sc, T* trueSol) {

	initialApproximation(n, x,b);


	std::cout << "\n\n              Метод релаксации при w = :"<<omega;
//	std::cout << "\nКубическая норма матрицы С: ||C||inf = " << normInf(n, C);
//	std::cout << "\nОктаэдрическая норма матрицы С: ||C||1 = " << norm1(n, C);

	T* x_new = new T[n];
	for (size_t i = 0; i < n; i++) {
		x_new[i] = x[i];
	}
	bool STOP = false;
	T* dx = new T[n];
	size_t k = 0;
	do {

		for (size_t i = 0; i < n; i++) {
			x[i] = x_new[i];
		}

		for (size_t i = 0; i < n; i++) {

			T sum1 = 0;
			T sum2 = 0;
			for (size_t j = 0; j < i-1; j++) {
				sum1 += A[i][j] * x_new[j]/A[i][i];
			}
			for (size_t j = i+1; j < n; j++) {
				sum2 += A[i][j] * x_new[j]/A[i][i];
			}

			x_new[i] = -omega * sum1 + (1 - omega) * x[i] - omega * sum2 + omega * b[i] / A[i][i];
		}
		k++;

		for (size_t i = 0; i < n; i++) {
			dx[i] = x_new[i] - x[i];
		}

		switch (sc) {
		case CRITERION_NORM:
			//if (norm1(n, dx) < (1 - normC) / normC * EPS) STOP = true;


		case CRITERION_EPS1:
			if (norm1(n, dx) < EPS) STOP = true;
			break;

		case CRITERION_TRUESOLVE:
			for (size_t i = 0; i < n; i++) {
				dx[i] = x_new[i] - trueSol[i];
			}
			if (norm1(n, dx) < EPS) STOP = true;
			break;


		case CRITERION_RESIDUAL:
			if (residual(n, A, b, x_new) < EPS) STOP = true;
			break;
		}
	} while (!STOP);


	std::cout << "\nКоличество итераций: " << k;

	for (size_t i = 0; i < n; i++) {
		x[i] = x_new[i];
	}


	delete[] x_new;
	delete[] dx;
	return 0;
}

int relaxationMethod_3d(const size_t& n, T* const a, T* const b, T* const c, T* const d, T* x, T omega, STOP_CRITERION sc, T* trueSol) {


	initialApproximation(n, x, d);


	T** B1 = new T * [n];// A=D+L+U
	//C=(E+wD^-1 L)^-1 (-D^-1 U)
	for (size_t i = 0; i < n; i++)	B1[i] = new T[n];
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < n; j++) {
			if (i == j+1) {
				B1[i][j] = -a[i] / b[i];
			}
			else {
				B1[i][j] = 0;
			}
		}
	}


	T** B2 = new T * [n];
	for (size_t i = 0; i < n; i++)	B2[i] = new T[n];
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < n; j++) {
			if (i == j-1) {
				B2[i][j] = -c[i] / b[i];
			}
			else {
				B2[i][j] = 0;
			}
		}
	}

	T** B = new T * [n];
	for (size_t i = 0; i < n; i++)	B[i] = new T[n];
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < n; j++) {

			B[i][j] = B1[i][j] + B2[i][j];
		}
	}
	T normB = normInf(n, B);

	T normB1 = normInf(n, B1);

	T normB2 = normInf(n,B2);

	std::cout << "\n ||B1|| = " << normB1;
	std::cout << "\n ||B2|| = " << normB2;
	std::cout << "\n ||B1||+||B2|| = " << normB1+normB2 << "\n";


	for (size_t i = 0; i < n; i++) {
		delete[] B1[i];
	}
	delete[] B1;

	for (size_t i = 0; i < n; i++) {
		delete[] B2[i];
	}
	delete[] B2;

	T* x_new = new T[n];
	for (size_t i = 0; i < n; i++) {
		x_new[i] = x[i];
	}
	bool STOP = false;
	T* dx = new T[n];
	size_t k = 0;
	do{

		for (size_t i = 0; i < n; i++) {
			x[i] = x_new[i];
		}

		x_new[0] = (1 - omega) * x[0] - omega * c[0] * x[1] / b[0] + omega * d[0] / b[0];
		x_new[1] = -omega * a[1] * x_new[0] / b[1] + (1 - omega) * x[1] - omega * c[1] * x[2] / b[1] + omega * d[1] / b[1];

		for (size_t i = 2; i < n - 1; i++) {
			x_new[i] = -omega * a[i] * x_new[i - 1] / b[i] + (1 - omega) * x[i] - omega * c[i] * x[i + 1] / b[i] + omega * d[i] / b[i];

		}

		x_new[n - 1] = -omega * a[n - 1] * x_new[n - 2] / b[n - 1] + (1 - omega) * x[n - 1] + omega * d[n - 1] / b[n - 1];

		k++;

		for (size_t i = 0; i < n; i++) {
			dx[i] = x_new[i] - x[i];
		}

		switch (sc) {
		case CRITERION_NORM:
			if (normInf(n, dx) < (1 - normB) / normB2 * EPS) STOP = true;
			break;

		case CRITERION_TRUESOLVE:
			for (size_t i = 0; i < n; i++) {
				dx[i] = x_new[i] - trueSol[i];
			}
			if (norm1(n, dx) < EPS) STOP = true;
			break;

		case CRITERION_EPS1:
			if (norm1(n, dx) < EPS) STOP = true;
			break;


		case CRITERION_RESIDUAL:
			for (size_t i = 1; i < n-1; i++) {
				dx[i] = a[i] * x_new[i - 1] + b[i] * x_new[i] + c[i] * x_new[i + 1] - d[i];
				dx[0] = b[0]*x_new[0] + c[0]*x_new[1] - d[0];
				dx[n-1] = a[n-1] * x_new[n-2] + b[n-1] * x_new[n-1] - d[n-1];

			}
			if (normInf(n, dx) < EPS) STOP = true;
			break;

		}
	} while (!STOP);

	std::cout << "\nНевязка решения: " << normInf(n, dx);
	for (size_t i = 0; i < n; i++) {
		dx[i] = x_new[i] - x[i];
	}
	for (size_t i = 0; i < n; i++) { dx[i] = x[i] - trueSol[i]; }
	std::cout << "\nНорма вектора ошибки: " << norm1(n, dx) << "\n";

	for (size_t i = 0; i < n; i++) {
		x[i] = x_new[i];
	}

	std::cout << "\nКоличество итераций: " << k<<"\n";

	delete[] x_new;
	delete[] dx;
	return 0;


}

//"=Это прогонка случайно написана, жалко удалять
//
//std::cout << "\n\n              Метод релаксации для трехдиагональной матрицы при w = :" << omega;
//
//size_t k = 0;
//T* alpha = new T[n];
//T* beta = new T[n];
//alpha[1] = c[0] / b[0];
//beta[1] = d[0] / b[0];
//for (size_t i = 1; i < n - 1; i++) {
//	T z = (b[i] - a[i] * alpha[i]);
//	alpha[i + 1] = c[i] / z;
//	beta[i + 1] = (d[i] + a[i] * beta[i]) / z;
//}
//
//x[n - 1] = (alpha[n - 1] + a[n - 1] * beta[n - 1]) / (b[n - 1] - a[n - 1] * alpha[n - 1]);
//for (size_t i = n - 2; i > 0; i--) {
//	x[i] = alpha[i + 1] * x[i + 1] + beta[i + 1];
//}
//
//
//std::cout << "\nКоличество итераций: " << k;
//
//
//delete[] alpha;
//delete[] beta;
//return 0;
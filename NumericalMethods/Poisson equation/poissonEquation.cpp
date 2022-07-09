#include"header.h"

//решение СЛАУ методом прогонки
void solve_3d(size_t n, T* a, T* b, T* c, T* d, T* x) {

	T* alpha = new T[n + 1];
	T* beta = new T[n + 1];
	alpha[1] = -c[0] * 1.0 / b[0];
	beta[1] = d[0] * 1.0 / b[0];
	for (size_t i = 1; i < n - 1; i++) {
		T z = 1. * (b[i] + a[i] * alpha[i]);
		alpha[i + 1] = -c[i] / z;
		beta[i + 1] = (d[i] - a[i] * beta[i]) / z;
	}

	x[n - 1] = (d[n - 1] - a[n - 1] * beta[n - 1]) / (b[n - 1] + a[n - 1] * alpha[n - 1]);
	for (size_t i = n - 2; i > 0; i--) {
		x[i] = alpha[i + 1] * x[i + 1] + beta[i + 1];
	}
	x[0] = alpha[1] * x[1] + beta[1];
	delete[] alpha;
	delete[] beta;
}

//Норма вектора
T norm1(size_t n, T* x) {
	T sum = 0;
	for (size_t i = 0; i < n; i++) {
		sum += abs(x[i]);
	}

	return sum;
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

void DSolvePoissonEquation(T L1, T L2, T(*f)(T, T),
							BOUNDARY_CONDITIONS bc1, T(*G1)(T),
							BOUNDARY_CONDITIONS bc2, T(*G2)(T),
							BOUNDARY_CONDITIONS bc3, T(*G3)(T),
							BOUNDARY_CONDITIONS bc4, T(*G4)(T)) {

	T h1 = 5e-2, h2 = 5e-2;					//шаги сетки
	T tmax = 50, tau = h1 * h2;
	size_t n1 = size_t(round(L1 / h1)) + 1;	//количество узлов
	size_t n2 = size_t(round(L2 / h2)) + 1;	
	size_t nt = size_t(round(tmax / tau)) + 1;

	T** phi = new T*[n1];
	for (size_t i = 0; i < n1; i++) {
		phi[i] = new T[n2];
		for (size_t j = 0; j < n2; j++) {
			phi[i][j] = f(i * h1, j * h2);
		}
	}


	T** y1 = new T * [n1]; //y^{k}
	for (size_t i = 0; i < n1; i++) {
		y1[i] = new T[n2];
	}

	T** ymid = new T * [n1]; //y^{k+1/2}
	for (size_t i = 0; i < n1; i++) {
		ymid[i] = new T[n2];
	}

	T** y2 = new T * [n1]; //y^{k+1}
	for (size_t i = 0; i < n1; i++) {
		y2[i] = new T[n2];
	}

	std::ofstream out;				// поток для записи
	out.open("D:\\УРАВНЕНИЕ ПУАССОНА.txt");



	//начальные данные U0
	for (size_t i = 0; i < n1; i++) {
		for (size_t j = 0; j < n2; j++) {
			y1[i][j] = 1;
		}
	}


	T* A1 = new T[n1];
	T* B1 = new T[n1];
	T* C1 = new T[n1];
	T* D1 = new T[n1];

	T* A2 = new T[n2];
	T* B2 = new T[n2];
	T* C2 = new T[n2];
	T* D2 = new T[n2];

	T t = 0; size_t count = 0;
	T* dy = new T [n1*n2]; //y^{k}

	for (size_t i = 0; i < n1*n2; i++) {
		dy[i]= 1;
	}
	while (normInf(n1*n2, dy) > EPS) {
		count++;
		T const1 = 1 / (h1 * h1); T const2 = 1 / (h2 * h2); T const3 = 1 / tau;	
		//std::cout << "                          t=" << t << "\n";

		for (size_t j = 1; j < n2-1; j++) {
			switch (bc3) {//x1=0
			case BC_1:
				A1[0] = 0; B1[0] = 1; C1[0] = 0; D1[0] = G3(j * h2);
				break;
			case BC_2:
				A1[0] = 0;
				B1[0] = -2 * (1. / (h1*h1) + 1. / tau);
				C1[0] = +2. / (h1*h1);
				D1[0] = -(2. / tau * y1[0][j] + phi[0][j] + (y1[0][j + 1] - 2 * y1[0][j] + y1[0][j - 1]) / (h2 * h2)
					- 2 * G3( j * h2) / h1);
				//std::cout << A1[0] << " " << B1[0] << " " << C1[0] << " " << D1[0] << "t="<<t<<"\n";
				//std::cout << "fuuuuuuuuuuuuuuuuuuck";
				break;
			}
			switch (bc4) {//x1=L1
			case BC_1:
				A1[n1 - 1] = 0; B1[n1 - 1] = 1; C1[n1 - 1] = 0; D1[n1 - 1] = G4(j * h2);
				break;
			case BC_2:
				A1[n1 - 1] = 2. / (h1 *h1);
				B1[n1 - 1] = -2 * (1. / (h1*h1)  + 1. / tau);
				C1[n1 - 1] = 0;
				D1[n1 - 1] = -(2. / tau * y1[n1-1][j] + phi[n1-1][j] + (y1[n1-1][j + 1] - 2 * y1[n1-1][j] + y1[n1-1][j - 1]) / (h2 * h2)
					+	2. * G4(j * h2) / h2);
				//std::cout << A1[n1 - 1] << " " << B1[n1 - 1] << " " << C1[n1 - 1] << " " << D1[n1 - 1] << "\n";
				//std::cout << "fuuuuuuuuuuuuuuuuuuck";
				break;
			}
			for (size_t k = 1; k < n1-1; k++) {
				A1[k] = const1;
				B1[k] = -2*(const1+const3);
				C1[k] = const1;
				D1[k] = -(2 * y1[k][j] / tau +((y1[k][j - 1] - 2 * y1[k][j] + y1[k][j + 1]) * const2) + phi[k][j]);
				//std::cout << A1[k] << " " << B1[k] << " " << C1[k] << " " << D1[k] << "\n";
			}
			T* vec = new T[n1];
			solve_3d(n1, A1, B1, C1, D1, vec);
			//std::cout << "\ncтрочка j\n";
			for (size_t it = 0; it < n1; it++) {
					ymid[it][j] = vec[it];
					//std::cout << ymid[it][j] << "  ";
			}
			delete[] vec;
		}
		//std::cout << "\ncтрочка iiiiiiiiii\n";
		for (size_t i = 1; i < n1-1; i++) {
			switch (bc1) {//x2=0
			case BC_1:
				A2[0] = 0; B2[0] = 1; C2[0] = 0; D2[0] = G1(i * h1);
				//std::cout << A2[0] << " " << B2[0] << " " << C2[0] << " " << D2[0] << "\n";
				break;
			case BC_2:
				A2[0] = 0;
				B2[0] = 2 * (1. / (h2) + 1. / tau);
				C2[0] = 2. / (h2 );
				D2[0] = 2. / tau * ymid[i][1]/* + phi[i][0]*/ + (ymid[i + 1][1] - ymid[i][1]) / (h1)
					- 2 * G1(i * h1) / h1;
				//A2[0] = 0;
				//B2[0] = 2 * (1. / h2 + 1. / tau);
				//C2[0] = -2. / (h2);
				//D2[0] = 2. / tau * ymid[i][0] /*+ f(0, j * h2)*/ + (ymid[i+1][0] - 2 * ymid[i][0] + ymid[i-1][0]) / (h1 * h1)
				//	- 2 * G3(i * h1) / h2;
				//std::cout << A2[0] << " " << B2[0] << " " << C2[0] << " " << D2[0] << "\n";
				std::cout << "fuuuuuuuuuuuuuuuuuuck";
				break;
			}
			switch (bc2) {//x2=L2
			case BC_1:
				A2[n2 - 1] = 0; B2[n2 - 1] = 1; C2[n2 - 1] = 0; D2[n2 - 1] = G2(i * h1);
				//std::cout << A2[n2 - 1] << " " << B2[n2 - 1] << " " << C2[n2 - 1] << " " << D2[n2 - 1] << "\n";
				break;
			case BC_2:
				A2[n2 - 1] = -2. / (h2);
				B2[n2 - 1] = 2 * (1. / (h2) + 1. / tau);
				C2[n2 - 1] = 0;
				D2[n2 - 1] = -(2. / tau * ymid[i][n2 - 2] /*+ phi[i][0] */+(ymid[i + 1][n2 - 2] - ymid[i][n2 - 2]) / (h1)
					+2 * G2(i * h1) / h1);
				//std::cout << "fuuuuuuuuuuuuuuuuuuck";
				break;
			}
			
			
			for (size_t k = 1; k < n2-1; k++) {
				A2[k] = const2;
				B2[k] = -2 * (const2 + const3);
				C2[k] = const2;
				D2[k] = -(2 * ymid[i][k] / tau + ((ymid[i-1][k] - 2 * ymid[i][k] + ymid[i+1][k]) * const1) * ymid[i][k] + phi[i][k]);
				//std::cout << A2[k]<<" "<<B2[k]<<" "<<C2[k]<<" "<<D2[k]<<"\n";
			}
			solve_3d(n2, A2, B2, C2, D2, y2[i]);

		}

		for (size_t j = 0; j < n2; j++) {
			y2[0][j] = G1(j * h2);
			y2[n1 - 1][j] = G2(j * h2);
		}
		size_t ind = 0;
		for (size_t i = 0; i < n1; i++) {
			for (size_t j = 0; j < n2; j++) {
				dy[ind]= y2[i][j]-y1[i][j];
				y1[i][j] = y2[i][j];
				ind++;
				//std::cout << y2[i][j] << "  ";
			}
		}





		if (abs(t - tmax / 2) < tau / 2 ) {
			std::cout << "\n50%";

		}
		else {
			if (abs(t - tmax / 4) < tau / 2) {
				std::cout << "\n25%";
			}
			else {
				if (abs(t - 3 * tmax / 4) < tau / 2) {
					std::cout << "\n75%";
				}
			}
		}


		t += tau;
	}
	std::cout << "\n100%\n";
	std::cout << count << "\n";

	for (size_t i = 0; i < n1; i++) {
		for (size_t j = 0; j < n2; j++) {
			out << std::setprecision(13) << y2[i][j]<<" ";
		}
		out << std::endl;
	}
	out.close();

	out.open("D:\\СЕТКА.txt");

	out << std::setprecision(13) << h1 << " " << h2 << " " << n1 << " " << n2;
	out.close();



	delete[] A1; delete[] B1; delete[] C1; delete[] D1; delete[] A2; delete[] B2; delete[] C2; delete[] D2;
	for (size_t i = 0; i < n1; i++) {
		delete phi[i];
	}
	delete[] phi;

	for (size_t i = 0; i < n1; i++) {
		delete y1[i];
		delete ymid[i];
		delete y2[i];
	}
	delete[] y1; delete[] ymid; delete[] y2;
	delete[] dy;
}





void showTable1(T L1, T L2, T(*f)(T, T),
	BOUNDARY_CONDITIONS bc1, T(*G1)(T),
	BOUNDARY_CONDITIONS bc2, T(*G2)(T),
	BOUNDARY_CONDITIONS bc3, T(*G3)(T),
	BOUNDARY_CONDITIONS bc4, T(*G4)(T)) {

	std::cout << "         n &          Шаг сетки tau_n       &        Норма невязки      &       отношение ошибок     &      Порядок сходимости p_n    \n\n";
	T normLast;
	T h1 = 0.5; T h2 = 0.5;
	T tau0 = h1 * h2;
	for (size_t kof = 1; kof <=5; kof++) {
		std::cout << std::setw(12) << kof << " & " << '\t';
		T tau = tau0 / pow(8, kof - 1); 
		h1 = h1 / 2.; h2 = h2/ 2.;					//шаги сетки

		size_t n1 = size_t(round(L1 / h1)) + 1;	//количество узлов
		size_t n2 = size_t(round(L2 / h2)) + 1;



		T** phi = new T * [n1];
		for (size_t i = 0; i < n1; i++) {
			phi[i] = new T[n2];
			for (size_t j = 0; j < n2; j++) {
				phi[i][j] = f(i * h1, j * h2);
			}
		}


		T** y1 = new T * [n1]; //y^{k}
		for (size_t i = 0; i < n1; i++) {
			y1[i] = new T[n2];
		}

		T** ymid = new T * [n1]; //y^{k+1/2}
		for (size_t i = 0; i < n1; i++) {
			ymid[i] = new T[n2];
		}

		T** y2 = new T * [n1]; //y^{k+1}
		for (size_t i = 0; i < n1; i++) {
			y2[i] = new T[n2];
		}


		//начальные данные U0
		for (size_t i = 0; i < n1; i++) {
			for (size_t j = 0; j < n2; j++) {
				y1[i][j] = 1;
			}
		}


		T* A1 = new T[n1];
		T* B1 = new T[n1];
		T* C1 = new T[n1];
		T* D1 = new T[n1];

		T* A2 = new T[n2];
		T* B2 = new T[n2];
		T* C2 = new T[n2];
		T* D2 = new T[n2];

		T t = 0; size_t count = 0;
		T* dy = new T[n1 * n2]; //y^{k}

		for (size_t i = 0; i < n1 * n2; i++) {
			dy[i] = 1;
		}
		while (normInf(n1 * n2, dy) > EPS) {
			count++;
			T const1 = 1 / (h1 * h1); T const2 = 1 / (h2 * h2); T const3 = 1 / tau;
			//std::cout << "                          t=" << t << "\n";

			for (size_t j = 1; j < n2 - 1; j++) {
				switch (bc3) {//x1=0
				case BC_1:
					A1[0] = 0; B1[0] = 1; C1[0] = 0; D1[0] = G3(j * h2);
					break;
				case BC_2:
					A1[0] = 0;
					B1[0] = -2 * (1. / (h1 * h1) + 1. / tau);
					C1[0] = +2. / (h1 * h1);
					D1[0] = -(2. / tau * y1[0][j] + phi[0][j] + (y1[0][j + 1] - 2 * y1[0][j] + y1[0][j - 1]) / (h2 * h2)
						- 2 * G3(j * h2) / h1);
					//std::cout << A1[0] << " " << B1[0] << " " << C1[0] << " " << D1[0] << "t="<<t<<"\n";
					//std::cout << "fuuuuuuuuuuuuuuuuuuck";
					break;
				}
				switch (bc4) {//x1=L1
				case BC_1:
					A1[n1 - 1] = 0; B1[n1 - 1] = 1; C1[n1 - 1] = 0; D1[n1 - 1] = G4(j * h2);
					break;
				case BC_2:
					A1[n1 - 1] = 2. / (h1 * h1);
					B1[n1 - 1] = -2 * (1. / (h1 * h1) + 1. / tau);
					C1[n1 - 1] = 0;
					D1[n1 - 1] = -(2. / tau * y1[n1 - 1][j] + phi[n1 - 1][j] + (y1[n1 - 1][j + 1] - 2 * y1[n1 - 1][j] + y1[n1 - 1][j - 1]) / (h2 * h2)
						+ 2. * G4(j * h2) / h2);
					//std::cout << A1[n1 - 1] << " " << B1[n1 - 1] << " " << C1[n1 - 1] << " " << D1[n1 - 1] << "\n";
					//std::cout << "fuuuuuuuuuuuuuuuuuuck";
					break;
				}
				for (size_t k = 1; k < n1 - 1; k++) {
					A1[k] = const1;
					B1[k] = -2 * (const1 + const3);
					C1[k] = const1;
					D1[k] = -(2 * y1[k][j] / tau + ((y1[k][j - 1] - 2 * y1[k][j] + y1[k][j + 1]) * const2) + phi[k][j]);
					//std::cout << A1[k] << " " << B1[k] << " " << C1[k] << " " << D1[k] << "\n";
				}
				T* vec = new T[n1];
				solve_3d(n1, A1, B1, C1, D1, vec);
				//std::cout << "\ncтрочка j\n";
				for (size_t it = 0; it < n1; it++) {
					ymid[it][j] = vec[it];
					//std::cout << ymid[it][j] << "  ";
				}
				delete[] vec;
			}
			//std::cout << "\ncтрочка iiiiiiiiii\n";
			for (size_t i = 1; i < n1 - 1; i++) {
				switch (bc1) {//x2=0
				case BC_1:
					A2[0] = 0; B2[0] = 1; C2[0] = 0; D2[0] = G1(i * h1);
					//std::cout << A2[0] << " " << B2[0] << " " << C2[0] << " " << D2[0] << "\n";
					break;
				case BC_2:
					A2[0] = 0;
					B2[0] = 2 * (1. / (h2)+1. / tau);
					C2[0] = 2. / (h2);
					D2[0] = 2. / tau * ymid[i][1]/* + phi[i][0]*/ + (ymid[i + 1][1] - ymid[i][1]) / (h1)
						-2 * G1(i * h1) / h1;
					//A2[0] = 0;
					//B2[0] = 2 * (1. / h2 + 1. / tau);
					//C2[0] = -2. / (h2);
					//D2[0] = 2. / tau * ymid[i][0] /*+ f(0, j * h2)*/ + (ymid[i+1][0] - 2 * ymid[i][0] + ymid[i-1][0]) / (h1 * h1)
					//	- 2 * G3(i * h1) / h2;
					//std::cout << A2[0] << " " << B2[0] << " " << C2[0] << " " << D2[0] << "\n";
					std::cout << "fuuuuuuuuuuuuuuuuuuck";
					break;
				}
				switch (bc2) {//x2=L2
				case BC_1:
					A2[n2 - 1] = 0; B2[n2 - 1] = 1; C2[n2 - 1] = 0; D2[n2 - 1] = G2(i * h1);
					//std::cout << A2[n2 - 1] << " " << B2[n2 - 1] << " " << C2[n2 - 1] << " " << D2[n2 - 1] << "\n";
					break;
				case BC_2:
					A2[n2 - 1] = -2. / (h2);
					B2[n2 - 1] = 2 * (1. / (h2)+1. / tau);
					C2[n2 - 1] = 0;
					D2[n2 - 1] = -(2. / tau * ymid[i][n2 - 2] /*+ phi[i][0] */ + (ymid[i + 1][n2 - 2] - ymid[i][n2 - 2]) / (h1)
						+2 * G2(i * h1) / h1);
					//std::cout << "fuuuuuuuuuuuuuuuuuuck";
					break;
				}


				for (size_t k = 1; k < n2 - 1; k++) {
					A2[k] = const2;
					B2[k] = -2 * (const2 + const3);
					C2[k] = const2;
					D2[k] = -(2 * ymid[i][k] / tau + ((ymid[i - 1][k] - 2 * ymid[i][k] + ymid[i + 1][k]) * const1) * ymid[i][k] + phi[i][k]);
					//std::cout << A2[k]<<" "<<B2[k]<<" "<<C2[k]<<" "<<D2[k]<<"\n";
				}
				solve_3d(n2, A2, B2, C2, D2, y2[i]);

			}

			for (size_t j = 0; j < n2; j++) {
				y2[0][j] = G1(j * h2);
				y2[n1 - 1][j] = G2(j * h2);
			}
			size_t ind = 0;
			for (size_t i = 0; i < n1; i++) {
				for (size_t j = 0; j < n2; j++) {
					dy[ind] = y2[i][j] - y1[i][j];
					y1[i][j] = y2[i][j];
					ind++;
					//std::cout << y2[i][j] << "  ";
				}
			}







			t += tau;
		}
		//std::cout << "\n100%\n";
		//std::cout << count << "\n";


		T* norma = new T[n1 * n2];
		T mx = 0; size_t ind = 0;
		for (size_t i = 0; i < n1; i++) {
			T my = 0;
			for (size_t j = 0; j < n2; j++) {
				norma[ind] = y2[i][j] - (mx * mx + my * my);
				ind++;
				my += h2;
			}
			mx += h1;
		}
		T normNew = normInf(n1 * n2, norma);

		//std::cout << "norma=" << y2[3][3];
		delete[] norma;
		std::cout << std::setw(20) << tau << " & " << "\t\t";
		if (kof == 1) {
			normLast = normNew;
			std::cout << std::setw(20) << normNew << " &    " << "         ---            " << "   &      " << "             ---      " << "\\\\\ \\hline\n";

		}
		else {
			T z = normNew / normLast;
			T pn = log10(z) / log10(0.5);
			std::cout << std::setw(20) << normNew << " & " << "\t\t";
			std::cout << std::setw(20) << z << " & " << "\t\t";
			std::cout << std::setw(17) << pn; std::cout << "\\\\\ \\hline\n";
			normLast = normNew;
		}

		delete[] A1; delete[] B1; delete[] C1; delete[] D1; delete[] A2; delete[] B2; delete[] C2; delete[] D2;
		for (size_t i = 0; i < n1; i++) {
			delete phi[i];
		}
		delete[] phi;

		for (size_t i = 0; i < n1; i++) {
			delete y1[i];
			delete ymid[i];
			delete y2[i];
		}
		delete[] y1; delete[] ymid; delete[] y2;
		delete[] dy;


	
	}


}
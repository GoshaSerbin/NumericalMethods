#include"header.h"

T Arctg(T x, T y) {
	T phi = 0;
	if (abs(x) < EPSZERO) {
		if (y > 0) phi = PI / 2;
		if (y < 0) phi = 3 * PI / 2;
	}
	else {
		if (x > 0 && y > 0) {
			phi = atan(y / x);
		}
		else {
			if (x > 0 && y < 0) {
				phi = atan(y / x) + 2 * PI;
			}
			else {
				if (x < 0 && y >= 0) {
					phi = atan(y / x) + PI;
				}
				else {
					if (x < 0 && y < 0) {
						phi = atan(y / x) + PI;
					}
				}
			}
		}

	}
	return phi;
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



void SolveQuadratures(T(*K)(T, T), T(*f)(T), T a, T b, T lambda) {
	T h = 1e-2;
	size_t N = size_t(round((b-a) / h)) + 1;		//количество узлов


	T** A = new T * [N];
	for (size_t i = 0; i < N; i++)	A[i] = new T[N];

	T* rpart = new T[N];
	T* y = new T[N]; //решение

	std::ofstream out;	
	out.open("D:\\МЕТОД КВАДРАТУР.txt");

	T x = a; 
	for (size_t i = 0; i < N; i++) {
		T s = a;
		for(size_t k = 0; k < N; k++) {
			T a;
			if (k == 0 || k == N - 1)	a = h / 2;
			else	a = h;

			A[i][k] = -lambda * a * K(x, s);
			if (i == k) A[i][k] += 1;
			s += h;
		}
		rpart[i] = f(x);
		x += h;
	}
	//std::cout << "hey";
	gauss_method(N, A, rpart, y);


	x = a;
	for (size_t i = 0; i < N; i++) {
		out << " " << std::setprecision(13) << x << " " << y[i] << std::endl;
		x += h;
	}

	out.close();

	//ошибка
	T* dy = new T[N];
	x = 0;
	for (size_t i = 0; i < N; i++) {
		dy[i] = y[i] - 1;
		x += h;
	}
	std::cout<< "Ошибка = "<<normInf(N, dy) << "\n";;
	delete[] dy;

	for (size_t i = 0; i < N; i++) {
		delete[] A[i];
	}
	delete[] A;
	delete[] rpart;
	delete[] y;
	
}

void SolveSimpleIteration(T(*K)(T, T), T(*f)(T), T a, T b, T lambda) {
	T h = 1e-2;
	size_t N = size_t(round((b - a) / h)) + 1;		//количество узлов
	T* y0 = new T[N]; 
	T* y1 = new T[N];
	size_t counter = 0;

	std::ofstream out;
	out.open("D:\\МЕТОД ПРОСТОЙ ИТЕРАЦИИ.txt");

	T* delta = new T[N];
	for (size_t i = 0; i < N; i++) {
		delta[i] = 1;
		y0[i] = 1;
	}

	while (normInf(N, delta) > EPS) {
		counter++;
		T x = a;//квадратурные формулы гаусса или нет
		//T h_int = h / 2.;
		//size_t M = 2 * N + 1;
		for (size_t i = 0; i < N; i++) {
			T integral = 0;
			T s = 0;
			for (size_t k = 0; k < N-1; k++) {
				integral += (y0[k] * K(x, s) + y0[k + 1] * K(x, s+h));
				s += h;
			}
			y1[i] = f(x) + lambda * h/2 * integral;
			x += h;
		}


		for (size_t i = 0; i < N; i++) {
			delta[i] = y1[i] - y0[i];
			y0[i] = y1[i];
		}
	}
	T x = a;
	for (size_t i = 0; i < N; i++) {
		out << " " << std::setprecision(13) << x << " " << y1[i] << std::endl;
		x += h;
	}
	std::cout << "counter=" << counter << "\n";

	out.close();

	//ошибка
	T* dy = new T[N];
	x = 0;
	for (size_t i = 0; i < N; i++) {
		dy[i] = y1[i] - (x * x);
		x += h;
	}
	std::cout << "Ошибка = " << normInf(N, dy) << "\n";
	delete[] dy;
	delete[] y0;
	delete[] y1;
}


void SolveFredgolm2(size_t m, T(*psi[10])(T), T(*phi[10])(T), T(*f)(T), T a, T b, T lambda) {
	
	T h = 0.005;
	size_t N = size_t(round((b - a) / h)) + 1;		//количество узлов
	T** alpha = new T * [m]; T* beta = new T [m]; 
	for (size_t i = 0; i < m; i++) {
		alpha[i] = new T[m];
	}

	std::ofstream out;
	out.open("D:\\ВЫРОЖДЕННОЕ ЯДРО.txt");

	T* y = new T[N];

	for (size_t i = 0; i < m; i++) {
		for (size_t j = 0; j < m; j++) {
			T integral1 = 0;
			T x = a;
			for (size_t k = 0; k < N-1; k++) {//2N+1
				integral1 += (phi[j](x)* psi[i](x)+ phi[j](x+h) * psi[i](x+h));
				x += h;
			}
			integral1 *= h/2.;
			alpha[i][j] = integral1;
		}
		T integral2 = 0;
		T x = a;
		for (size_t k = 0; k < N-1; k++) {//2N+1
			integral2 += (f(x) * psi[i](x) + f(x + h) * psi[i](x + h));
			x += h;
		}
		integral2 *= h / 2.;
		beta[i] = integral2;
	}



	T* C = new T[m];
	T** A = new T * [m];
	for (size_t i = 0; i < m; i++)	A[i] = new T[m];

	for (size_t i = 0; i < m; i++) {
		for (size_t j = 0; j < m; j++) {
			A[i][j] = -lambda * alpha[i][j];
			if (i == j) A[i][j] += 1;
		}
	}
	gauss_method(m, A, beta, C);

	T x = a;
	for (size_t i = 0; i < N; i++) {
		T sum = 0;
		for (size_t j = 0; j < m; j++) {
			sum += C[j] * phi[j](x);
		}
		y[i] = f(x) + lambda * sum;
		out << " " << std::setprecision(13) << x << " " <<  y[i] << std::endl;
		x += h;
	}

	out.close();

	T* dy = new T[N];
	x = 0;
	for (size_t i = 0; i < N; i++) {
		dy[i] = y[i] - 1;
		x += h;
	}
	std::cout << "Ошибка = " << normInf(N, dy) << "\n";

	delete[] C;
	delete[] y;
	delete[] beta;
	for (size_t i = 0; i < m; i++) {
		delete[] alpha[i];
		delete[] A[i];
	}
	delete[] A;
	delete[] alpha;
}


void SolveFredgolm2(T(*K)(T, T), T(*f)(T), T a, T b, T lambda, SOLVE_TYPE st) {
	switch (st) {
	case ST_QUADRATURES:
		SolveQuadratures(K, f, a, b, lambda);

		break;
	case ST_SIMPLE_ITERATION:
		SolveSimpleIteration(K, f, a, b, lambda);

	}
}


T Qx(T x, T y) {
	return -y / (2 * PI) / (x * x + y * y);
}

T Qy(T x, T y) {
	return x / (2 * PI) / (x * x + y * y);
}

void SolveSingularIntegralEquation(T(*f)(T, T)) {
	size_t N = 160;		//количество узлов
	T dl = 2 * PI / N;
	T** normalVector = new T*[N];
	for (size_t i = 0; i < N; i++) {
		normalVector[i] = new T[2];
		normalVector[i][0] = cos(dl * (i+1 - 0.5));
		normalVector[i][1] = sin(dl * (i+1 - 0.5));
	}

	std::ofstream out;
	out.open("D:\\СИНГУЛЯРНОЕ УРАВНЕНИЕ.txt");

	T* g = new T[N+1];

	T** A = new T * [N+1];
	for (size_t i = 0; i < N+1; i++)	A[i] = new T[N+1];

	T* rpart = new T[N + 1];
	for (size_t i = 0; i < N; i++) {
		rpart[i] = f(cos(dl * (i+1 - 0.5)), sin(dl * (i+1 - 0.5)));//Вообще k_i=n_i
	}
	rpart[N] = 0;

	for (size_t i = 0; i < N; i++) {
		for (size_t j = 0; j < N; j++) {
			T sum_x = Qx(cos(dl * (i+1 - 0.5))- cos(dl * (j+1 - 1)), sin(dl * (i+1 - 0.5))- sin(dl * (j+1 - 1)));
			T sum_y = Qy(cos(dl * (i+1 - 0.5)) - cos(dl * (j+1 - 1)), sin(dl * (i+1 - 0.5)) - sin(dl * (j+1 - 1)));
			A[i][j] = normalVector[i][0] * sum_x + normalVector[i][1] * sum_y;
			A[i][j] *= dl;
		}
		A[i][N] = 1;
	}
	for (size_t j = 0; j < N; j++) {
		A[N][j] = 1;
	}
	A[N][N] = 0;
	//for (size_t i = 0; i < N + 1; i++) {
	//	for (size_t j = 0; j < N + 1; j++) {
	//		std::cout << A[i][j] << " ";
	//	}
	//	std::cout << rpart[i]<<"\n";
	//}


	gauss_method(N+1, A, rpart, g);

	T phi = 0;
	for (size_t i = 0; i < N; i++) {
		out << " " << std::setprecision(13) << phi << " " << g[i] << std::endl;
		phi += dl;
	}
	std::cout << "\nR=" << g[N]<<"\n";
	out.close();


	//std::ofstream outE;
	//outE.open("D:\\фи.txt");
	//T x = -1; T y = -sqrt(1 - x * x);
	//for (size_t i = 0; i < 2*N; i++) {
	//	y = -sqrt(1 - x * x);
	//	outE << " " << std::setprecision(13) << x << " " << y << " " << Arctg(x,y) << std::endl;
	//	x += dl/2;
	//}
	//outE.close();

	delete[] g;
	delete[] normalVector[0]; delete[] normalVector[1]; delete[] normalVector;
	for (size_t i = 0; i < N+1; i++) {	
		delete[] A[i];

	}
	delete[] A;
	delete[] rpart;
}


//
//void DSolveWaveEquation(T L, T tmax, T a2, T(*f)(T), T(*ddf)(T), T(*g)(T), T(*phi)(T), T(*psi)(T), DIFFERENTIAL_TYPE dt) {
//
//	T h = 0.005;	T tau = h;				//шаги сетки
//	size_t nx = size_t(round(L / h)) + 1;		//количество узлов
//	size_t nt = size_t(round(tmax / tau)) + 1;
//
//	T* y1 = new T[nx];							//y^{j-1}
//	T* y2 = new T[nx];							//y^{j}
//	T* y3 = new T[nx];							//y^{j+1}
//
//	std::ofstream out;							// поток для записи
//	out.open("D:\\ВОЛНОВОЕ УРАВНЕНИЕ 1.txt");
//	std::ofstream outE;
//	outE.open("D:\\ЭНЕРГИЯ.txt");
//
//
//	T x = 0;
//	for (size_t i = 0; i < nx; i++) {
//		y1[i] = f(x);
//		x += h;
//	}
//	out << std::setprecision(13) << 0;
//	for (size_t i = 0; i < nx; i++) {
//		out << " " << std::setprecision(13) << y1[i];
//	}
//	out << std::endl;
//	T dE = 0;
//	outE << std::setprecision(13) << 0 << " " << std::setprecision(13) << dE << std::endl;
//
//
//	x = h;
//	y2[0] = phi(tau); y2[nx - 1] = psi(tau);
//	if (dt == DT_ANALYTIC) {
//		for (size_t i = 1; i < nx - 1; i++) {
//			y2[i] = y1[i] + tau * g(x) + a2 * tau * tau / 2 * ddf(x);
//			x += h;
//		}
//	}
//	else {
//		for (size_t i = 1; i < nx - 1; i++) {
//			T fxx = (f(x + h) - 2 * f(x) + f(x - h)) / h / h;
//			y2[i] = y1[i] + tau * g(x) + a2 * tau * tau / 2 * fxx;
//			x += h;
//		}
//	}
//
//	out << std::setprecision(13) << tau;
//	for (size_t i = 0; i < nx; i++) {
//		out << " " << std::setprecision(13) << y2[i];
//	}
//	out << std::endl;
//
//	bool check = false;
//	T t = 2 * tau;
//	std::cout << "\n0%";
//	while (t <= tmax + EPS) {
//
//		y3[0] = phi(t); y3[nx - 1] = psi(t);
//
//		for (size_t i = 1; i < nx - 1; i++) {
//			y3[i] = pow(tau / h, 2) * a2 * (y2[i + 1] - 2 * y2[i] + y2[i - 1]) - y1[i] + 2 * y2[i];
//		}
//
//
//		T S1 = 0; T S2 = 0;
//		for (size_t i = 1; i < nx - 1; i++) {
//			S1 += (y2[i] - y1[i]) * (y3[i] - 2 * y2[i] + y1[i]);
//			S2 += (y2[i] - y2[i - 1]) * ((y2[i] - y1[i]) - (y2[i - 1] - y1[i - 1]));
//		}
//		S2 += (y2[nx - 1] - y2[nx - 1 - 1]) * ((y2[nx - 1] - y1[nx - 1]) - (y2[nx - 1 - 1] - y1[nx - 1 - 1]));
//		dE = h * S1 / a2 / (pow(tau, 3)) + S2 / h / tau + (-(y2[nx - 1] - y2[nx - 2]) * (y2[nx - 1] - y1[nx - 1]) + (y2[1] - y2[0]) * (y2[0] - y1[0])) / (h * tau);
//		//std::cout << " " << /*S1 << " " << S2 << " " <<*/ dE << "\n";
//		outE << std::setprecision(13) << t << " " << std::setprecision(13) << dE << std::endl;
//
//		for (size_t i = 0; i < nx; i++) {
//			y1[i] = y2[i];
//			y2[i] = y3[i];
//		}
//
//		out << std::setprecision(13) << t;
//		for (size_t i = 0; i < nx; i++) {
//			out << " " << std::setprecision(13) << y3[i];
//		}
//		out << std::endl;
//
//		t += tau;
//
//		if (abs(t - tmax / 2) < tau / 2 && check) {
//			std::cout << "\n50%";
//
//		}
//		else {
//			if (abs(t - tmax / 4) < tau / 2) {
//				std::cout << "\n25%";
//			}
//			else {
//				if (abs(t - 3 * tmax / 4) < tau / 2) {
//					std::cout << "\n75%";
//				}
//			}
//		}
//	}
//	std::cout << "\n100%\n";
//	out.close();
//	outE.close();
//
//	out.open("D:\\СЕТКА.txt");
//	x = 0;
//	for (size_t i = 0; i < nx; i++) {
//		out << std::setprecision(13) << x << std::endl;
//		x += h;
//	}
//
//	delete[] y1; delete[] y2; delete[] y3;
//}

void MakeTable(T L, T tmax, T a2, T(*f)(T), T(*df)(T), T(*g)(T), T(*phi)(T), T(*psi)(T)) {


	T tau0 = 2e-2;
	T h0 = 2e-1;
	size_t nx = size_t(round(L / h0)) + 1;
	size_t nt = size_t(round(tmax / tau0)) + 1;

	std::cout << "         n   &    Шаг сетки tau_n    &  Норма ошибки         &     Отношение w_n     & Порядок сходимости & оценка ошибки err \n\n";
	T normLast;

	for (size_t n = 1; n <= 5; n++) {
		T tau = tau0 / pow(2, n - 1);	T h = h0 / pow(2, n - 1);

		T* ySol1 = new T[nx];
		T* ySol2 = new T[size_t(round(L / (h / 2))) + 1];

		std::cout << std::setw(12) << n << " & " << '\t';

		//решаем для tau_n
		{
			nx = size_t(round(L / h)) + 1;
			nt = size_t(round(tmax / tau)) + 1;

			T* y1 = new T[nx];							//y^{j-1}
			T* y2 = new T[nx];							//y^{j}
			T* y3 = new T[nx];							//y^{j+1}

			T x = 0;
			for (size_t i = 0; i < nx; i++) {
				y1[i] = f(x);
				x += h;
			}

			x = h;
			y2[0] = phi(tau); y2[nx - 1] = psi(tau);
			for (size_t i = 1; i < nx - 1; i++) {
				y2[i] = y1[i] + tau * g(x) + a2 * tau * tau / 2 * df(x);
				x += h;
			}

			T t = 2 * tau;
			while (t <= tmax + EPS) {
				y3[0] = phi(t); y3[nx - 1] = psi(t);

				for (size_t i = 1; i < nx - 1; i++) {
					y3[i] = pow(tau / h, 2) * a2 * (y2[i + 1] - 2 * y2[i] + y2[i - 1]) - y1[i] + 2 * y2[i];
				}

				for (size_t i = 0; i < nx; i++) {
					y1[i] = y2[i];
					y2[i] = y3[i];
				}

				t += tau;

			}

			for (size_t i = 0; i < nx; i++) {
				ySol1[i] = y2[i];
			}
			delete[] y1; delete[] y2; delete[] y3;
		}
		T nx1 = nx;

		//решаем для tau_{n+1}
		{
			tau = tau / 2; h = h / 2;
			nx = size_t(round(L / h)) + 1;
			nt = size_t(round(tmax / tau)) + 1;

			T* y1 = new T[nx];							//y^{j-1}
			T* y2 = new T[nx];							//y^{j}
			T* y3 = new T[nx];							//y^{j+1}

			T x = 0;
			for (size_t i = 0; i < nx; i++) {
				y1[i] = f(x);
				x += h;
			}

			x = h;
			y2[0] = phi(tau); y2[nx - 1] = psi(tau);
			for (size_t i = 1; i < nx - 1; i++) {
				y2[i] = y1[i] + tau * g(x) + a2 * tau * tau / 2 * df(x);
				x += h;
			}

			T t = 2 * tau;
			while (t <= tmax + EPS) {
				y3[0] = phi(t); y3[nx - 1] = psi(t);

				for (size_t i = 1; i < nx - 1; i++) {
					y3[i] = pow(tau / h, 2) * a2 * (y2[i + 1] - 2 * y2[i] + y2[i - 1]) - y1[i] + 2 * y2[i];
				}

				for (size_t i = 0; i < nx; i++) {
					y1[i] = y2[i];
					y2[i] = y3[i];
				}

				t += tau;

			}

			for (size_t i = 0; i < nx; i++) {
				ySol2[i] = y2[i];
			}
			delete[] y1; delete[] y2; delete[] y3;
		}
		T nx2 = nx;

		T* dy = new T[nx1];

		for (size_t i = 0; i < nx1; i++) {
			dy[i] = ySol2[2 * i] - ySol1[i];
		}

		T normNew = normInf(nx1, dy);
		delete[] dy;

		T w;
		T pn;
		std::cout << std::setw(20) << tau << " & " << "\t";

		if (n == 1) {
			std::cout << " ------                            ----- \\\\\ \\hline\n";
			normLast = normNew;
		}
		else {
			w = normLast / normNew;
			pn = log10(w) / log10(2);
			std::cout << std::setw(20) << normLast << " & " << "\t";
			std::cout << std::setw(20) << w << " & " << "\t";
			std::cout << std::setw(17) << pn << " & " << "\t";
			std::cout << std::setw(13) << normLast / (1 - pow(2, -pn)); std::cout << "\\\\\ \\hline\n";
			normLast = normNew;
		}



	}

	std::cout << "----------------------------------------------------------------------------------------\n\n";



}

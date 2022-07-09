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
	//std::cout << "коснтанты:\n";
	//for (size_t i = 0; i < n; i++) {
	//	std::cout << a[i] << "\t";
	//}
	//std::cout << "\n";
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

T integral(T* y, size_t nx, T a, T b) {
	T I = 0;
	for (size_t i = 0; i < nx-1; i++) {
		I += y[i]+y[i+1];
	}
	I *= (b - a) / (1.*nx-1)/2;
	return I;
}


void DSolveHeatEquation(T L, T tmax, T c, T ro, T(*K)(T), THERMAL_CONDUCTIVITY_K tc, T(*u0)(T), T (*ubound)(T), T(*P)(T, T, T), T t0, T Q) {

	T tau = 2e-4; T h = 2e-1;					//шаги сетки
	size_t nx = size_t(round(L / h)) + 1;	//количество узлов
	size_t nt = size_t(round(tmax / tau)) + 1;

	T* y = new T[nx];						//y^j
	T* y1 = new T[nx];					//y^{j+1}

	T sigma = 0.5;

	T Tmax = 0;

	T maxU = 0;

	if (tc == TC_X) {

		std::ofstream out;				// поток для записи
		out.open("D:\\УРАВНЕНИЕ ТЕПЛОПРОВОДНОСТИ 1.txt");

		T x = 0;
		for (size_t i = 0; i < nx; i++) {
			y[i] = u0(x);
			x += h;
			//if (x > 1) std::cout << x << " ";
		}



		T* a = new T[nx];
		a[0] = 0;//добавляем чтобы совпадало с A[1] <-> a[1]
		x = 0;
		for (size_t i = 1; i < nx; i++) {//методом трапеций самое то
			size_t k = 20;
			T dh = h / k;
			T I = 0;
			for (size_t j = 0; j < k; j++) {
				I += 1 / K(x + dh / 2);
				x += dh;
			}
			I *= dh;
			//считаем интеграл => ai
			//не зависят от t?
			a[i] = h / I;
			//std::cout << a[i]<<"\n";
			x += h;
		}

		T* A = new T[nx];
		T* B = new T[nx];
		T* C = new T[nx];
		T* F = new T[nx];
		
		T t = 0;
		while (t <= tmax + EPS) {

			out << std::setprecision(13) << t;
			for (size_t i = 0; i < nx; i++) {
				out << " " << std::setprecision(13) << y[i];
			}
			out << std::endl;

			T temperature = y[nx / 2];
			if (temperature > maxU) {
				Tmax = t;
				maxU = temperature;
			}

			T* wminus = new T[nx];
			T* wplus = new T[nx];
			wminus[0] = 0; wplus[0] = 0;
			for (size_t i = 1; i < nx - 1; i++) {
				wminus[i] = a[i] * (y[i] - y[i - 1]) / h;
				wplus[i] = a[i+1] * (y[i + 1] - y[i]) / h;//a[i+1]????
			//	std::cout << wminus[i] << " " << wplus[i] << "\n";
			}
			//std::cout << "end";
			wminus[nx - 1] = a[nx - 1] * (y[nx - 1] - y[nx - 1 - 1]) / h;


			T kappa = (sigma * a[nx - 1] / h) / (c * ro * h / (2 * tau) + sigma * a[nx - 1] / h);
			T mu = (c * ro * y[nx - 1] * h / (2 * tau) + sigma * P(t + tau, t0, Q) + (1 - sigma) * (P(t, t0, Q) - wminus[nx - 1])) / (c * ro * h / (2 * tau) + sigma * a[nx - 1] / h);



			A[0] = 0; C[0] = 1; B[0] = 0; F[0] = ubound(t+tau);
			for (size_t i = 1; i < nx - 1; i++) {
				A[i] = sigma / h * a[i];
				B[i] = sigma / h * a[i + 1];
				C[i] = -(sigma / h * (a[i] + a[i + 1]) + c * ro * h / tau);
				F[i] = -(c * ro * h / tau * y[i] + (1 - sigma) * (wplus[i] - wminus[i]));
				//std::cout << A[i]<<" "<<B[i]<<" "<<C[i]<<" "<<F[i]<<"\n";
			}//опечатка в методичке..
			//std::cout << "kdfjdkfjdkfjkdfjd\n";
			A[nx - 1] = -kappa; C[nx - 1] = 1; B[nx - 1] = 0; F[nx - 1] = mu;

			//std::cout << kappa << " "<<mu<<"\n";
			//std::cout << (P(t, t0, Q) - wminus[nx - 1]) <<" ";

			//здесь надо чтобы рассматривал какие поставлены условия

			solve_3d(nx, A, C, B, F, y1);
			//y1[nx - 1] = kappa * y1[nx - 2] + mu;


			for (size_t i = 0; i < nx; i++) {
				y[i] = y1[i];
			}

			t += tau;
		}

		out.close();

		out.open("D:\\СЕТКА.txt");
		x = 0;
		for (size_t i = 0; i < nx; i++) {
			out << std::setprecision(13) << x << std::endl;
			x += h;
		}

		delete[] A; delete[] B; delete[] C; delete[] F;
		delete[] a;


	}
	else {
		sigma = 0.5;
		std::ofstream out;				// поток для записи
		out.open("D:\\УРАВНЕНИЕ ТЕПЛОПРОВОДНОСТИ КВАЗИЛИНЕЙНОЕ 1.txt");

		T x = 0;
		for (size_t i = 0; i < nx; i++) {
			y[i] = u0(x);
			x += h;
			//if (x > 1) std::cout << x << " ";
		}

		T* A = new T[nx];
		T* B = new T[nx];
		T* C = new T[nx];
		T* F = new T[nx];

		T t = 0;
		while (t <= tmax + EPS) {

			out << std::setprecision(13) << t;
			for (size_t i = 0; i < nx; i++) {
				out << " " << std::setprecision(13) << y[i];
			}
			out << std::endl;

			T temperature = y[nx / 2];
			if (temperature > maxU) {
				Tmax = t;
				maxU = temperature;
			}


			T* a = new T[nx];
			a[0] = 0;//добавляем чтобы совпадало с A[1] <-> a[1]
			for (size_t i = 1; i < nx; i++) {//методом трапеций самое то
				a[i] = (K(y[i]) + K(y[i - 1]))/2.;
			//std::cout << a[i]<<"\n";
			}


			T* wminus = new T[nx];
			T* wplus = new T[nx];
			wminus[0] = 0; wplus[0] = 0;
			for (size_t i = 1; i < nx - 1; i++) {
				wminus[i] = a[i] * (y[i] - y[i - 1]) / h;
				wplus[i] = a[i+1] * (y[i + 1] - y[i]) / h;//iiiiii????
			//	std::cout << wminus[i] << " " << wplus[i] << "\n";
			}
			//std::cout << "end";
			wminus[nx - 1] = a[nx - 1] * (y[nx - 1] - y[nx - 1 - 1]) / h;

			T* ynew = new T[nx];
			T* y0 = new T[nx];
			for (size_t j = 0; j < nx; j++) {
				ynew[j] = y[j];
				y0[j] = y[j];
			}
			T* dy = new T[nx]; dy[0] = 1;
			size_t count = 1;
			//for (size_t k = 0; k < 2; k++) {
			while(norm1(nx,dy)>0.00001) {
				A[0] = 0; B[0] = 1; C[0] = 0; F[0] = ubound(t+tau);
				count++;
				for (size_t i = 1; i < nx - 1; i++) {
					A[i] = -(K(y0[i]) + K(y0[i - 1])) / 2 / h / h;
					B[i] = c * ro / tau + (K(y0[i + 1]) + K(y0[i]) + K(y0[i]) + K(y0[i - 1])) / 2 / h / h;
					C[i] = -(K(y0[i+1]) + K(y0[i])) / 2 / h / h;
					F[i] = c * ro * y[i] / tau;
				//	std::cout << A[i]<<" "<<B[i]<<" "<<C[i]<<" "<<F[i]<<"\n";
				}//опечатка в методичке..
				//std::cout << "kdfjdkfjdkfjkdfjd\n";
				T kappa = (sigma * a[nx - 1] / h) / (c * ro * h / (2 * tau) + sigma * a[nx - 1] / h);
				T mu = (c * ro * y[nx - 1] * h / (2 * tau) + sigma * P(t + tau, t0, Q) + (1 - sigma) * (P(t, t0, Q) - wminus[nx - 1])) / (c * ro * h / (2 * tau) + sigma * a[nx - 1] / h);
				A[nx - 1] = -kappa; B[nx - 1] = 1; C[nx - 1] = 0; F[nx - 1] = mu;
				//std::cout << kappa << " "<<mu<<"\n";
				//std::cout << (P(t, t0, Q) - wminus[nx - 1]) <<" ";

				//здесь надо чтобы рассматривал какие поставлены условия

				solve_3d(nx, A, B, C, F, ynew);
				for (size_t r = 0; r < nx; r++) {
					dy[r] = y0[r] - ynew[r];
				}
				//y1[nx - 1] = kappa * y1[nx - 2] + mu;
				for (size_t j = 0; j < nx; j++) {
					y0[j] = ynew[j];
				}
			}
			//std::cout << count << " \n";
			for (size_t j = 0; j < nx; j++) {
				y1[j] = ynew[j];
			}
			delete[] y0;
			delete[] ynew;


			for (size_t i = 0; i < nx; i++) {
				y[i] = y1[i];
			}

			t += tau;
			delete[] a;
		}

	out.close();

	out.open("D:\\СЕТКА.txt");
	x = 0;
	for (size_t i = 0; i < nx; i++) {
		out << std::setprecision(13) << x << std::endl;
		x += h;
	}

	delete[] A; delete[] B; delete[] C; delete[] F;

	}

	std::cout << "T^*=" << Tmax << "\n";
}

void DSolveHeatEquation(T L, T tmax, T c, T ro, T(*K)(T), THERMAL_CONDUCTIVITY_K tc, T(*u0)(T), T(*P)(T, T, T), T t0, T Q, T(*ubound)(T)) {

	T tau = 1e-2; T h = 1e-2;				//шаги сетки
	size_t nx = size_t(round(L / h)) + 1;	//количество узлов
	size_t nt = size_t(round(tmax / tau)) + 1;

	T* y = new T[nx];						//y^j
	T* y1 = new T[nx];					//y^{j+1}

	T sigma = 0.5;

	if (tc == TC_X) {

		std::ofstream out;				// поток для записи
		out.open("D:\\УРАВНЕНИЕ ТЕПЛОПРОВОДНОСТИ 2.txt");

		T x = 0;
		for (size_t i = 0; i < nx; i++) {
			y[i] = u0(x);
			x += h;
			//if (x > 1) std::cout << x << " ";
		}



		T* a = new T[nx];
		a[0] = 0;//добавляем чтобы совпадало с A[1] <-> a[1]
		x = 0;
		for (size_t i = 1; i < nx; i++) {//методом трапеций самое то
			size_t k = 20;
			T dh = h / k;
			T I = 0;
			for (size_t j = 0; j < k; j++) {
				I += 1 / K(x + dh / 2);
				x += dh;
			}
			I *= dh;
			//считаем интеграл => ai
			//не зависят от t?
			a[i] = h / I;
			//std::cout << a[i]<<"\n";
			x += h;
		}



		T* A = new T[nx];
		T* B = new T[nx];
		T* C = new T[nx];
		T* F = new T[nx];

		T t = 0;
		while (t <= tmax + EPS) {

			out << std::setprecision(13) << t;
			for (size_t i = 0; i < nx; i++) {
				out << " " << std::setprecision(13) << y[i];
			}
			out << std::endl;


			T* wminus = new T[nx];
			T* wplus = new T[nx];
			wminus[0] = 0; wplus[0] = 0;
			for (size_t i = 1; i < nx - 1; i++) {
				wminus[i] = a[i] * (y[i] - y[i - 1]) / h;
				wplus[i] = a[i + 1] * (y[i + 1] - y[i]) / h;//a[i+1]????
			//	std::cout << wminus[i] << " " << wplus[i] << "\n";
			}
			//std::cout << "end";
			wminus[nx - 1] = a[nx - 1] * (y[nx - 1] - y[nx - 1 - 1]) / h;


			T kappa = (sigma * a[1] / h) / (c * ro * h / (2 * tau) + sigma * a[1] / h);
			T mu = (c * ro * y[0] * h / (2 * tau) + sigma * P(t + tau, t0, Q) + (1 - sigma) * (P(t, t0, Q) + wplus[1])) / (c * ro * h / (2 * tau) + sigma * a[1] / h);


			//тут меняется только эта строчка но мне лень следовать общепринятым духовным ценностям в программировании
			A[nx - 1] = 0; C[nx - 1] = 1; B[nx - 1] = 0; F[nx - 1] = ubound(t + tau);
			for (size_t i = 1; i < nx - 1; i++) {
				A[i] = sigma / h * a[i];
				B[i] = sigma / h * a[i + 1];
				C[i] = -(sigma / h * (a[i] + a[i + 1]) + c * ro * h / tau);
				F[i] = -(c * ro * h / tau * y[i] + (1 - sigma) * (wplus[i] - wminus[i]));
				//std::cout << A[i]<<" "<<B[i]<<" "<<C[i]<<" "<<F[i]<<"\n";
			}//опечатка в методичке..
			//std::cout << "kdfjdkfjdkfjkdfjd\n";
			//A[0] = 0; C[0] = kappa; B[0] =1 ; F[0] = mu;
			A[0] = 0; C[0] = 1; B[0] = -kappa; F[0] = mu;
			//std::cout << kappa << " "<<mu<<"\n";
			//std::cout << (P(t, t0, Q) - wminus[nx - 1]) <<" ";

			//здесь надо чтобы рассматривал какие поставлены условия

			solve_3d(nx, A, C, B, F, y1);
			//y1[nx - 1] = kappa * y1[nx - 2] + mu;

			//std::cout << y[nx - 1] << "\n";
			for (size_t i = 0; i < nx; i++) {
				y[i] = y1[i];
			}

			t += tau;
		}

		out.close();

		out.open("D:\\СЕТКА.txt");
		x = 0;
		for (size_t i = 0; i < nx; i++) {
			out << std::setprecision(13) << x << std::endl;
			x += h;
		}

		delete[] A; delete[] B; delete[] C; delete[] F;
		delete[] a;


	}
	else {
		sigma = 1;
		std::ofstream out;				// поток для записи
		out.open("D:\\УРАВНЕНИЕ ТЕПЛОПРОВОДНОСТИ КВАЗИЛИНЕЙНОЕ 2.txt");

		T x = 0;
		for (size_t i = 0; i < nx; i++) {
			y[i] = u0(x);
			x += h;
			//if (x > 1) std::cout << x << " ";
		}

		T* A = new T[nx];
		T* B = new T[nx];
		T* C = new T[nx];
		T* F = new T[nx];

		T t = 0;
		while (t <= tmax + EPS) {

			out << std::setprecision(13) << t;
			for (size_t i = 0; i < nx; i++) {
				out << " " << std::setprecision(13) << y[i];
			}
			out << std::endl;

			T* a = new T[nx];
			a[0] = 0;//добавляем чтобы совпадало с A[1] <-> a[1]
			for (size_t i = 1; i < nx; i++) {//методом трапеций самое то
				a[i] = (K(y[i]) + K(y[i - 1])) / 2.;
				//std::cout << a[i]<<"\n";
			}


			T* wminus = new T[nx];
			T* wplus = new T[nx];
			wminus[0] = 0; wplus[0] = 0;
			for (size_t i = 1; i < nx - 1; i++) {
				wminus[i] = a[i] * (y[i] - y[i - 1]) / h;
				wplus[i] = a[i + 1] * (y[i + 1] - y[i]) / h;//iiiiii????
			//	std::cout << wminus[i] << " " << wplus[i] << "\n";
			}
			//std::cout << "end";
			wminus[nx - 1] = a[nx - 1] * (y[nx - 1] - y[nx - 1 - 1]) / h;

			T* ynew = new T[nx];
			T* y0 = new T[nx];
			for (size_t j = 0; j < nx; j++) {
				ynew[j] = y[j];
				y0[j] = y[j];
			}

			T kappa = (sigma * a[1] / h) / (c * ro * h / (2 * tau) + sigma * a[1] / h);
			T mu = (c * ro * y[0] * h / (2 * tau) + sigma * P(t + tau, t0, Q) + (1 - sigma) * (P(t, t0, Q) - wminus[1])) / (c * ro * h / (2 * tau) + sigma * a[1] / h);

			for (size_t k = 0; k < 2; k++) {
				A[nx - 1] = 0; B[nx - 1] = 1; C[nx - 1] = 0; F[nx - 1] = ubound(t + tau);

				for (size_t i = 1; i < nx - 1; i++) {
					A[i] = -(K(y0[i]) + K(y0[i - 1])) / 2 / h / h;
					B[i] = c * ro / tau + (K(y0[i + 1]) + K(y0[i]) + K(y0[i]) + K(y0[i - 1])) / 2 / h / h;
					C[i] = -(K(y0[i + 1]) + K(y0[i])) / 2 / h / h;
					F[i] = c * ro * y[i] / tau;
					//	std::cout << A[i]<<" "<<B[i]<<" "<<C[i]<<" "<<F[i]<<"\n";
				}//опечатка в методичке..
				//std::cout << "kdfjdkfjdkfjkdfjd\n";
				A[0] = -kappa; B[0] = 1; C[0] = 0; F[0] = mu;
				//std::cout << kappa << " "<<mu<<"\n";
				//std::cout << (P(t, t0, Q) - wminus[nx - 1]) <<" ";

				//здесь надо чтобы рассматривал какие поставлены условия

				solve_3d(nx, A, B, C, F, ynew);
				//y1[nx - 1] = kappa * y1[nx - 2] + mu;
				for (size_t j = 0; j < nx; j++) {
					y0[j] = ynew[j];
				}
			}

			for (size_t j = 0; j < nx; j++) {
				y1[j] = ynew[j];
			}
			delete[] y0;
			delete[] ynew;


			for (size_t i = 0; i < nx; i++) {
				y[i] = y1[i];
			}

			t += tau;
			delete[] a;
		}

		out.close();

		out.open("D:\\СЕТКА.txt");
		x = 0;
		for (size_t i = 0; i < nx; i++) {
			out << std::setprecision(13) << x << std::endl;
			x += h;
		}

		delete[] A; delete[] B; delete[] C; delete[] F;

	}

}

void DSolveHeatEquation(T L, T tmax, T c, T ro, T(*K)(T), THERMAL_CONDUCTIVITY_K tc, T(*u0)(T), T(*P1)(T, T, T), T(*P2)(T, T, T), T t0, T Q) {


	T tau = 1e-3; T h = 1e-3;				//шаги сетки
	size_t nx = size_t(round(L / h)) + 1;	//количество узлов
	size_t nt = size_t(round(tmax / tau)) + 1;

	T* y = new T[nx];						//y^j
	T* y1 = new T[nx];					//y^{j+1}

	T sigma = 0.5;

	if (tc == TC_X) {

		std::ofstream out;				// поток для записи
		out.open("D:\\УРАВНЕНИЕ ТЕПЛОПРОВОДНОСТИ 3.txt");

		T x = 0;
		for (size_t i = 0; i < nx; i++) {
			y[i] = u0(x);
			x += h;
			//if (x > 1) std::cout << x << " ";
		}



		T* a = new T[nx];
		a[0] = 0;//добавляем чтобы совпадало с A[1] <-> a[1]
		x = 0;
		for (size_t i = 1; i < nx; i++) {//методом трапеций самое то
			size_t k = 20;
			T dh = h / k;
			T I = 0;
			for (size_t j = 0; j < k; j++) {
				I += 1 / K(x + dh / 2);
				x += dh;
			}
			I *= dh;
			//считаем интеграл => ai
			//не зависят от t?
			a[i] = h / I;
			//std::cout << a[i]<<"\n";
			x += h;
		}



		T* A = new T[nx];
		T* B = new T[nx];
		T* C = new T[nx];
		T* F = new T[nx];

		T t = 0;
		while (t <= tmax + EPS) {

			out << std::setprecision(13) << t;
			for (size_t i = 0; i < nx; i++) {
				out << " " << std::setprecision(13) << y[i];
			}
			out << std::endl;


			T* wminus = new T[nx];
			T* wplus = new T[nx];
			wminus[0] = 0; wplus[0] = 0;
			for (size_t i = 1; i < nx - 1; i++) {
				wminus[i] = a[i] * (y[i] - y[i - 1]) / h;
				wplus[i] = a[i + 1] * (y[i + 1] - y[i]) / h;//a[i+1]????
			//	std::cout << wminus[i] << " " << wplus[i] << "\n";
			}
			//std::cout << "end";
			wminus[nx - 1] = a[nx - 1] * (y[nx - 1] - y[nx - 1 - 1]) / h;


			T kappa1 = (sigma * a[1] / h) / (c * ro * h / (2 * tau) + sigma * a[1] / h);
			T mu1 = (c * ro * y[0] * h / (2 * tau) + sigma * P1(t + tau, t0, Q) + (1 - sigma) * (P1(t, t0, Q) + wplus[1])) / (c * ro * h / (2 * tau) + sigma * a[1] / h);

			T kappa2 = (sigma * a[nx - 1] / h) / (c * ro * h / (2 * tau) + sigma * a[nx - 1] / h);
			T mu2 = (c * ro * y[nx - 1] * h / (2 * tau) + sigma * P2(t + tau, t0, Q) + (1 - sigma) * (P2(t, t0, Q) - wminus[nx - 1])) / (c * ro * h / (2 * tau) + sigma * a[nx - 1] / h);


			//тут меняется только эта строчка но мне лень следовать общепринятым духовным ценностям в программировании
			A[nx - 1] = -kappa2; C[nx - 1] = 1; B[nx - 1] = 0; F[nx - 1] = mu2;
			for (size_t i = 1; i < nx - 1; i++) {
				A[i] = sigma / h * a[i];
				B[i] = sigma / h * a[i + 1];
				C[i] = -(sigma / h * (a[i] + a[i + 1]) + c * ro * h / tau);
				F[i] = -(c * ro * h / tau * y[i] + (1 - sigma) * (wplus[i] - wminus[i]));
				//std::cout << A[i]<<" "<<B[i]<<" "<<C[i]<<" "<<F[i]<<"\n";
			}//опечатка в методичке..
			//std::cout << "kdfjdkfjdkfjkdfjd\n";
			//A[0] = 0; C[0] = kappa; B[0] =1 ; F[0] = mu;
			A[0] = 0; C[0] = 1; B[0] = -kappa1; F[0] = mu1;
			//std::cout << kappa << " "<<mu<<"\n";
			//std::cout << (P(t, t0, Q) - wminus[nx - 1]) <<" ";

			//здесь надо чтобы рассматривал какие поставлены условия

			solve_3d(nx, A, C, B, F, y1);
			//y1[nx - 1] = kappa * y1[nx - 2] + mu;

			//std::cout << y[nx - 1] << "\n";
			for (size_t i = 0; i < nx; i++) {
				y[i] = y1[i];
			}

			std::cout<<integral(y, nx, 0, 1)<<" \n";


			t += tau;
		}

		out.close();

		out.open("D:\\СЕТКА.txt");
		x = 0;
		for (size_t i = 0; i < nx; i++) {
			out << std::setprecision(13) << x << std::endl;
			x += h;
		}

		delete[] A; delete[] B; delete[] C; delete[] F;
		delete[] a;


	}
	else {
		sigma = 0.5;
		std::ofstream out;				// поток для записи
		out.open("D:\\УРАВНЕНИЕ ТЕПЛОПРОВОДНОСТИ КВАЗИЛИНЕЙНОЕ 3.txt");

		T x = 0;
		for (size_t i = 0; i < nx; i++) {
			y[i] = u0(x);
			x += h;
			//if (x > 1) std::cout << x << " ";
		}

		T* A = new T[nx];
		T* B = new T[nx];
		T* C = new T[nx];
		T* F = new T[nx];

		T t = 0;
		while (t <= tmax + EPS) {

			out << std::setprecision(13) << t;
			for (size_t i = 0; i < nx; i++) {
				out << " " << std::setprecision(13) << y[i];
			}
			out << std::endl;

			T* a = new T[nx];
			a[0] = 0;//добавляем чтобы совпадало с A[1] <-> a[1]
			for (size_t i = 1; i < nx; i++) {//методом трапеций самое то
				a[i] = (K(y[i]) + K(y[i - 1])) / 2.;
				//std::cout << a[i]<<"\n";
			}


			T* wminus = new T[nx];
			T* wplus = new T[nx];
			wminus[0] = 0; wplus[0] = 0;
			for (size_t i = 1; i < nx - 1; i++) {
				wminus[i] = a[i] * (y[i] - y[i - 1]) / h;
				wplus[i] = a[i + 1] * (y[i + 1] - y[i]) / h;//iiiiii????
			//	std::cout << wminus[i] << " " << wplus[i] << "\n";
			}
			//std::cout << "end";
			wminus[nx - 1] = a[nx - 1] * (y[nx - 1] - y[nx - 1 - 1]) / h;

			T* ynew = new T[nx];
			T* y0 = new T[nx];
			for (size_t j = 0; j < nx; j++) {
				ynew[j] = y[j];
				y0[j] = y[j];
			}

			T kappa1 = (sigma * a[1] / h) / (c * ro * h / (2 * tau) + sigma * a[1] / h);
			T mu1 = (c * ro * y[0] * h / (2 * tau) + sigma * P1(t + tau, t0, Q) + (1 - sigma) * (P1(t, t0, Q) - wminus[1])) / (c * ro * h / (2 * tau) + sigma * a[1] / h);
			T kappa2 = (sigma * a[nx - 1] / h) / (c * ro * h / (2 * tau) + sigma * a[nx - 1] / h);
			T mu2 = (c * ro * y[nx - 1] * h / (2 * tau) + sigma * P2(t + tau, t0, Q) + (1 - sigma) * (P2(t, t0, Q) - wminus[nx - 1])) / (c * ro * h / (2 * tau) + sigma * a[nx - 1] / h);

			for (size_t k = 0; k < 2; k++) {
				A[nx - 1] = -kappa2; B[nx - 1] = 1; C[nx - 1] = 0; F[nx - 1] = mu2;

				for (size_t i = 1; i < nx - 1; i++) {
					A[i] = -(K(y0[i]) + K(y0[i - 1])) / 2 / h / h;
					B[i] = c * ro / tau + (K(y0[i + 1]) + K(y0[i]) + K(y0[i]) + K(y0[i - 1])) / 2 / h / h;
					C[i] = -(K(y0[i + 1]) + K(y0[i])) / 2 / h / h;
					F[i] = c * ro * y[i] / tau;
					//	std::cout << A[i]<<" "<<B[i]<<" "<<C[i]<<" "<<F[i]<<"\n";
				}//опечатка в методичке..
				//std::cout << "kdfjdkfjdkfjkdfjd\n";
				A[0] = -kappa1; B[0] = 1; C[0] = 0; F[0] = mu1;
				//std::cout << kappa << " "<<mu<<"\n";
				//std::cout << (P(t, t0, Q) - wminus[nx - 1]) <<" ";

				//здесь надо чтобы рассматривал какие поставлены условия

				solve_3d(nx, A, B, C, F, ynew);
				//y1[nx - 1] = kappa * y1[nx - 2] + mu;
				for (size_t j = 0; j < nx; j++) {
					y0[j] = ynew[j];
				}
			}

			for (size_t j = 0; j < nx; j++) {
				y1[j] = ynew[j];
			}
			delete[] y0;
			delete[] ynew;


			for (size_t i = 0; i < nx; i++) {
				y[i] = y1[i];
			}

			t += tau;
			delete[] a;
		}

		out.close();

		out.open("D:\\СЕТКА.txt");
		x = 0;
		for (size_t i = 0; i < nx; i++) {
			out << std::setprecision(13) << x << std::endl;
			x += h;
		}

		delete[] A; delete[] B; delete[] C; delete[] F;

	}

}


void table(T L, T tmax, T c, T ro, T(*K)(T), THERMAL_CONDUCTIVITY_K tc, T(*u0)(T), T(*ubound1)(T), T(*ubound2)(T)) {

	 T h0 = 1e-1;		T tau0 = h0*h0/15;		//шаги сетки
	size_t nx = size_t(round(L / h0)) + 1;	//количество узлов
	size_t nt = size_t(round(tmax / tau0)) + 1;



	T sigma = 0;

	std::ofstream out;				// поток для записи
	out.open("D:\\УРАВНЕНИЕ ТЕПЛОПРОВОДНОСТИ 1.txt");

	std::cout << "         n   &        Шаг сетки tau_n       &        Норма невязки      &       отношение ошибок     &      Порядок сходимости p_n    \n\n";
	T normLast;
	for (size_t n = 1; n <= 5; n++) {
		T tau = tau0 / pow(4, n - 1); h0 = sqrt(15 * tau);
		nx = size_t(round(L / h0)) + 1;	//количество узлов
		nt = size_t(round(tmax / tau)) + 1;
		T* y = new T[nx];						//y^j
		T* y1 = new T[nx];					//y^{j+1}
		std::cout << std::setw(12) << n << " & " << '\t';
		T x = 0;
		T* yTrue = new T[nx];
		for (size_t i = 0; i < nx; i++) {
			yTrue[i] =4* exp(-7*PI*PI*tmax)*sin(PI*x)+2*x+3 ;
			x += h0;
		}


		x = 0;
		for (size_t i = 0; i < nx; i++) {
			y[i] = u0(x);
			x += h0;
			//if (x > 1) std::cout << x << " ";
		}



		T* a = new T[nx];
		a[0] = 0;//добавляем чтобы совпадало с A[1] <-> a[1]
		x = h0;
		for (size_t i = 1; i < nx; i++) {//методом трапеций самое то
			//size_t k = 20;
			//T dh = h0 / k;
			//T I = 0;
			//for (size_t j = 0; j < k; j++) {
			//	I += 1 / K(x + dh / 2);
			//	x += dh;
			//}
			//I *= dh;
			//считаем интеграл => ai
			//не зависят от t?
			//a[i] = h0 / I;
			//std::cout << a[i]<<"\n";
			a[i] = 0.5*(K(x) + K(x - h0));
			x += h0;
			
		}



		T* A = new T[nx];
		T* B = new T[nx];
		T* C = new T[nx];
		T* F = new T[nx];

		T t = 0;
		while (t <= tmax + EPS) {
			if (n == 5) {
				out << std::setprecision(13) << t;
				for (size_t i = 0; i < nx; i++) {
					out << " " << std::setprecision(13) << y[i];
				}
				out << std::endl;
			}



			T* wminus = new T[nx];
			T* wplus = new T[nx];
			wminus[0] = 0; wplus[0] = 0;
			for (size_t i = 1; i < nx - 1; i++) {
				wminus[i] = a[i] * (y[i] - y[i - 1]) / h0;
				wplus[i] = a[i + 1] * (y[i + 1] - y[i]) / h0;//a[i+1]????
			//	std::cout << wminus[i] << " " << wplus[i] << "\n";
			}
			//std::cout << "end";
			wminus[nx - 1] = a[nx - 1] * (y[nx - 1] - y[nx - 1 - 1]) / h0;

			if (sigma > EPS) {
				//тут меняется только эта строчка но мне лень следовать общепринятым духовным ценностям в программировании
				A[nx - 1] = 0; C[nx - 1] = 1; B[nx - 1] = 0; F[nx - 1] = ubound2(t + tau);
				for (size_t i = 1; i < nx - 1; i++) {
					A[i] = sigma / h0 * a[i];
					B[i] = sigma / h0 * a[i + 1];
					C[i] = -(sigma / h0 * (a[i] + a[i + 1]) + c * ro * h0 / tau);
					F[i] = -(c * ro * h0 / tau * y[i] + (1 - sigma) * (wplus[i] - wminus[i]));
					//std::cout << A[i] << " " << B[i] << " " << C[i] << " " << F[i] << "\n";
				}//опечатка в методичке..
				//std::cout << "kdfjdkfjdkfjkdfjd\n";
				//A[0] = 0; C[0] = kappa; B[0] =1 ; F[0] = mu;
				A[0] = 0; C[0] = 1; B[0] = 0; F[0] = ubound1(t + tau);
				//std::cout << kappa << " "<<mu<<"\n";
				//std::cout << (P(t, t0, Q) - wminus[nx - 1]) <<" ";

				//здесь надо чтобы рассматривал какие поставлены условия

				solve_3d(nx, A, C, B, F, y1);
			}
			else {
				for (size_t i = 1; i < nx-1; i++) {
					C[i] = -(sigma / h0 * (a[i] + a[i + 1]) + c * ro * h0 / tau);
					F[i] = -(c * ro * h0 / tau * y[i] + (wplus[i] - wminus[i]));
					y1[i] = y[i] + (wplus[i] - wminus[i]) / (c * ro * h0 / tau);
				//	std::cout << y1[i]<<"\n";
				}
				y1[0]= ubound1(t + tau);
				y1[nx-1]=ubound2(t + tau);
			//std::cout << "fdfdf\n";
			}
			
			//y1[nx - 1] = kappa * y1[nx - 2] + mu;

			//std::cout << y[nx - 1] << "\n";
			for (size_t i = 0; i < nx; i++) {
				y[i] = y1[i];
			//	std::cout << y[i] << "\n";
			}

			t += tau;
		}

		delete[] A; delete[] B; delete[] C; delete[] F;
		delete[] a;




		T* dy = new T[nx];
		for (size_t i = 0; i < nx; i++) {
			dy[i] = y[i] - yTrue[i];
			//	std::cout << "y=" << y[i]<< "u=" << yTrue[i]<<"\n";
		}
		T normNew = normInf(nx,dy);
		
		delete[] dy;
		T z;
		T pn;
		//считаем и выводим

		std::cout << std::setw(20) << tau << " & " << "\t\t";
		if (n == 1) {
			normLast = normNew;
			std::cout << std::setw(20) << normNew << " &    " << "         ---            " << "   &      " << "             ---      " << "\\\\\ \\hline\n";

		}
		else {
			z = normNew / normLast;
			pn = log10(z) / log10(0.5);
			std::cout << std::setw(20) << normNew << " & " << "\t\t";
			std::cout << std::setw(20) << z << " & " << "\t\t";
			std::cout << std::setw(17) << pn; std::cout << "\\\\\ \\hline\n";
			normLast = normNew;
		}

		delete[] yTrue;



	}
	out.close();
	out.open("D:\\СЕТКА.txt");
	T x = 0;
	for (size_t i = 0; i < nx; i++) {
		out << std::setprecision(13) << x << std::endl;
		x += h0;
	}


}


void table2(T L, T tmax, T c, T ro, T(*K)(T), THERMAL_CONDUCTIVITY_K tc, T(*u0)(T), T(*ubound1)(T), T(*ubound2)(T)) {

	T h = 0.2;		T tau0 = 10e-4;		//шаги сетки
	size_t nx = size_t(round(L / h)) + 1;	//количество узлов
	size_t nt = size_t(round(tmax / tau0)) + 1;



	T sigma = 0.5;

	std::ofstream out;				// поток для записи
	out.open("D:\\УРАВНЕНИЕ ТЕПЛОПРОВОДНОСТИ 1.txt");

	std::cout << "         n   &        Шаг сетки tau_n       &        Норма невязки      &       отношение ошибок     &      Порядок сходимости p_n    \n\n";
	T normLast;
	for (size_t n = 1; n <= 5; n++) {
		T tau = tau0 / pow(2, n - 1); T h0 = h / pow(2, n - 1);
		nx = size_t(round(L / h0)) + 1;	//количество узлов
		nt = size_t(round(tmax / tau)) + 1;
		T* y = new T[nx];						//y^j
		T* y1 = new T[nx];					//y^{j+1}
		std::cout << std::setw(12) << n << " & " << '\t';
		T x = 0;
		T* yTrue = new T[nx];
		for (size_t i = 0; i < nx; i++) {
			yTrue[i] = 4 * exp(-7 * PI * PI * tmax) * sin(PI * x) + 2 * x + 3;
			x += h0;
		}


		x = 0;
		for (size_t i = 0; i < nx; i++) {
			y[i] = u0(x);
			x += h0;
			//if (x > 1) std::cout << x << " ";
		}



		T* a = new T[nx];
		a[0] = 0;//добавляем чтобы совпадало с A[1] <-> a[1]
		x = h0;
		for (size_t i = 1; i < nx; i++) {//методом трапеций самое то
			//size_t k = 20;
			//T dh = h0 / k;
			//T I = 0;
			//for (size_t j = 0; j < k; j++) {
			//	I += 1 / K(x + dh / 2);
			//	x += dh;
			//}
			//I *= dh;
			////считаем интеграл => ai
			////не зависят от t?
			//a[i] = h0 / I;
			//std::cout << a[i]<<"\n";
			//a[i] = 0.5 * (K(x) + K(x - h0));
			a[i] = K(x);
			x += h0;

		}



		T* A = new T[nx];
		T* B = new T[nx];
		T* C = new T[nx];
		T* F = new T[nx];

		T t = 0;
		while (t <= tmax + EPS) {
			if (n == 5) {
				out << std::setprecision(13) << t;
				for (size_t i = 0; i < nx; i++) {
					out << " " << std::setprecision(13) << y[i];
				}
				out << std::endl;
			}



			T* wminus = new T[nx];
			T* wplus = new T[nx];
			wminus[0] = 0; wplus[0] = 0;
			for (size_t i = 1; i < nx - 1; i++) {
				wminus[i] = a[i] * (y[i] - y[i - 1]) / h0;
				wplus[i] = a[i + 1] * (y[i + 1] - y[i]) / h0;//a[i+1]????
			//	std::cout << wminus[i] << " " << wplus[i] << "\n";
			}
			//std::cout << "end";
			wminus[nx - 1] = a[nx - 1] * (y[nx - 1] - y[nx - 1 - 1]) / h0;

			if (sigma > EPS) {
				//тут меняется только эта строчка но мне лень следовать общепринятым духовным ценностям в программировании
				A[nx - 1] = 0; C[nx - 1] = 1; B[nx - 1] = 0; F[nx - 1] = ubound2(t + tau);
				for (size_t i = 1; i < nx - 1; i++) {
					A[i] = sigma / h0 * a[i];
					B[i] = sigma / h0 * a[i + 1];
					C[i] = -(sigma / h0 * (a[i] + a[i + 1]) + c * ro * h0 / tau);
					F[i] = -(c * ro * h0 / tau * y[i] + (1 - sigma) * (wplus[i] - wminus[i]));
					//std::cout << A[i] << " " << B[i] << " " << C[i] << " " << F[i] << "\n";
				}//опечатка в методичке..
				//std::cout << "kdfjdkfjdkfjkdfjd\n";
				//A[0] = 0; C[0] = kappa; B[0] =1 ; F[0] = mu;
				A[0] = 0; C[0] = 1; B[0] = 0; F[0] = ubound1(t + tau);
				//std::cout << kappa << " "<<mu<<"\n";
				//std::cout << (P(t, t0, Q) - wminus[nx - 1]) <<" ";

				//здесь надо чтобы рассматривал какие поставлены условия

				solve_3d(nx, A, C, B, F, y1);
				y1[0] = ubound1(t + tau);
				y1[nx - 1] = ubound2(t + tau);
			}
			else {
				for (size_t i = 1; i < nx - 1; i++) {
					C[i] = -(sigma / h0 * (a[i] + a[i + 1]) + c * ro * h0 / tau);
					F[i] = -(c * ro * h0 / tau * y[i] + (wplus[i] - wminus[i]));
					y1[i] = y[i] + (wplus[i] - wminus[i]) / (c * ro * h0 / tau);
					//	std::cout << y1[i]<<"\n";
				}
				y1[0] = ubound1(t + tau);
				y1[nx - 1] = ubound2(t + tau);
				//std::cout << "fdfdf\n";
			}

			//y1[nx - 1] = kappa * y1[nx - 2] + mu;

			//std::cout << y[nx - 1] << "\n";
			for (size_t i = 0; i < nx; i++) {
				y[i] = y1[i];
				//	std::cout << y[i] << "\n";
			}

			t += tau;
		}

		delete[] A; delete[] B; delete[] C; delete[] F;
		delete[] a;




		T* dy = new T[nx];
		T fuck = 0;
		for (size_t i = 0; i < nx; i++) {
			dy[i] = y[i] - yTrue[i];
			if (abs(dy[i]) > fuck) {
				fuck = abs(dy[i]);
			}
			//	std::cout << "y=" << y[i]<< "u=" << yTrue[i]<<"\n";
		}
		T normNew = fuck; normInf(nx, dy);

		delete[] dy;
		T z;
		T pn;
		//считаем и выводим

		std::cout << std::setw(20) << tau << " & " << "\t\t";
		if (n == 1) {
			normLast = normNew;
			std::cout << std::setw(20) << normNew << " &    " << "         ---            " << "   &      " << "             ---      " << "\\\\\ \\hline\n";

		}
		else {
			z = normNew / normLast;
			pn = log10(z) / log10(0.5);
			std::cout << std::setw(20) << normNew << " & " << "\t\t";
			std::cout << std::setw(20) << z << " & " << "\t\t";
			std::cout << std::setw(17) << pn; std::cout << "\\\\\ \\hline\n";
			normLast = normNew;
		}

		delete[] yTrue;



	}
	out.close();
	out.open("D:\\СЕТКА.txt");
	T x = 0;
	for (size_t i = 0; i < nx; i++) {
		out << std::setprecision(13) << x << std::endl;
		x += h;
	}


}


//void DSolveHeatEquation(T L, T tmax, T c, T ro, T(*K)(T), THERMAL_CONDUCTIVITY_K tc, T(*u0)(T), T(*P1)(T, T, T), T t0, T Q, T(*P2)(T, T, T)) {
//
//	T tau = 1e-3; T h = 1e-3;				//шаги сетки
//	size_t nx = size_t(round(L / h)) + 1;	//количество узлов
//	size_t nt = size_t(round(tmax / tau)) + 1;
//
//	T* y = new T[nx];						//y^j
//	T* y1 = new T[nx];					//y^{j+1}
//
//	T sigma = 0.5;
//
//	if (tc == TC_X) {
//
//		std::ofstream out;				// поток для записи
//		out.open("D:\\УРАВНЕНИЕ ТЕПЛОПРОВОДНОСТИ.txt");
//
//		T x = 0;
//		for (size_t i = 0; i < nx; i++) {
//			y[i] = u0(x);
//			x += h;
//			//if (x > 1) std::cout << x << " ";
//		}
//
//
//
//		T* a = new T[nx];
//		a[0] = 0;//добавляем чтобы совпадало с A[1] <-> a[1]
//		x = 0;
//		for (size_t i = 1; i < nx; i++) {//методом трапеций самое то
//			size_t k = 20;
//			T dh = h / k;
//			T I = 0;
//			for (size_t j = 0; j < k; j++) {
//				I += 1 / K(x + dh / 2);
//				x += dh;
//			}
//			I *= dh;
//			//считаем интеграл => ai
//			//не зависят от t?
//			a[i] = h / I;
//			//std::cout << a[i]<<"\n";
//			x += h;
//		}
//
//
//
//		T* A = new T[nx];
//		T* B = new T[nx];
//		T* C = new T[nx];
//		T* F = new T[nx];
//
//		T t = 0;
//		while (t <= tmax + EPS) {
//
//			out << std::setprecision(13) << t;
//			for (size_t i = 0; i < nx; i++) {
//				out << " " << std::setprecision(13) << y[i];
//			}
//			out << std::endl;
//
//
//			T* wminus = new T[nx];
//			T* wplus = new T[nx];
//			wminus[0] = 0; wplus[0] = 0;
//			for (size_t i = 1; i < nx - 1; i++) {
//				wminus[i] = a[i] * (y[i] - y[i - 1]) / h;
//				wplus[i] = a[i + 1] * (y[i + 1] - y[i]) / h;//a[i+1]????
//			//	std::cout << wminus[i] << " " << wplus[i] << "\n";
//			}
//			//std::cout << "end";
//			wminus[nx - 1] = a[nx - 1] * (y[nx - 1] - y[nx - 1 - 1]) / h;
//
//
//			T kappa = (sigma * a[nx - 1] / h) / (c * ro * h / (2 * tau) + sigma * a[nx - 1] / h);
//			T mu = (c * ro * y[nx - 1] * h / (2 * tau) + sigma * P(t + tau, t0, Q) + (1 - sigma) * (P(t, t0, Q) - wminus[nx - 1])) / (c * ro * h / (2 * tau) + sigma * a[nx - 1] / h);
//
//
//
//			A[0] = 0; C[0] = 1; B[0] = 0; F[0] = ubound(t + tau);
//			for (size_t i = 1; i < nx - 1; i++) {
//				A[i] = sigma / h * a[i];
//				B[i] = sigma / h * a[i + 1];
//				C[i] = -(sigma / h * (a[i] + a[i + 1]) + c * ro * h / tau);
//				F[i] = -(c * ro * h / tau * y[i] + (1 - sigma) * (wplus[i] - wminus[i]));
//				//std::cout << A[i]<<" "<<B[i]<<" "<<C[i]<<" "<<F[i]<<"\n";
//			}//опечатка в методичке..
//			//std::cout << "kdfjdkfjdkfjkdfjd\n";
//			A[nx - 1] = -kappa; C[nx - 1] = 1; B[nx - 1] = 0; F[nx - 1] = mu;
//
//			//std::cout << kappa << " "<<mu<<"\n";
//			//std::cout << (P(t, t0, Q) - wminus[nx - 1]) <<" ";
//
//			//здесь надо чтобы рассматривал какие поставлены условия
//
//			solve_3d(nx, A, C, B, F, y1);
//			//y1[nx - 1] = kappa * y1[nx - 2] + mu;
//
//
//			for (size_t i = 0; i < nx; i++) {
//				y[i] = y1[i];
//			}
//
//			t += tau;
//		}
//
//		out.close();
//
//		out.open("D:\\СЕТКА.txt");
//		x = 0;
//		for (size_t i = 0; i < nx; i++) {
//			out << std::setprecision(13) << x << std::endl;
//			x += h;
//		}
//
//		delete[] A; delete[] B; delete[] C; delete[] F;
//		delete[] a;
//
//
//	}
//	else {
//		sigma = 1;
//		std::ofstream out;				// поток для записи
//		out.open("D:\\УРАВНЕНИЕ ТЕПЛОПРОВОДНОСТИ КВАЗИЛИНЕЙНОЕ.txt");
//
//		T x = 0;
//		for (size_t i = 0; i < nx; i++) {
//			y[i] = u0(x);
//			x += h;
//			//if (x > 1) std::cout << x << " ";
//		}
//
//		T* A = new T[nx];
//		T* B = new T[nx];
//		T* C = new T[nx];
//		T* F = new T[nx];
//
//		T t = 0;
//		while (t <= tmax + EPS) {
//
//			out << std::setprecision(13) << t;
//			for (size_t i = 0; i < nx; i++) {
//				out << " " << std::setprecision(13) << y[i];
//			}
//			out << std::endl;
//
//			T* a = new T[nx];
//			a[0] = 0;//добавляем чтобы совпадало с A[1] <-> a[1]
//			for (size_t i = 1; i < nx; i++) {//методом трапеций самое то
//				a[i] = (K(y[i]) + K(y[i - 1])) / 2.;
//				//std::cout << a[i]<<"\n";
//			}
//
//
//			T* wminus = new T[nx];
//			T* wplus = new T[nx];
//			wminus[0] = 0; wplus[0] = 0;
//			for (size_t i = 1; i < nx - 1; i++) {
//				wminus[i] = a[i] * (y[i] - y[i - 1]) / h;
//				wplus[i] = a[i + 1] * (y[i + 1] - y[i]) / h;//iiiiii????
//			//	std::cout << wminus[i] << " " << wplus[i] << "\n";
//			}
//			//std::cout << "end";
//			wminus[nx - 1] = a[nx - 1] * (y[nx - 1] - y[nx - 1 - 1]) / h;
//
//			T* ynew = new T[nx];
//			T* y0 = new T[nx];
//			for (size_t j = 0; j < nx; j++) {
//				ynew[j] = y[j];
//				y0[j] = y[j];
//			}
//
//			for (size_t k = 0; k < 2; k++) {
//				A[0] = 0; B[0] = 1; C[0] = 0; F[0] = ubound(t + tau);
//
//				for (size_t i = 1; i < nx - 1; i++) {
//					A[i] = -(K(y0[i]) + K(y0[i - 1])) / 2 / h / h;
//					B[i] = c * ro / tau + (K(y0[i + 1]) + K(y0[i]) + K(y0[i]) + K(y0[i - 1])) / 2 / h / h;
//					C[i] = -(K(y0[i + 1]) + K(y0[i])) / 2 / h / h;
//					F[i] = c * ro * y[i] / tau;
//					//	std::cout << A[i]<<" "<<B[i]<<" "<<C[i]<<" "<<F[i]<<"\n";
//				}//опечатка в методичке..
//				//std::cout << "kdfjdkfjdkfjkdfjd\n";
//				T kappa = (sigma * a[nx - 1] / h) / (c * ro * h / (2 * tau) + sigma * a[nx - 1] / h);
//				T mu = (c * ro * y[nx - 1] * h / (2 * tau) + sigma * P(t + tau, t0, Q) + (1 - sigma) * (P(t, t0, Q) - wminus[nx - 1])) / (c * ro * h / (2 * tau) + sigma * a[nx - 1] / h);
//				A[nx - 1] = -kappa; B[nx - 1] = 1; C[nx - 1] = 0; F[nx - 1] = mu;
//				//std::cout << kappa << " "<<mu<<"\n";
//				//std::cout << (P(t, t0, Q) - wminus[nx - 1]) <<" ";
//
//				//здесь надо чтобы рассматривал какие поставлены условия
//
//				solve_3d(nx, A, B, C, F, ynew);
//				//y1[nx - 1] = kappa * y1[nx - 2] + mu;
//				for (size_t j = 0; j < nx; j++) {
//					y0[j] = ynew[j];
//				}
//			}
//
//			for (size_t j = 0; j < nx; j++) {
//				y1[j] = ynew[j];
//			}
//			delete[] y0;
//			delete[] ynew;
//
//
//			for (size_t i = 0; i < nx; i++) {
//				y[i] = y1[i];
//			}
//
//			t += tau;
//			delete[] a;
//		}
//
//		out.close();
//
//		out.open("D:\\СЕТКА.txt");
//		x = 0;
//		for (size_t i = 0; i < nx; i++) {
//			out << std::setprecision(13) << x << std::endl;
//			x += h;
//		}
//
//		delete[] A; delete[] B; delete[] C; delete[] F;
//
//	}
//
//}

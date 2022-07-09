#include"header.h"

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


void DSolveWaveEquation(T L, T tmax, T a2, T(*f)(T), T(*ddf)(T), T(*g)(T), T(*phi)(T), T(*psi)(T), DIFFERENTIAL_TYPE dt) {

	T h = 0.005;	T tau = h;				//шаги сетки
	size_t nx = size_t(round(L / h)) + 1;		//количество узлов
	size_t nt = size_t(round(tmax / tau)) + 1;

	T* y1 = new T[nx];							//y^{j-1}
	T* y2 = new T[nx];							//y^{j}
	T* y3 = new T[nx];							//y^{j+1}

	std::ofstream out;							// поток для записи
	out.open("D:\\ВОЛНОВОЕ УРАВНЕНИЕ 1.txt");
	std::ofstream outE;
	outE.open("D:\\ЭНЕРГИЯ.txt");


	T x = 0;
	for (size_t i = 0; i < nx; i++) {
		y1[i] = f(x);
		x += h;
	}
	out << std::setprecision(13) << 0;
	for (size_t i = 0; i < nx; i++) {
		out << " " << std::setprecision(13) << y1[i];
	}
	out << std::endl;
	T dE = 0;
	outE << std::setprecision(13) << 0 << " " << std::setprecision(13) << dE << std::endl;


	x = h;
	y2[0] = phi(tau); y2[nx-1] = psi(tau);
	if (dt==DT_ANALYTIC) {
		for (size_t i = 1; i < nx - 1; i++) {
			y2[i] = y1[i] + tau * g(x) + a2 * tau * tau / 2 * ddf(x);
			x += h;
		}
	}
	else {
		for (size_t i = 1; i < nx - 1; i++) {
			T fxx = (f(x + h) - 2 * f(x) + f(x - h)) / h / h;
			y2[i] = y1[i] + tau * g(x) + a2 * tau * tau / 2 * fxx;
			x += h;
		}
	}

	out << std::setprecision(13) << tau;
	for (size_t i = 0; i < nx; i++) {
		out << " " << std::setprecision(13) << y2[i];
	}
	out << std::endl;

	bool check = false;
	T t = 2*tau;
	std::cout << "\n0%";
	while (t <= tmax + EPS) {

		y3[0] = phi(t); y3[nx-1] = psi(t);

		for (size_t i = 1; i < nx - 1; i++) {
			y3[i] = pow(tau / h, 2) * a2 * (y2[i + 1] - 2 * y2[i] + y2[i - 1]) - y1[i] + 2 * y2[i];
		}


		T S1 = 0; T S2 = 0;
		for (size_t i = 1; i < nx - 1; i++) {
			S1 += (y2[i] - y1[i]) * (y3[i] - 2 * y2[i] + y1[i]);
			S2 += (y2[i] - y2[i - 1]) * ((y2[i] - y1[i]) - (y2[i - 1] - y1[i - 1]));
		}
		S2 += (y2[nx - 1] - y2[nx - 1 - 1]) * ((y2[nx - 1] - y1[nx - 1]) - (y2[nx - 1 - 1] - y1[nx - 1 - 1]));
		dE = h*S1 / a2 / (pow(tau, 3)) + S2 / h / tau + (-(y2[nx - 1] - y2[nx - 2]) * (y2[nx - 1] - y1[nx - 1]) + (y2[1] - y2[0]) * (y2[0] - y1[0])) / (h * tau);
		//std::cout << " " << /*S1 << " " << S2 << " " <<*/ dE << "\n";
		outE << std::setprecision(13) << t << " " << std::setprecision(13) << dE << std::endl;

		for (size_t i = 0; i < nx; i++) {
			y1[i] = y2[i];
			y2[i] = y3[i];
		}

		out << std::setprecision(13) << t;
		for (size_t i = 0; i < nx; i++) {
			out << " " << std::setprecision(13) << y3[i];
		}
		out << std::endl;

		t += tau;

		if (abs(t-tmax / 2) <tau/2 && check) {
			std::cout << "\n50%";

		}
		else {
			if (abs(t - tmax / 4) < tau / 2) {
				std::cout << "\n25%";
			}
			else {
				if (abs(t - 3*tmax / 4) < tau / 2) {
					std::cout << "\n75%";
				}
			}
		}
	}
	std::cout << "\n100%\n";
	out.close();
	outE.close();

	out.open("D:\\СЕТКА.txt");
	x = 0;
	for (size_t i = 0; i < nx; i++) {
		out << std::setprecision(13) << x << std::endl;
		x += h;
	}

	delete[] y1; delete[] y2; delete[] y3;
}

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
		T* ySol2 = new T[size_t(round(L / (h/2))) + 1];

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
			dy[i] = ySol2[2*i] - ySol1[i];
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

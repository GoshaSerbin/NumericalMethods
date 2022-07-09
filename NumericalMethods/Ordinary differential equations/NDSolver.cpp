#include"header.h"

void NDSolveRungeKutta2(T* (*f)(T, T*, size_t, T*), T* initialConditions, size_t N, T t0, T tmax, T tau0 = TAU0) {

	//Коэффициенты таблицы Бутчера
	T b21 = 0.5;
	T a2 = 0.5;
	T sigma1 = 0;
	T sigma2 = 1;
	

	T* y = new T[N]; //y_{n+1} для шага tau
	T* y_n = new T[N];
	T* y_m = new T[N];//y между y_n и y_{n+1}
	T* k1 = new T[N];
	T* k2 = new T[N];

	T tau = tau0;
	size_t n = 1;
	T t_n = t0;
	for (size_t i = 0; i < N; i++) {
		y_n[i] = initialConditions[i];//y_i(t_0);
	}

	std::ofstream out;          // поток для записи
	out.open("D:\\РУНГЕ-КУТТА2.txt");
	while (t_n <= tmax+EPS) {

		out << std::setprecision(13)<< t_n;
		for (size_t i = 0; i < N; i++) {
			out << " " << std::setprecision(13) << y_n[i];
		}
		out << std::endl;

		{
			f(t_n, y_n, N, k1);
			for (size_t i = 0; i < N; i++) {
				y_m[i] = y_n[i] + tau * b21 * k1[i];
			}
			f(t_n + a2 * tau, y_m, N, k2);

			for (size_t i = 0; i < N; i++) {
				y[i] = y_n[i] + tau * (sigma1 * k1[i] + sigma2 * k2[i]);
			}
		}

		for (size_t i = 0; i < N; i++) {
			y_n[i] = y[i];
		}


		t_n = t_n + tau;
	}

	delete[] y;
	delete[] y_n;
	delete[] y_m;
	delete[] k1; delete[] k2;
}


void NDSolveRungeKutta4(T* (*f)(T, T*, size_t, T*), T* initialConditions, size_t N, T t0, T tmax, T tau0=TAU0) {

	//Коэффициенты таблицы Бутчера
	T b21 = 0.5;
	T b31 = 0; T b32 = 0.5;
	T b41 = 0; T b42 = 0; T b43 = 1;

	T a2 = 0.5; T a3 = 0.5; T a4 = 1;
	T s1 = 1. / 6; T s2 = 2. / 6; T s3 = 2. / 6; T s4 = 1. / 6;

	T* y = new T[N]; //y_{n+1}
	T* y_n = new T[N];
	T* y_n2 = new T[N];
	T* y_n3 = new T[N];
	T* y_n4 = new T[N];

	T tau = tau0;
	size_t n = 1;
	T t_n = t0;
	for (size_t i = 0; i < N; i++) {
		y_n[i] = initialConditions[i];//y_i(t_0);
	}

	T* k1 = new T[N];//вектор для сохранения промежуточных значений
	T* k2 = new T[N];
	T* k3 = new T[N];
	T* k4 = new T[N];

	std::ofstream out;          // поток для записи
	out.open("D:\\РУНГЕ-КУТТА4.txt");

	while (t_n <= tmax + EPS) {

		out << std::setprecision(13) << t_n;
		for (size_t i = 0; i < N; i++) {
			out << " " << std::setprecision(13) << y_n[i];
		}
		out << std::endl;

		f(t_n, y_n, N, k1);

		for (size_t i = 0; i < N; i++) {
			y_n2[i] = y_n[i] + tau * b21 * k1[i];
		}
		f(t_n + a2 * tau, y_n2, N, k2);

		for (size_t i = 0; i < N; i++) {
			y_n3[i] = y_n[i] + tau * (b31 * k1[i] + b32 * k2[i]);
		}
		f(t_n + a3 * tau, y_n3, N, k3);


		for (size_t i = 0; i < N; i++) {
			y_n4[i] = y_n[i] + tau * (b41 * k1[i] + b42 * k2[i] + b43 * k3[i]);
		}
		f(t_n + a4 * tau, y_n4, N, k4);

		for (size_t i = 0; i < N; i++) {
			y[i] = y_n[i] + tau * (s1 * k1[i] + s2 * k2[i] + s3 * k3[i] + s4 * k4[i]);
		}

		t_n+=tau;
		n++;

		for (size_t i = 0; i < N; i++) {
			y_n[i] = y[i];
		}
	}


	delete[] y;
	delete[] y_n;
	delete[] y_n2; delete[] y_n3; delete[] y_n4;
}

void NDSolveRungeKutta4_M(T* (*f)(T, T*, size_t, T*), T* initialConditions, size_t N, T t0, T tmax, T tau0 = TAU0) {

	T facmin = 0.55;
	T facmax = 1.3;//Надо менять ниже в коде тоже
	T fac = 0.8;

	//Коэффициенты таблицы Бутчера
	T b21 = 0.5;
	T b31 = 0; T b32 = 0.5;
	T b41 = 0; T b42 = 0; T b43 = 1;

	T a2 = 0.5; T a3 = 0.5; T a4 = 1;
	T s1 = 1./6; T s2 = 2. / 6; T s3 = 2. / 6; T s4 = 1. / 6;

	T p = 4;//порядок сходимости
	T tol = 1e-10;//локальная точность

	T* y2 = new T[N]; //y_{n+1} для шага 2*tau
	T* y11 = new T[N]; T* y12 = new T[N]; //y_{n+1} для шага tau
	T* y_n = new T[N];
	T* y_n2 = new T[N];
	T* y_n3 = new T[N];
	T* y_n4 = new T[N];

	T tau = tau0;
	size_t n = 1;
	T t_n = t0;
	for (size_t i = 0; i < N; i++) {
		y_n[i] = initialConditions[i];//y_i(t_0);
	}
	std::ofstream out;          // поток для записи
	out.open("D:\\РУНГЕ-КУТТА4М.txt");
	out << std::setprecision(13) << t_n;
	for (size_t i = 0; i < N; i++) {
		out << " " << std::setprecision(13) << y_n[i];
	}
	out << std::endl;
	T* k1 = new T[N];
	T* k2 = new T[N];
	T* k3 = new T[N];
	T* k4 = new T[N];


	std::ofstream outtau;          // поток для записи
	outtau.open("D:\\ТАУ.txt");
	while (t_n <= tmax + EPS) {

		outtau << std::setprecision(13)<< t_n <<" "<< std::setprecision(13) << tau <<std::endl;
		T tau2 = 2 * tau;
		
		//2*tau
		{
			f(t_n, y_n, N, k1);

			for (size_t i = 0; i < N; i++) {
				y_n2[i] = y_n[i] + tau2 * b21 * k1[i];
			}
			f(t_n + a2 * tau2, y_n2, N, k2);

			for (size_t i = 0; i < N; i++) {
				y_n3[i] = y_n[i] + tau2 * (b31 * k1[i] + b32 * k2[i]);
			}
			f(t_n + a3 * tau2, y_n3, N, k3);


			for (size_t i = 0; i < N; i++) {
				y_n4[i] = y_n[i] + tau2 * (b41 * k1[i] + b42 * k2[i] + b43 * k3[i]);
			}
			f(t_n + a4 * tau2, y_n4, N, k4);

			for (size_t i = 0; i < N; i++) {
				y2[i] = y_n[i] + tau2 * (s1 * k1[i] + s2 * k2[i] + s3 * k3[i] + s4 * k4[i]);
			}

		}


		//tau
		{
			f(t_n, y_n, N, k1);

			for (size_t i = 0; i < N; i++) {
				y_n2[i] = y_n[i] + tau * b21 * k1[i];
			}
			f(t_n + a2 * tau, y_n2, N, k2);

			for (size_t i = 0; i < N; i++) {
				y_n3[i] = y_n[i] + tau * (b31 * k1[i] + b32 * k2[i]);
			}
			f(t_n + a3 * tau, y_n3, N, k3);


			for (size_t i = 0; i < N; i++) {
				y_n4[i] = y_n[i] + tau * (b41 * k1[i] + b42 * k2[i] + b43 * k3[i]);
			}
			f(t_n + a4 * tau, y_n4, N, k4);

			for (size_t i = 0; i < N; i++) {
				y11[i] = y_n[i] + tau * (s1 * k1[i] + s2 * k2[i] + s3 * k3[i] + s4 * k4[i]);
			}

			T t_n2 = t_n + tau;

			f(t_n2, y11, N, k1);

			for (size_t i = 0; i < N; i++) {
				y_n2[i] = y11[i] + tau * b21 * k1[i];
			}
			f(t_n2 + a2 * tau, y_n2, N, k2);

			for (size_t i = 0; i < N; i++) {
				y_n3[i] = y11[i] + tau * (b31 * k1[i] + b32 * k2[i]);
			}
			f(t_n2 + a3 * tau, y_n3, N, k3);


			for (size_t i = 0; i < N; i++) {
				y_n4[i] = y11[i] + tau * (b41 * k1[i] + b42 * k2[i] + b43 * k3[i]);
			}
			f(t_n2 + a4 * tau, y_n4, N, k4);

			for (size_t i = 0; i < N; i++) {
				y12[i] = y11[i] + tau * (s1 * k1[i] + s2 * k2[i] + s3 * k3[i] + s4 * k4[i]);
			}
		}

		T* dy = new T[N];

		for (size_t i = 0; i < N; i++) {
			dy[i] = y2[i] - y12[i];
		}
		//T err = abs(norm1(dy, N)/(1-pow(2,p)));
		T err = abs(abs((y2[0] - y12[0])) / (1 - pow(2, p)));
	//	std::cout << err <<" ";
	//	std::cout << "r" << err << " ";
		delete[] dy;

		T tau_new;
		if (facmin > fac * pow(tol / err, 1. / (1 + p))) {
			tau_new = facmin * tau;
		}
		else {
			if (facmax < fac * pow(tol / err, 1. / (1 + p))) {
				tau_new = facmax * tau;
			}
			else {
				tau_new = fac * pow(tol / err, 1. / (1 + p)) * tau;
			}
		}
		if (err > tol) {
			//std::cout << "new=" << tau_new<<" ";
			tau = tau_new;
			facmax = 1;
			continue;
		}
		else { facmax = 1.3; }//      ТУТ
		for (size_t i = 0; i < N; i++) {
			y_n[i] = y12[i];
		}
		//std::cout << t_n << " \n";
		t_n += tau;
		out << std::setprecision(13) << t_n;
		for (size_t i = 0; i < N; i++) {
			out << " " << std::setprecision(13) << y11[i];
		}
		out << std::endl;
		t_n += tau;
		out << std::setprecision(13) << t_n;
		for (size_t i = 0; i < N; i++) {
			out << " " << std::setprecision(13) << y12[i];
		}
		out << std::endl;

		tau = tau_new;

		n++;
	}


	delete[] y11; delete[] y12; delete[] y2;
	delete[] y_n;
	delete[] y_n2; delete[] y_n3; delete[] y_n4;
	delete[] k1; delete[] k2; delete[] k3; delete[] k4;
}

void NDSolveExplicitEuler(T* (*f)(T, T*, size_t, T*), T* initialConditions, size_t N, T t0, T tmax, T tau0 = TAU0) {


	T* y = new T[N]; //y_{n+1}
	T* y_n = new T[N];

	T tau = tau0;
	size_t n = 1;
	T t_n = t0;
	for (size_t i = 0; i < N; i++) {
		y_n[i] = initialConditions[i];//y_i(t_0);
	}
	T* k = new T[N];
	std::ofstream out;          // поток для записи
	out.open("D:\\ЯВНЫЙ ЭЙЛЕР.txt");// окрываем файл для записи
	while (t_n <= tmax + EPS) {

		out << std::setprecision(13) << t_n;
		for (size_t i = 0; i < N; i++) {
			out << " " << std::setprecision(13) << y_n[i];
		}
		out << std::endl;

		
		f(t_n, y_n, N, k);
		for (size_t i = 0; i < N; i++) {
			y[i] = y_n[i] + tau * k[i];
		}

		t_n += tau;
		n++;

		for (size_t i = 0; i < N; i++) {
			y_n[i] = y[i];
		}
	}

	delete[] k;
	delete[] y;
	delete[] y_n;
}

void NDSolveImplicitEuler(T* (*f)(T, T*, size_t, T*), T* initialConditions, size_t N, T t0, T tmax, T tau0 = TAU0) {


	T* y = new T[N]; //y_{n+1}
	T* y_n = new T[N];
	T* vec = new T[N];//вектор для сохранения промежуточных значений

	T tau = tau0;
	size_t n = 1;
	T t = t0;
	for (size_t i = 0; i < N; i++) {
		y_n[i] = initialConditions[i];//y_i(t_0);
	}
	T* d1 = new T[N];
	T* d2 = new T[N];
	for (size_t i = 0; i < N; i++) {
		d1[i] = -100;
		d2[i] = 100;
	}

	std::ofstream out;          // поток для записи
	out.open("D:\\НЕЯВНЫЙ ЭЙЛЕР.txt");// окрываем файл для записи
	while (t <= tmax + EPS) {

		out << std::setprecision(13) << t;
		for (size_t i = 0; i < N; i++) {
			out << " " << std::setprecision(13) << y_n[i];
		}
		out << std::endl;

		t += tau;
		n++;

		if (NSolveNewton(systemForImplicitEuler, f, t, y_n, tau, d1, d2, y_n, N, y) == -1) break;

		for (size_t i = 0; i < N; i++) {
			y_n[i] = y[i];
		}
	}


	delete[] y;
	delete[] y_n;
	delete[] vec;
	delete[] d1;
	delete[] d2;
}

void NDSolveSymmetric(T* (*f)(T, T*, size_t, T*), T* initialConditions, size_t N, T t0, T tmax, T tau0 = TAU0) {


	T* y = new T[N]; //y_{n+1}
	T* y_n = new T[N];
	T* vec = new T[N];//вектор для сохранения промежуточных значений

	T tau = tau0;
	size_t n = 1;
	T t_n = t0;
	for (size_t i = 0; i < N; i++) {
		y_n[i] = initialConditions[i];//y_i(t_0);
	}
	T* d1 = new T[N];
	T* d2 = new T[N];
	for (size_t i = 0; i < N; i++) {
		d1[i] = -100;
		d2[i] = 100;
	}

	std::ofstream out;          // поток для записи
	out.open("D:\\СИММЕТРИЧНАЯ СХЕМА.txt");// окрываем файл для записи
	while (t_n <= tmax + EPS) {

		out << std::setprecision(13) << t_n;
		for (size_t i = 0; i < N; i++) {
			out << " " << std::setprecision(13) << y_n[i];
		}
		out << std::endl;


		if (NSolveNewton(systemForSymmetric, f, t_n, y_n, tau, d1, d2, y_n, N, y) == -1) break;

		for (size_t i = 0; i < N; i++) {
			y_n[i] = y[i];
		}


		t_n += tau;
		n++;
	}


	delete[] y;
	delete[] y_n;
	delete[] vec;
	delete[] d1;
	delete[] d2;
}


void initRungeKutta(T* (*f)(T, T*, size_t, T*), T* initialConditions, size_t N, size_t k, T t0, T tau, T** fval) {
	//Коэффициенты таблицы Бутчера
	T b21 = 0.5;
	T b31 = 0; T b32 = 0.5;
	T b41 = 0; T b42 = 0; T b43 = 1;

	T a2 = 0.5; T a3 = 0.5; T a4 = 1;
	T s1 = 1. / 6; T s2 = 2. / 6; T s3 = 2. / 6; T s4 = 1. / 6;

	T* y = new T[N]; //y_{n+1}
	T* y_n = new T[N];
	T* y_n2 = new T[N];
	T* y_n3 = new T[N];
	T* y_n4 = new T[N];

	size_t n = 0;
	T t_n = t0;
	for (size_t i = 0; i < N; i++) {
		y_n[i] = initialConditions[i];//y_i(t_0);
	}

	for (size_t i = 0; i < N; i++) {
		fval[0][i] = y_n[i];
	}

	T* k1 = new T[N];//вектор для сохранения промежуточных значений
	T* k2 = new T[N];
	T* k3 = new T[N];
	T* k4 = new T[N];

	while (n< k) {

		f(t_n, y_n, N, k1);

		for (size_t i = 0; i < N; i++) {
			y_n2[i] = y_n[i] + tau * b21 * k1[i];
		}
		f(t_n + a2 * tau, y_n2, N, k2);

		for (size_t i = 0; i < N; i++) {
			y_n3[i] = y_n[i] + tau * (b31 * k1[i] + b32 * k2[i]);
		}
		f(t_n + a3 * tau, y_n3, N, k3);


		for (size_t i = 0; i < N; i++) {
			y_n4[i] = y_n[i] + tau * (b41 * k1[i] + b42 * k2[i] + b43 * k3[i]);
		}
		f(t_n + a4 * tau, y_n4, N, k4);

		for (size_t i = 0; i < N; i++) {
			y[i] = y_n[i] + tau * (s1 * k1[i] + s2 * k2[i] + s3 * k3[i] + s4 * k4[i]);
		}

		t_n = t0 + (n+1) * tau;
		n++;

		for (size_t i = 0; i < N; i++) {
			y_n[i] = y[i];
		}
		//std::cout << n << "";
		for (size_t i = 0; i < N; i++) {
			fval[n][i] = y_n[i];
		}
	}


	delete[] y;
	delete[] y_n;
	delete[] y_n2; delete[] y_n3; delete[] y_n4;
}

void NDSolveAdams(T* (*f)(T, T*, size_t, T*), T* initialConditions, size_t N, T t0, T tmax, T tau0 = TAU0) {


	size_t k = 3;//size_t k = 4;

	T* y = new T[N]; //y_{n+1}
	T* y_n = new T[N];


	T tau = tau0;
	size_t n = 1;
	T t_n = t0;
	for (size_t i = 0; i < N; i++) {
		y_n[i] = initialConditions[i];//y_i(t_0);
	}

	std::ofstream out;          // поток для записи
	out.open("D:\\АДАМС.txt");

	T* mesh = new T[k + 1];//возможно заменить mesh на просто x
	for (size_t i = 0; i <= k; i++) {
		mesh[i] = t0 + i * tau;
	}

	T** yval = new T*[k+1]; //f[i][j]: i -- номер функции, j=0..k -- k+1 предыдущих значений ее
	for (size_t i = 0; i < k+1; i++) {
		 yval[i] = new T [N];
	}

	initRungeKutta(f, initialConditions, N, k, t0, tau, yval);
	for (size_t i = 0; i < k; i++) {
		out << std::setprecision(13) << t_n;
		for (size_t j = 0; j < N; j++) {
			out << " " << std::setprecision(13) << yval[i][j];
		}
		out << std::endl;
		t_n += tau;
	}

	for (size_t i = 0; i < N; i++) {
		y_n[i] = yval[k][i];
	}
	//находим первые 4 значения (одно -- начальное) c помощью рунге кутты лучше наверно сделать
	//отдельную функцию ну и можно даже на шаги тау и 2тау тоже отдельные функции

	T* k1 = new T[N];
	T* k2 = new T[N];
	T* k3 = new T[N];
	T* k4 = new T[N];
	//T* k5 = new T[N];
	while (t_n <= tmax + EPS) {


		out << std::setprecision(13) << t_n;
		for (size_t i = 0; i < N; i++) {
			out << " " << std::setprecision(13) << y_n[i];
		}
		out << std::endl;


		for (size_t i = 0; i < N; i++) {
			f(mesh[0], yval[k - 3], N, k1);
			f(mesh[1], yval[k - 2], N, k2);
			f(mesh[2], yval[k - 1], N, k3);
			f(mesh[3], yval[k], N, k4);
			//f(mesh[0], yval[k - 4], N, k1);
			//f(mesh[1], yval[k - 3], N, k2);
			//f(mesh[2], yval[k - 2], N, k3);
			//f(mesh[3], yval[k - 1], N, k4);
			//f(mesh[4], yval[k], N, k5);
			y[i] = yval[k][i] + tau * (55 * k4[i] - 59 * k3[i] +
				37 * k2[i] - 9 * k1[i]) / 24;
			//y[i] = yval[k][i] + tau * (1901 / 720 *k5[i]  - 1387 / 720 * k4[i] +
			//	109 / 30 * k3[i] - 637 / 360 * k2[i] + 251 / 720 *k1[i]);

		}

		t_n += tau;
		n++;

		for (size_t i = 0; i < N; i++) {
			y_n[i] = y[i];
		}

		for (size_t i = 0; i < k; i++) {
			mesh[i] = mesh[i + 1];//ну это мы сохраняем, то есть надо просто перемещать их
		}
		mesh[k] = mesh[k - 1] + tau;
		for (size_t i = 0; i < k; i++) {
			for (size_t j = 0; j < N; j++) {
				yval[i][j] = yval[i+1][j];//ну это мы сохраняем, то есть надо просто перемещать их
			}
		}
		for (size_t j = 0; j < N; j++) {
			yval[k][j] = y_n[j];
		}

	}

	delete[] yval;
	delete[] mesh;
	delete[] y;
	delete[] y_n;
	delete[] k1; delete[] k2; delete[] k3; delete[] k4;// delete[] k5;
}

void NDSolveCorrection(T* (*f)(T, T*, size_t, T*), T* initialConditions, size_t N, T t0, T tmax, T tau0 = TAU0) {

	size_t k = 3;

	T* y = new T[N]; //y_{n+1}
	T* y_n = new T[N];


	T tau = tau0;
	size_t n = 1;
	T t_n = t0;
	for (size_t i = 0; i < N; i++) {
		y_n[i] = initialConditions[i];//y_i(t_0);
	}

	std::ofstream out;          // поток для записи
	out.open("D:\\КОРЕКЦИЯ.txt");

	T* mesh = new T[k + 1];//возможно заменить mesh на просто x
	for (size_t i = 0; i <= k; i++) {
		mesh[i] = t0 + i * tau;
	}

	T** yval = new T * [k + 1]; //f[i][j]: i -- номер функции, j=0..k -- k+1 предыдущих значений ее
	for (size_t i = 0; i < k + 1; i++) {
		yval[i] = new T[N];
	}

	initRungeKutta(f, initialConditions, N, k, t0, tau, yval);
	for (size_t i = 0; i < k; i++) {
		out << std::setprecision(13) << t_n;
		for (size_t j = 0; j < N; j++) {
			out << " " << std::setprecision(13) << yval[i][j];
		}
		out << std::endl;
		t_n += tau;
	}

	for (size_t i = 0; i < N; i++) {
		y_n[i] = yval[k][i];
	}
	//находим первые 4 значения (одно -- начальное) c помощью рунге кутты лучше наверно сделать
	//отдельную функцию ну и можно даже на шаги тау и 2тау тоже отдельные функции

	T* k1 = new T[N];
	T* k2 = new T[N];
	T* k3 = new T[N];
	T* k4 = new T[N];
	while (t_n <= tmax + EPS) {


		out << std::setprecision(13) << t_n;
		for (size_t i = 0; i < N; i++) {
			out << " " << std::setprecision(13) << y_n[i];
		}
		out << std::endl;


		for (size_t i = 0; i < N; i++) {
			f(mesh[0], yval[k - 3], N, k1);
			f(mesh[1], yval[k - 2], N, k2);
			f(mesh[2], yval[k - 1], N, k3);
			f(mesh[3], yval[k], N, k4);
			y[i] = yval[k][i] + tau * (55 * k4[i] - 59 * k3[i] +
				37* k2[i] - 9 * k1[i]) / 24;

		}
		t_n +=tau;
		f(t_n, y, N, k1);

		for (size_t i = 0; i < N; i++) {
			y[i]=y_n[i]+ tau * (9 * k1[i] + 19 * k4[i] -
				5 * k3[i] + 1 * k2[i]) / 24;
		}
		n++;

		for (size_t i = 0; i < N; i++) {
			y_n[i] = y[i];
		}

		for (size_t i = 0; i < k; i++) {
			mesh[i] = mesh[i + 1];//ну это мы сохраняем, то есть надо просто перемещать их
		}
		mesh[k] = mesh[k - 1] + tau;
		for (size_t i = 0; i < k; i++) {
			for (size_t j = 0; j < N; j++) {
				yval[i][j] = yval[i + 1][j];//ну это мы сохраняем, то есть надо просто перемещать их
			}
		}
		for (size_t j = 0; j < N; j++) {
			yval[k][j] = y_n[j];
		}

	}

	delete[] yval;
	delete[] mesh;
	delete[] y;
	delete[] y_n;
	delete[] k1; delete[] k2; delete[] k3; delete[] k4; 
}

void NDSolve(T* (*f)(T, T*, size_t, T*), T* initialConditions, size_t N, T t0, T tmax, SOLVE_TYPE st, T tau0) {


	switch (st) {
	case ST_RUNGEKUTTA2:
		NDSolveRungeKutta2(f, initialConditions, N, t0, tmax, tau0);
		break;

	case ST_RUNGEKUTTA4:
		NDSolveRungeKutta4(f, initialConditions, N, t0, tmax, tau0);
		break;

	case ST_RUNGEKUTTA4_M:
		NDSolveRungeKutta4_M(f, initialConditions, N, t0, tmax, tau0);
		break;

	case ST_EXPLICIT_EULER:
		NDSolveExplicitEuler(f, initialConditions, N, t0, tmax, tau0);
		break;

	case ST_IMPLICIT_EULER:
		NDSolveImplicitEuler(f, initialConditions, N, t0, tmax, tau0);
		break;

	case ST_SYMMETRIC_SCHEME:
		NDSolveSymmetric(f, initialConditions, N, t0, tmax, tau0);
		break;
	case ST_ADAMS:
		NDSolveAdams(f, initialConditions, N, t0, tmax, tau0);
		break;
	case ST_CORRECTION:
		NDSolveCorrection(f, initialConditions, N, t0, tmax, tau0);
		break;
	}

}





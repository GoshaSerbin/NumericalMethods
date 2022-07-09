#include"header.h"
#include <iostream>
#include <fstream>
#include <sstream>

void makeTable1(T* (*f)(T, T*, size_t, T*), T* initialConditions, size_t N, T t0, T tmax, T* trueSolve, T tau0) {
	T* u_k = new T[N];//можно передавать
	std::cout << "ST_RUNGEKUTTA2\n";
	u_k[0] = trueSolve[0]; u_k[1] = trueSolve[1];//sin в 2пи * n
	//Рунге кутта 2-го порядка
	{


		std::cout << "         n   &    Шаг сетки tau_n    &  Норма ошибки || y_k-u_k ||   &      Отношение ошибок z_n     &       Порядок сходимости p_n    \n\n";
		T normLast;
		for (size_t n = 1; n <= 5; n++) {
			T tau = tau0 / pow(2, n - 1);
			std::cout << std::setw(12) << n << " & " << '\t';

			NDSolve(f, initialConditions, N, t0, tmax, ST_RUNGEKUTTA2, tau);
			T* y_k = new T[N];
			std::ifstream in("D:\\РУНГЕ-КУТТА2.txt"); // окрываем файл для чтения
			if (in.is_open())
			{
				T t;
				size_t count = floor((tmax - t0) / tau);
				std::string l1;
				for (size_t i = 0; i < count; i++) {
					std::getline(in, l1);                      // Read the current line
				}
				in >> t;
				//std::cout << "t=" << t << "  ";
				for (size_t i = 0; i < N; i++) {
					in >> y_k[i];
				}
				in.close();

			}
			//std::cout << "yyy" << y_k[0];

			T* dy = new T[N];
			for (size_t i = 0; i < N; i++) {
				dy[i] = y_k[i] - u_k[i];
				//	std::cout << "y=" << y_k[i]<< "u=" << u_k[i];
			}
			T normNew = normInf(dy, N);

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

			delete[] y_k;



		}
		std::cout << "----------------------------------------------------------------------------------------\n\n";


	}


	std::cout << "ST_RUNGEKUTTA4\n";
	//Рунге кутта 4-го порядка
	{
		std::cout << "         n   &    Шаг сетки tau_n    &  Норма ошибки || y_k-u_k ||   &      Отношение ошибок z_n     &       Порядок сходимости p_n    \n\n";
		T normLast;
		for (size_t n = 1; n <= 5; n++) {
			T tau = tau0 / pow(2, n - 1);
			std::cout << std::setw(12) << n << " & " << '\t';

			NDSolve(f, initialConditions, N, t0, tmax, ST_RUNGEKUTTA4, tau);
			T* y_k = new T[N];
			std::ifstream in("D:\\РУНГЕ-КУТТА4.txt"); // окрываем файл для чтения
			if (in.is_open())
			{
				T t;
				size_t count = floor((tmax - t0) / tau);
				std::string l1;
				for (size_t i = 0; i < count; i++) {
					std::getline(in, l1);                      // Read the current line
				}
				in >> t;
				//std::cout << "t=" << t << "  ";
				for (size_t i = 0; i < N; i++) {
					in >> y_k[i];
				}
				in.close();

			}
			//std::cout << "yyy" << y_k[0];

			T* dy = new T[N];
			for (size_t i = 0; i < N; i++) {
				dy[i] = y_k[i] - u_k[i];
				//	std::cout << "y=" << y_k[i]<< "u=" << u_k[i];
			}
			T normNew = normInf(dy, N);

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

			delete[] y_k;



		}
		std::cout << "----------------------------------------------------------------------------------------\n\n";


	}


	std::cout << "ST_EXPLICIT_EULER\n";
	//ST_EXPLICIT_EULER
	{
		std::cout << "         n   &    Шаг сетки tau_n    &  Норма ошибки || y_k-u_k ||   &      Отношение ошибок z_n     &       Порядок сходимости p_n    \n\n";
		T normLast;
		for (size_t n = 1; n <= 5; n++) {
			T tau = tau0 / pow(2, n - 1);
			std::cout << std::setw(12) << n << " & " << '\t';

			NDSolve(f, initialConditions, N, t0, tmax, ST_EXPLICIT_EULER, tau);
			T* y_k = new T[N];
			std::ifstream in("D:\\ЯВНЫЙ ЭЙЛЕР.txt"); // окрываем файл для чтения
			if (in.is_open())
			{
				T t;
				size_t count = floor((tmax - t0) / tau);
				std::string l1;
				for (size_t i = 0; i < count; i++) {
					std::getline(in, l1);                      // Read the current line
				}
				in >> t;
				//std::cout << "t=" << t << "  ";
				for (size_t i = 0; i < N; i++) {
					in >> y_k[i];
				}
				in.close();

			}
			//std::cout << "yyy" << y_k[0];

			T* dy = new T[N];
			for (size_t i = 0; i < N; i++) {
				dy[i] = y_k[i] - u_k[i];
				//	std::cout << "y=" << y_k[i]<< "u=" << u_k[i];
			}
			T normNew = norm1(dy, N);

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

			delete[] y_k;



		}
		std::cout << "----------------------------------------------------------------------------------------\n\n";


	}


	std::cout << "ST_IMPLICIT_EULER\n";
	//ST_IMPLICIT_EULER
	{
		std::cout << "         n   &    Шаг сетки tau_n    &  Норма ошибки || y_k-u_k ||   &      Отношение ошибок z_n     &       Порядок сходимости p_n    \n\n";
		T normLast;
		for (size_t n = 1; n <= 5; n++) {
			T tau = tau0 / pow(2, n - 1);
			std::cout << std::setw(12) << n << " & " << '\t';

			NDSolve(f, initialConditions, N, t0, tmax, ST_IMPLICIT_EULER, tau);
			T* y_k = new T[N];
			std::ifstream in("D:\\НЕЯВНЫЙ ЭЙЛЕР.txt"); // окрываем файл для чтения
			if (in.is_open())
			{
				T t;
				size_t count = floor((tmax - t0) / tau);
				std::string l1;
				for (size_t i = 0; i < count; i++) {
					std::getline(in, l1);                      // Read the current line
				}
				in >> t;
				//std::cout << "t=" << t << "  ";
				for (size_t i = 0; i < N; i++) {
					in >> y_k[i];
				}
				in.close();

			}
			//std::cout << "yyy" << y_k[0];

			T* dy = new T[N];
			for (size_t i = 0; i < N; i++) {
				dy[i] = y_k[i] - u_k[i];
				//	std::cout << "y=" << y_k[i]<< "u=" << u_k[i];
			}
			T normNew = normInf(dy, N);

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

			delete[] y_k;



		}
		std::cout << "----------------------------------------------------------------------------------------\n\n";


	}

	std::cout << "ST_SYMMETRIC_SCHEME\n";
	{
		std::cout << "         n   &    Шаг сетки tau_n    &  Норма ошибки || y_k-u_k ||   &      Отношение ошибок z_n     &       Порядок сходимости p_n    \n\n";
		T normLast;
		for (size_t n = 1; n <= 5; n++) {
			T tau = tau0 / pow(2, n - 1);
			std::cout << std::setw(12) << n << " & " << '\t';

			NDSolve(f, initialConditions, N, t0, tmax, ST_SYMMETRIC_SCHEME, tau);
			T* y_k = new T[N];
			std::ifstream in("D:\\СИММЕТРИЧНАЯ СХЕМА.txt"); // окрываем файл для чтения
			if (in.is_open())
			{
				T t;
				size_t count = floor((tmax - t0) / tau);
				std::string l1;
				for (size_t i = 0; i < count; i++) {
					std::getline(in, l1);                      // Read the current line
				}
				in >> t;
				//std::cout << "t=" << t << "  ";
				for (size_t i = 0; i < N; i++) {
					in >> y_k[i];
				}
				in.close();

			}
			//std::cout << "yyy" << y_k[0];

			T* dy = new T[N];
			for (size_t i = 0; i < N; i++) {
				dy[i] = y_k[i] - u_k[i];
				//	std::cout << "y=" << y_k[i]<< "u=" << u_k[i];
			}
			T normNew = normInf(dy, N);

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

			delete[] y_k;



		}
		std::cout << "----------------------------------------------------------------------------------------\n\n";

	}

	std::cout << "ST_ADAMS\n";
	{
		std::cout << "         n   &    Шаг сетки tau_n    &  Норма ошибки || y_k-u_k ||   &      Отношение ошибок z_n     &       Порядок сходимости p_n    \n\n";
		T normLast;
		for (size_t n = 1; n <= 5; n++) {
			T tau = tau0 / pow(2, n - 1);
			std::cout << std::setw(12) << n << " & " << '\t';

			NDSolve(f, initialConditions, N, t0, tmax, ST_ADAMS, tau);
			T* y_k = new T[N];
			std::ifstream in("D:\\АДАМС.txt"); // окрываем файл для чтения
			if (in.is_open())
			{
				T t;
				size_t count = floor((tmax - t0) / tau);
				std::string l1;
				for (size_t i = 0; i < count; i++) {
					std::getline(in, l1);                      // Read the current line
				}
				in >> t;
				//std::cout << "t=" << t << "  ";
				for (size_t i = 0; i < N; i++) {
					in >> y_k[i];
				}
				in.close();

			}
			//std::cout << "yyy" << y_k[0];

			T* dy = new T[N];
			for (size_t i = 0; i < N; i++) {
				dy[i] = y_k[i] - u_k[i];
				//	std::cout << "y=" << y_k[i]<< "u=" << u_k[i];
			}
			T normNew = normInf(dy, N);

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

			delete[] y_k;



		}
		std::cout << "----------------------------------------------------------------------------------------\n\n";

	}
	
	std::cout << "ST_CORRECTION\n";
	{

		std::cout << "         n   &    Шаг сетки tau_n    &  Норма ошибки || y_k-u_k ||   &      Отношение ошибок z_n     &       Порядок сходимости p_n    \n\n";
		T normLast;
		for (size_t n = 1; n <= 5; n++) {
			T tau = tau0 / pow(2, n - 1);
			std::cout << std::setw(12) << n << " & " << '\t';

			NDSolve(f, initialConditions, N, t0, tmax, ST_CORRECTION, tau);
			T* y_k = new T[N];
			std::ifstream in("D:\\КОРЕКЦИЯ.txt"); // окрываем файл для чтения
			if (in.is_open())
			{
				T t;
				size_t count = floor((tmax - t0) / tau);
				std::string l1;
				for (size_t i = 0; i < count; i++) {
					std::getline(in, l1);                      // Read the current line
				}
				in >> t;
				//std::cout << "t=" << t << "  ";
				for (size_t i = 0; i < N; i++) {
					in >> y_k[i];
				}
				in.close();

			}
			//std::cout << "yyy" << y_k[0];

			T* dy = new T[N];
			for (size_t i = 0; i < N; i++) {
				dy[i] = y_k[i] - u_k[i];
				//	std::cout << "y=" << y_k[i]<< "u=" << u_k[i];
			}
			T normNew = normInf(dy, N);

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

			delete[] y_k;



		}
		std::cout << "----------------------------------------------------------------------------------------\n\n";

	}
	
}

void makeTable2(T* (*f)(T, T*, size_t, T*), T* initialConditions, size_t N, T t0, T tmax, T tau0) {

	std::cout << "ST_RUNGEKUTTA2\n";
	//Рунге кутта 2-го порядка
	{
		T* y_k1 = new T[N];
		T* y_k2 = new T[N];
		T* y_k3 = new T[N];
		std::cout << "         n   &    Шаг сетки tau_n    &  Норма ошибки         &     Отношение w_n     & Порядок сходимости & оценка ошибки err \n\n";
		T normLast;
		for (size_t n = 1; n <= 5; n++) {
			T tau = tau0 / pow(2, n - 1);
			std::cout << std::setw(12) << n << " & " << '\t';
			NDSolve(f, initialConditions, N, t0, tmax, ST_RUNGEKUTTA2, tau);
			std::ifstream in("D:\\РУНГЕ-КУТТА2.txt"); // окрываем файл для чтения
			if (in.is_open())
			{
				T t;
				size_t count = floor((tmax - t0) / tau);
				std::string l1;
				for (size_t i = 0; i < count; i++) {
					std::getline(in, l1);                      // Read the current line
				}
				in >> t;
				//std::cout << "t=" << t << "  ";
				for (size_t i = 0; i < N; i++) {
					in >> y_k1[i];
				}
				in.close();
			}

			tau = tau0 / pow(2, n);
			NDSolve(f, initialConditions, N, t0, tmax, ST_RUNGEKUTTA2, tau);
			in.open("D:\\РУНГЕ-КУТТА2.txt");
			if (in.is_open())
			{
				T t;
				size_t count = floor((tmax - t0) / tau);
				std::string l1;
				for (size_t i = 0; i < count; i++) {
					std::getline(in, l1);                      // Read the current line
				}
				in >> t;
				//std::cout << "t=" << t << "  ";
				for (size_t i = 0; i < N; i++) {
					in >> y_k2[i];
				}
				in.close();

			}

			tau = tau0 / pow(2, n+1);
			NDSolve(f, initialConditions, N, t0, tmax, ST_RUNGEKUTTA2, tau);
			in.open("D:\\РУНГЕ-КУТТА2.txt");
			if (in.is_open())
			{
				T t;
				size_t count = floor((tmax - t0) / tau);
				std::string l1;
				for (size_t i = 0; i < count; i++) {
					std::getline(in, l1);                      // Read the current line
				}
				in >> t;
				//std::cout << "t=" << t << "  ";
				for (size_t i = 0; i < N; i++) {
					in >> y_k3[i];
				}
				in.close();

			}


			//std::cout << "yyy" << y_k1[0]<<" "<<y_k2[0]<<" "<< y_k3[0];

			T* dy1 = new T[N];
			for (size_t i = 0; i < N; i++) {
				dy1[i] = y_k2[i] - y_k1[i];
				//	std::cout << "y=" << y_k[i]<< "u=" << u_k[i];
			}

			T* dy2 = new T[N];
			for (size_t i = 0; i < N; i++) {
				dy2[i] = y_k3[i] - y_k2[i];
				//	std::cout << "y=" << y_k[i]<< "u=" << u_k[i];
			}

			T normNew1 = normInf(dy1, N);
			T normNew2 = normInf(dy2, N);

			delete[] dy1; delete[] dy2;
			T w;
			T pn;
			//считаем и выводим

			std::cout << std::setw(20) << tau << " & " << "\t";
			if (n <= 5) {
				w = normNew1 / normNew2;
				pn = log10(w) / log10(2);
				std::cout << std::setw(20) << normNew1 << " & " << "\t";
				std::cout << std::setw(20) << w << " & " << "\t";
				std::cout << std::setw(17) << pn << " & " << "\t";
				std::cout << std::setw(13) << normNew1/(1-pow(2,-pn)); std::cout << "\\\\\ \\hline\n";

			}






		}
		std::cout << "----------------------------------------------------------------------------------------\n\n";
		delete[] y_k1; delete[] y_k2; delete[] y_k3;
	}


}

void makePic1(T* (*f)(T, T*, size_t, T*), T* restPoint, T t0, T tau0) {
	T N = 2;
	T x1 = -2; T x2 = 2;
	T y1 = -2; T y2 = 2;
	T h = (x2 - x1) / 10.;

	T* fval1 = new T[N];
	T* fval2 = new T[N];
	T* X = new T[N];
	T** J = new T * [N];
	for (size_t i = 0; i < N; i++) {
		J[i] = new T[N];
	}
	X[0] = restPoint[0]; X[1] = restPoint[1];
	for (size_t i = 0; i < N; i++) {
		for (size_t j = 0; j < N; j++) {
			T* Xn = new T[N];
			for (size_t k = 0; k < N; k++) {
				Xn[k] = X[k];
				if (k == j) {
					Xn[k] += eps;
				}
			}
			J[i][j] = (f(t0, Xn, N, fval1)[i] - f(t0, X, N, fval2)[i]) / eps;//Вместо функции возвращающей вектор логичней бы было сделать вектор из функций
			delete[] Xn;
		}
	}
	std::cout << J[0][0] << " " << J[0][1] << " " << J[1][0] << " " << J[1][1];
	T* k1 = new T[N];
	T* k2 = new T[N];
	T* X_m = new T[N];
	T* XNew = new T[N];

	T b21 = 0.5;
	T a2 = 0.5;
	T sigma1 = 0;
	T sigma2 = 1;

	std::ofstream out;          // поток для записи
	out.open("D:\\РИСУНОК.txt"); // окрываем файл для записи
	for (T x = x1; x <= x2+EPS; x += h) {
		for (T y = y1; y <= y2+EPS; y += h) {
			X[0] = x; X[1] = y;
			{
				for (size_t i = 0; i < N; i++) {
					X_m[i] = X[i] + tau0 * b21 * (J[i][0] * X[0] + J[i][1] * X[1]);
				}

				for (size_t i = 0; i < N; i++) {
					XNew[i] = X[i] + tau0 * (sigma1 * (J[i][0] * X[0] + J[i][1] * X[1]) + sigma2 * (J[i][0] * X_m[0] + J[i][1] * X_m[1]));
				}
			}

			out << std::setprecision(13) << x << " " << std::setprecision(13) << y << " ";
			for (size_t i = 0; i < N; i++) {
				out << std::setprecision(13) << XNew[i] - X[i]<<" ";
			}
			out << std::endl;

			//в отрисовке далее можно умножить на какой-нибудь коээффциент
		}
	}
	delete[] XNew;
	delete[] k1; delete[] k2;
	delete[] X;
	for (size_t i = 0; i < N; i++) {
		delete[] J[i];
	}
	delete[] J;
	delete[] fval1; delete[] fval2;
}
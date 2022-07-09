#include"header.h"

int readSettings(const char* const fileName, size_t& n, char* const matrixName, char* const vectorName, char* const solveName) {
	std::ifstream ifile;
	ifile.open(fileName);
	if (!ifile.is_open()) {
		std::cerr << " Error : file with settings is not open !\n";
		return 1;
	}
	ifile >> n;
	std::string str;
	getline(ifile, str);
	ifile >> matrixName;
	getline(ifile, str);
	ifile >> vectorName;
	getline(ifile, str);
	ifile >> solveName;
	ifile.close();
	return 0;

}

T* testSystem(T t, T* y, size_t N, T* vec) {
	
	vec[0] = y[1];
	vec[1] = -y[0];
	return vec;
}

T* testf(T t, T* y, size_t N, T* fval) {

	fval[0] = y[1];
	fval[1] = 4*sin(t)-5*y[0];
	return fval;
}

T* testf2(T t, T* y, size_t N, T* fval) {

	fval[0] = 2 * y[0] +y[1] * y[1] - 1;
	fval[1] = 6 * y[0] - y[1] * y[1] + 1;
	return fval;
}

//T* testf28(T t, T* y, size_t N, T* fval) {
//
//	fval[0] = 2 * y[0] + y[1] * y[1] - 1;
//	fval[1] = 6 * y[0] - y[1] * y[1] + 1;
//	return fval;
//
//
//		double* y = new double[dim];
//
//		double temp;
//
//		double p1, p2;
//		p1 = pow((y[0] + mu) * (y[0] + mu) + u[2] * u[2] + u[4] * u[4], 3 / 2);
//		p2 = pow((u[0] - (1 - mu)) * (u[0] - (1 - mu)) + u[2] * u[2] + u[4] * u[4], 3 / 2);
//		y[0] = u[1];
//		y[1] = 2 * u[3] + u[0] - (1 - mu) * (u[0] + mu) / p1 - mu * (u[0] - (1 - mu)) / p2;
//		y[2] = u[3];
//		y[3] = -2 * u[1] + u[2] - (1 - mu) * u[2] / p1 - mu * u[2] / p2;
//		y[4] = u[5];
//		y[5] = -4 * PI * PI * (1 - mu) * u[4] / p1 - 4 * PI * PI * mu * u[4] / p2;
//
//}



T* testf3(T t, T* y, size_t N, T* fval) {

	fval[0] = 2 * y[0];
	fval[1] = y[1] * 1;
	return fval;
}

T* testf7(T t, T* y, size_t N, T* fval) {


	T k1 = 0.4 + 0.05 * sin(0.1 * t); T k2 = 0.1 + 0.01 * cos(0.01 * t);
	T a1 = 0.5; T a2 = 0.7;
	fval[0] = k1 * y[0] - a1 * y[0] * y[1]; 
	fval[1] = -k2 * y[1] + a2 * y[0] * y[1];
	return fval;
}//Н.У. 0.4 2.4    от 0 до 100

T* testf22(T t, T* y, size_t N, T* fval) {
	T delta = 0.02;
	T a = 2.5;
	T b = -1.;
	T F = 0.5;
	T nu = sqrt(a) + 0.15;
	fval[0] = y[1];
	fval[1] = F * cos(nu * t) - a * y[0] - b * pow(y[0], 3) - 2 * delta * y[1]; 
	return fval;
}// 0.1 0.1    от 0 до 25

int main() {

#ifdef CHECK_MEMORY
	_CrtMemState _ms;
	_CrtMemCheckpoint(&_ms);
#endif

	setlocale(LC_ALL, "RUS");
	std::cout << "\n                                            Лабораторная работа №1 \n";


	size_t N = 2;//количество уравнений в системе
	T* (*f)(T,T*,size_t, T*) = testf;
	T t0 = 0; T tmax = 100;
	T* trueSolve = new T[N]; 
	trueSolve[0] = sin(tmax); trueSolve[1] = cos(tmax);
	T* initCond = new T[N];
	initCond[0] = 0; initCond[1] = 1.0;
	std::cout << "\n ТЕСТ 1.1 \n";
	//NDSolve(f, initCond, N, t0, tmax, ST_RUNGEKUTTA2);
	std::cout << "\n ТЕСТ 1.2 \n";
	//NDSolve(f, initCond, N, t0, tmax, ST_EXPLICIT_EULER);
	std::cout << "\n ТЕСТ 1.3 \n";
	//NDSolve(f, initCond, N, t0, tmax, ST_IMPLICIT_EULER);
	std::cout << "\n ТЕСТ 1.4 \n";
	//NDSolve(f, initCond, N, t0, tmax, ST_SYMMETRIC_SCHEME);
	std::cout << "\n ТЕСТ 1.5 \n";
	//NDSolve(f, initCond, N, t0, tmax, ST_RUNGEKUTTA4_M);
	std::cout << "\n ТЕСТ 1.6 \n";
	//NDSolve(f, initCond, N, t0, tmax, ST_ADAMS);
	std::cout << "\n ТЕСТ 1.7 \n";
	//NDSolve(f, initCond, N, t0, tmax, ST_CORRECTION);
	std::cout << "                                                  Таблица 1\n";
	makeTable1(f, initCond, N, t0, tmax, trueSolve);
	std::cout << "\n\n                                                  Таблица 2\n";
	makeTable2(f, initCond, N, t0, tmax, 0.05);
	//makePic1(f, initCond, t0);

	delete[] initCond;
#ifdef CHECK_MEMORY
	_CrtMemDumpAllObjectsSince(&_ms);
#endif

	system("pause");
	return 0;
}
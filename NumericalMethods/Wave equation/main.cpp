#include"header.h"

T f1(T x) {
	//return sin(PI * x);
	return (x + 0.5) * (x + 0.5);
}

T ddf1(T x) {
	//return -PI*PI*sin(PI * x);
	return 2;
}

T f2(T x) {
	if (x < -1 / 3. + 2 || x >2 + 1 / 3.) return 0;
	else return 1;
}

T ddf2(T x) {
	return 0;
}

T g1(T x) {
	return 0;
	return (x + 1) * sin(x);
}

T phi1(T t) {
	return 0;
	return 0.5*(0.5 + t);
}

T psi1(T t) {
	return 0;
	return 2.25;
}

int main() {

#ifdef CHECK_MEMORY
	_CrtMemState _ms;
	_CrtMemCheckpoint(&_ms);
#endif

	setlocale(LC_ALL, "RUS");
	std::cout << "\n                                            Лабораторная работа №2 \n";
	T tmax = 4;			//max время
	T L = 4;				//длина струны
	T a2 = 1;				//скорость распространения малых возмущений

	T(*f)(T) = f2;	T(*ddf)(T) = ddf2;			//начальная форма струны
	T(*g)(T) = g1;			//начальная скорость струны
	T(*phi)(T) = phi1;	T(*psi)(T) = psi1;			//законы движения концов струны




	std::cout << "\n             ТЕСТ 1 \n";
	DSolveWaveEquation(L, tmax, a2, f, ddf, g, phi, psi, DT_DISCRETE);
	std::cout << "\n					Оценка порядка сходимости:\n";
	//MakeTable(L, tmax, a2, f, df, g, phi, psi);
	//std::cout << "                                                  Таблица 1\n";
	//makeTable1(f, initCond, N, t0, tmax, trueSolve);
	//std::cout << "\n\n                                                  Таблица 2\n";
	//makeTable2(f, initCond, N, t0, tmax, 0.05);
	//makePic1(f, initCond, t0);

#ifdef CHECK_MEMORY
	_CrtMemDumpAllObjectsSince(&_ms);
#endif

	system("pause");
	return 0;
}
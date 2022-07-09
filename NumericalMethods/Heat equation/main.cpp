#include"header.h"

T f(T x) {
	return 2 * x + 3 + 4 * sin(PI * x);
	//return 0;
	//return 0.1 + x*(1 - x);
	//return 0.1;
}

T boundF(T t) {
	//return 0.1;
	T u0 =10 ;
	return u0 * pow(t, 0.5);
}
T boundF1(T t) {
	return 3;

}
T boundF2(T t) {
	return 5;

}

T P1(T t,T t0,T Q) {
	if (t < t0) return Q;
	else return 0;
}

T P2(T t, T t0, T Q) {
	if (t < t0) return 2*Q*t;
	else return 0;
}

T P3(T t, T t0, T Q) {
	if (t < t0) return 2 * Q * (t0 - t);
	else return 0;
}

T P4(T t, T t0, T Q) {
	if (t < 0.5 * t0) return 2 * Q * t;
	else if (t >= t0) return 0;
	else return 2 * Q * (t0 - t);
}

//P=0
T PThermallyInsulted(T t, T t0, T Q) {
	return 0;
}

T Kx(T x) {
	return 7;
	T x1=2/5.; T x2=1/2.;
	T k1=1.5; T k2=0.1;
	if (x < x1) return k1;
	else if (x < x2) return k1*(x-x2)/(x1-x2)+k2* (x - x1) / (x2 - x1);
	else return k2;
}

T Ku(T u) {
	//T alpha = 0.5; T beta = 2; T gamma = 1.5;

	//return alpha + beta * pow(u, gamma);
	return 0.5 * u * u;
}




//испоользовать разные параметры для одного названия и решить проблему с возможными гран усл

//вообще можно оптимизировать и добавить отдельный случай дял сигма =1 или 0

int main() {

#ifdef CHECK_MEMORY
	_CrtMemState _ms;
	_CrtMemCheckpoint(&_ms);
#endif

	setlocale(LC_ALL, "RUS");
	std::cout << "\n                                            Лабораторная работа №2 \n";
	T tmax = 0.2;			//max время
	T L = 1;				//длина стержня
	T c = 1;				//удельная теплоемкость
	T ro = 1;				//плотность

	T(*K)(T) = Kx;			//коэфф-т теплопроводности

	T(*u0)(T) = f;			//начальные условия u(x,0)

	T t0=0.5; T Q=10;			//параметры теплового потока
	T (*P)(T, T, T) = P2;	//тепловой поток

	T (*ubound)(T) = boundF;			//температура на конце стержня


	
	std::cout << "\n             ТЕСТ 1 \n";
	//DSolveHeatEquation(L, tmax, c, ro, K, TC_X, u0, ubound, P, t0, Q);

	std::cout << "\n             ТЕСТ 2 \n";//там где ubound по сути не надо, тк достаточно u(x,0)
	K = Ku;
	//DSolveHeatEquation(L, tmax, c, ro, K, TC_U, u0, ubound, P, t0, Q);

	std::cout << "\n             ТЕСТ 3 \n";
	K = Kx;
	//DSolveHeatEquation(L, tmax, c, ro, K, TC_X, u0, P, t0, Q, ubound);

	//std::cout << "\n             ТЕСТ термоизолированные концы \n";
	//DSolveHeatEquation(L, tmax, c, ro, K, TC_X, u0, PThermallyInsulted, PThermallyInsulted, t0, Q);

	std::cout << "\n             ТЕСТ из методы №3 \n";
	K = Ku;

	//DSolveHeatEquation(L, tmax, c, ro, K, TC_U, u0, ubound, PThermallyInsulted, t0, Q);
	std::cout << "\n             ТЕСТ 5 ПОРЯДОК АППРОКСИМАЦИИ\n";
	K = Kx;
	T(*ubound1)(T) = boundF1; T(*ubound2)(T) = boundF2;
	table2(L, tmax, c, ro, K, TC_X, u0, ubound1, ubound2);
	//DSolveHeatEquation(L, tmax, c, ro, K, TC_X, u0, P, t0, Q, PThermallyInsulted);

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
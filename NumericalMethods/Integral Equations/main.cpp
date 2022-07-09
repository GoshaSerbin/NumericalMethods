#include"header.h"

T f1(T x) {
	return x+exp(-x);
}

T ftest(T x) {
	return (1+sin(x))*0.5;
}

T f2(T x) {
	return log(x+1)+x/2;
}

T ftest2(T x) {
	return x * x - x * x * x;
}

T fvec(T x, T y) {
	size_t N = 17;

	return sin((N + 1) / 2 * Arctg(x, y)) ;
}

T Ktest(T x, T s) {
	return (1 - x * cos(x * s));
}

T Ktest2(T x, T s) {
	return 3*x*x*x;
}

T K1(T x, T s) { //x2=0
	return x*(exp(-x*s)-1);
	//return 1 - x1;
}

T x1(T x) {
	return 1;
}

T x2(T x) {
	return -x;
}

T s3(T s) {
	return s*s/2.;
}

T x3(T x) {
	return pow(x, 3);
}

T x4(T x) {
	return -pow(x, 5);
}

T s4(T s) {
	return pow(s, 4) / 24.;
}

T x5(T x) {
	return pow(x, 7);
}

T s5(T s) {
	return pow(s, 6) / 720.;
}

T x6(T x) {
	return -pow(x, 9);
}

T s6(T s) {
	return pow(s, 8) / 40320.;
}

T x7(T x) {
	return pow(x, 11);
}

T s7(T s) {
	return pow(s, 10) / 3628800.;
}

int main() {

#ifdef CHECK_MEMORY
	_CrtMemState _ms;
	_CrtMemCheckpoint(&_ms);
#endif

	setlocale(LC_ALL, "RUS");

	//std::cout << Arctg(0,-1);
	std::cout << "\n                                            Ëàáîðàòîðíàÿ ðàáîòà ¹5 \n";
	T a = 0;
	T b = 1;
	T lambda = 0.5;


	std::cout << "\n             ÒÅÑÒ 1 \n";
	SolveFredgolm2(Ktest, ftest, a, b, lambda, ST_QUADRATURES);
	std::cout << "\n             ÒÅÑÒ 2 \n";
	//SolveFredgolm2(Ktest2, ftest2, a, b, lambda, ST_SIMPLE_ITERATION);
	std::cout << "\n             ÒÅÑÒ 3 \n";
	T(*psi[7])(T); T(*phi[7])(T);
	psi[0] = x1; psi[1] = x1; psi[2] = s3; psi[3] = s4; psi[4] = s5; psi[5] = s6; psi[6] = s7;
	phi[0] = x1; phi[1] = x2; phi[2] = x3; phi[3] = x4; phi[4] = x5; phi[5] = x6; phi[6] = x7;
	//SolveFredgolm2(7, psi, phi, ftest, a, b, lambda);
	std::cout << "\n             ÒÅÑÒ 4 \n";
	SolveSingularIntegralEquation(fvec);

#ifdef CHECK_MEMORY
	_CrtMemDumpAllObjectsSince(&_ms);
#endif

	system("pause");
	return 0;
}
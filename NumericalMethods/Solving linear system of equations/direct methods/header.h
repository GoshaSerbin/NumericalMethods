#include<fstream>
#include<iostream>
#include<string>
#include<math.h>
#include<iomanip>
#include<random>

typedef  float T;

#ifndef CONSTANTS_H
#define CONSTANTS_H
	const size_t countOfTimeLayers = 6;
	const size_t maxLengthOfFileName = 64;
	const T disturbancePercent = 0.01;
	const size_t numOfDisturbances = 10000;
	const T EPS = 1e-8;
#else
	extern const size_t countOfTimeLayers;
	extern const size_t maxLengthOfFileName;
	extern const double epsT;
#endif


#define CHECK_MEMORY
#ifdef CHECK_MEMORY
	#include<crtdbg.h>
	#define _CRTDBG_MAP_ALLOC
	#define new new (_NORMAL_BLOCK, __FILE__, __LINE__)
#endif


	void showLinearSystem(const size_t& n, T** const A, T* const b);

	T** inverse(size_t n, T** const A);

	void showMatrix(const size_t& n, T** const A);

	void showMatrixMultiplication(size_t n, T** const A, T** const B);

	T norm1(size_t n, T** A);

	T norm1(size_t n, T* x);

	T normInf(size_t n, T** A);

	T normInf(size_t n, T* x);

	void transpose(size_t n, T** const A);

	int gauss_method(const size_t& n, T** const A, T* const b, T* x);

	void gauss_reverse(const size_t& n, T** const A, T* const b, T* x);

	int QR_method(const size_t& n, T** const A, T* const b, T* x);

	void QR_decomposition(size_t n, T** A, T** Q, T** R);

	void QR_solve(size_t n, T** Q, T** R, T* b, T* x);

	void disturbance(size_t n, T* b);

	T residual(const size_t& n, T** const A, T* const b, T* x);





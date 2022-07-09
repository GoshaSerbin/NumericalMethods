#include<fstream>
#include<iostream>
#include<string>
#include<math.h>
#include<iomanip>
#include<random>

typedef  double T;

#ifndef CONSTANTS_H
#define CONSTANTS_H

enum STOP_CRITERION {
	CRITERION_NORM,
	CRITERION_TRUESOLVE,
	CRITERION_EPS1,
	CRITERION_RESIDUAL
};

const size_t maxLengthOfFileName = 64;
const size_t sizeOf3dMatrix = 200;
const size_t numOfDisturbances = 10000;
const T EPS = 1e-7;
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

void MatrixMultiplication(size_t n, T** const A, T** const B, T** C);

T norm1(size_t n, T** A);

T norm1(size_t n, T* x);

T normInf(size_t n, T** A);

T normInf(size_t n, T* x);

void transpose(size_t n, T** const A);

int simpleIterationMethod(const size_t& n, T** const A, T* const b, T* x, T tau, STOP_CRITERION sc, T* trueSol);

int JacobiMethod(const size_t& n, T** const A, T* const b, T* x, STOP_CRITERION sc, T* trueSol);

int relaxationMethod(const size_t& n, T** const A, T* const b, T* x, T omega, STOP_CRITERION sc, T* trueSol);

int SeidelsMethod(const size_t& n, T** const A, T* const b, T* x, STOP_CRITERION sc, T* trueSol);

//int relaxationMethod(const size_t& n, T** const A, T* const b, T* x, T omega = 1);

int relaxationMethod_3d(const size_t& n, T* const a, T* const b, T* const c, T* const d, T* x, T omega, STOP_CRITERION sc, T* trueSol);


T residual(const size_t& n, T** const A, T* const b, T* x);


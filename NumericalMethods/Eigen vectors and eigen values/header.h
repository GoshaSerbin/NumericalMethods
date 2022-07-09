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
const T EPS = 1e-6;
const T EPSZERO = 1e-15;
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

T normF(size_t n, T* x);

void toHessenbergMatrix(size_t n, T** const A, T** const H);

void transpose(size_t n, T** const A);

void QR_decomposition(size_t n, T** A, T** Q, T** R);

void QR_decForHessenburg(size_t n, T** A, T** Q, T** R);

void QR_solve(size_t n, T** Q, T** R, T* b, T* x);

void eigenvaluesSearch(size_t n, T** const A, bool withShift, bool withHessenburg, T* eigenvalues);

void eigenvectorsSearch(size_t n, T** const A, T* eigenvalues, T** eigenvectors);

void eigenSearch(size_t n, T** const A, T* eigenvalues, T** eigenvectors);

T residual(const size_t& n, T** const A, T* const b, T* x);



#pragma once
#include<fstream>
#include<iostream>
#include<string>
#include<math.h>
#include<iomanip>
#include<random>

typedef  double T;

#ifndef CONSTANTS_H
#define CONSTANTS_H

enum SOLVE_TYPE {
	ST_RUNGEKUTTA2,
	ST_RUNGEKUTTA4,
	ST_RUNGEKUTTA4_M,
	ST_EXPLICIT_EULER,
	ST_IMPLICIT_EULER,
	ST_SYMMETRIC_SCHEME,
	ST_ADAMS,
	ST_CORRECTION
};

enum DIFFERENTIAL_TYPE {
	DT_ANALYTIC,
	DT_BACK,
	DT_FORWARD
};

const size_t maxLengthOfFileName = 64;
const T difference = 1e-8;//1e-5 было для порядка сходимости
const T eps = 1e-7;
const T EPSZERO = 1e-5;
const T EPS = 1e-8;
const T TAU0 = 5e-2;
const T PI = 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651328230664709384;
#else
extern const size_t maxLengthOfFileName;
#endif


#define CHECK_MEMORY
#ifdef CHECK_MEMORY
#include<crtdbg.h>
#define _CRTDBG_MAP_ALLOC
#define new new (_NORMAL_BLOCK, __FILE__, __LINE__)
#endif

void transpose(size_t n, T** const A);

int QR_method(const size_t& n, T** const A, T* const b, T* x);

T* systemForImplicitEuler(T* (*f)(T, T*, size_t, T*), T t, T* y_n, T tau, T* y, size_t N, T* vec);

size_t NSolveNewton(T* (*system)(T* (*f)(T, T*, size_t, T*), T, T*, T, T*, size_t, T*), T* (*f)(T, T*, size_t, T*), T t, T* y_n, T tau, T* d1, T* d2, T* initY, size_t N, T* YNew);

T* systemForSymmetric(T* (*f)(T, T*, size_t, T*), T t_n, T* y_n, T tau, T* y, size_t N, T* vec);

void makeTable1(T* (*f)(T, T*, size_t, T*), T* initialConditions, size_t N, T t0, T tmax, T* trueSolve, T tau0 = TAU0);

void makeTable2(T* (*f)(T, T*, size_t, T*), T* initialConditions, size_t N, T t0, T tmax, T tau0 = TAU0);

void makePic1(T* (*f)(T, T*, size_t, T*), T* restPoint, T t0, T tau0 = 1e-5*TAU0);

void NDSolve(T* (*f)(T, T*, size_t, T*), T* initialConditions, size_t N, T t0, T tmax, SOLVE_TYPE st = ST_RUNGEKUTTA2, T tau0=TAU0);

T norm1(T* x, size_t n);

T normInf(T* x, size_t n);

T normInf(T** A, size_t n);

T norm1(T** A, size_t n);


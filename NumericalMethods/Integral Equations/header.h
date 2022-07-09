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
	ST_QUADRATURES,
	ST_SIMPLE_ITERATION,
};

const size_t maxLengthOfFileName = 64;
const T difference = 1e-8;//1e-5 было для порядка сходимости
const T eps = 1e-7;
const T EPSZERO = 1e-11;
const T EPS = 1e-5;
const T TAU0 = 5e-3;
const T PI = 3.14159265358979323846264338327950288419716939937510582;
#else
extern const size_t maxLengthOfFileName;
#endif


#define CHECK_MEMORY
#ifdef CHECK_MEMORY
#include<crtdbg.h>
#define _CRTDBG_MAP_ALLOC
#define new new (_NORMAL_BLOCK, __FILE__, __LINE__)
#endif


T Arctg(T x, T y);

int gauss_method(const size_t& n, T** const A, T* const b, T* x);

void SolveFredgolm2(T(*K)(T,T), T(*f)(T), T a, T b, T lambda, SOLVE_TYPE st = ST_QUADRATURES);

void SolveFredgolm2(size_t m, T(*psi[10])(T), T(*phi[10])(T), T(*f)(T), T a, T b, T lambda);

void SolveSingularIntegralEquation(T(*f)(T,T));

void makeTable1(T(*K)(T, T), T(*f)(T), T a, T b, T lambda);

T norm1(size_t n, T* x);

T normInf(size_t n, T* x);




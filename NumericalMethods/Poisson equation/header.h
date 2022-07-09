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

//вид граничных условий
enum BOUNDARY_CONDITIONS {
	BC_1,//1-го рода
	BC_2,//2-го рода
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
const T EPS = 1e-4;
const T TAU0 = 5e-3;
const T PI = 3.141592653589793238462643383279502884197169399375105820974944592307816406286;
#else
extern const size_t maxLengthOfFileName;
#endif


#define CHECK_MEMORY
#ifdef CHECK_MEMORY
#include<crtdbg.h>
#define _CRTDBG_MAP_ALLOC
#define new new (_NORMAL_BLOCK, __FILE__, __LINE__)
#endif


void DSolvePoissonEquation(T L1, T L2, T(*f)(T, T),
							BOUNDARY_CONDITIONS bc1, T (*G1)(T),
							BOUNDARY_CONDITIONS bc2, T(*G2)(T),
							BOUNDARY_CONDITIONS bc3, T(*G3)(T),
							BOUNDARY_CONDITIONS bc4, T(*G4)(T));
void showTable1(T L1, T L2, T(*f)(T, T),
	BOUNDARY_CONDITIONS bc1, T(*G1)(T),
	BOUNDARY_CONDITIONS bc2, T(*G2)(T),
	BOUNDARY_CONDITIONS bc3, T(*G3)(T),
	BOUNDARY_CONDITIONS bc4, T(*G4)(T));
T norm1(size_t n, T* x);

T normInf(size_t n, T* x);

T normInf(T** A, size_t n);

T norm1(T** A, size_t n);



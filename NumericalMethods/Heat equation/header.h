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

//коэффициент теплопроводности
enum THERMAL_CONDUCTIVITY_K {
	TC_X,//зависит от x
	TC_U//зависит от u
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
const T TAU0 = 5e-3;
const T PI = 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651328230664709384;
const T T0 = 1;
const T Q0 = 1;
#else
extern const size_t maxLengthOfFileName;
#endif


#define CHECK_MEMORY
#ifdef CHECK_MEMORY
#include<crtdbg.h>
#define _CRTDBG_MAP_ALLOC
#define new new (_NORMAL_BLOCK, __FILE__, __LINE__)
#endif

//в конце по умолчанию можно добавить шаги сетки
void DSolveHeatEquation(T L, T tmax, T с, T ro, T(*K)(T), THERMAL_CONDUCTIVITY_K tc, T(*u0)(T), T (*ubound)(T), T(*P)(T, T, T), T t0, T Q);

void DSolveHeatEquation(T L, T tmax, T с, T ro, T(*K)(T), THERMAL_CONDUCTIVITY_K tc, T(*u0)(T), T(*P)(T, T, T), T t0, T Q, T(*ubound)(T));

void DSolveHeatEquation(T L, T tmax, T с, T ro, T(*K)(T), THERMAL_CONDUCTIVITY_K tc, T(*u0)(T), T(*P0)(T, T, T), T(*P)(T, T, T), T t0, T Q);

void table(T L, T tmax, T с, T ro, T(*K)(T), THERMAL_CONDUCTIVITY_K tc, T(*u0)(T), T(*ubound1)(T), T(*ubound2)(T));

void table2(T L, T tmax, T с, T ro, T(*K)(T), THERMAL_CONDUCTIVITY_K tc, T(*u0)(T), T(*ubound1)(T), T(*ubound2)(T));


T norm1(size_t n, T* x);

T normInf(size_t n, T* x);

T normInf(T** A, size_t n);

T norm1(T** A, size_t n);


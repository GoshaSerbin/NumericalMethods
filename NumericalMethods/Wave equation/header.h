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

enum DIFFERENTIAL_TYPE {
	DT_ANALYTIC,
	DT_DISCRETE,
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

void DSolveWaveEquation(T L, T tmax, T a2,  T(*f)(T), T(*df)(T), T(*g)(T), T(*phi)(T), T(*psi)(T), DIFFERENTIAL_TYPE dt=DT_ANALYTIC);


void MakeTable(T L, T tmax, T a2, T(*f)(T), T(*ddf)(T), T(*g)(T), T(*phi)(T), T(*psi)(T));

T norm1(size_t n, T* x);

T normInf(size_t n, T* x);




#include<fstream>
#include<iostream>
#include<string>
#include<math.h>
#include<iomanip>
#include<random>

typedef  double T;

#ifndef CONSTANTS_H
#define CONSTANTS_H

enum MESH_TYPE {
	MT_UNIFORM,
	MT_CHEBYSHEV
};

const size_t maxLengthOfFileName = 64;
const size_t numOfPointsForGraphics = 512;//количество значений интерполянта для построения графика;;
const T PI = 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651328230664709384460955058223172535940812848111745028410270193852110555964462294895493038196442881097566593344612847564823378678316527120190914564856692;
const T EPSZERO = 1e-15;
#else
extern const size_t maxLengthOfFileName;
#endif


#define CHECK_MEMORY
#ifdef CHECK_MEMORY
#include<crtdbg.h>
#define _CRTDBG_MAP_ALLOC
#define new new (_NORMAL_BLOCK, __FILE__, __LINE__)
#endif

void makeMeshGrid(size_t n, T a, T b, MESH_TYPE mt, T* mesh);

size_t makeMeshGrid(T h, T a, T b, T* mesh);

void solve_3d(size_t n, T* a, T* b, T* c, T* d, T* x);

void LagrangeInterpolation(T* mesh, T* func, T A, T B, size_t num, T* interpolant, size_t m);

void splineInterpolation(T* mesh, T* func,T A, T B, size_t num, T* interpolant, size_t m);

void showVector(size_t num, T* mesh);

int output(size_t n, T* x, T* y, char const* const fileName);

void showTable1(T(*func)(T), T a, T b);

void showTable2(T(*func)(T), T a, T b, size_t N, size_t q, void (*method)(T*, T*,T,T, size_t, T*, size_t), MESH_TYPE mt = MT_UNIFORM);

void makeGraphic(T(*func)(T), T a, T b, size_t N, void(*method)(T*, T*,T,T, size_t, T*, size_t), MESH_TYPE mt = MT_UNIFORM,
	const char* const outFileName = "D:\\out_test.txt");
//Норма ошибки интерполяции функции func на отрезке [a,b] с использованием num точек для проверки
T errorNorm(T(*func)(T), T a, T b, T* interpolant, size_t num);

T norm1(size_t n, T** A);

T norm1(size_t n, T* x);

T normInf(size_t n, T** A);

T normInf(size_t n, T* x);

T normF(size_t n, T* x);


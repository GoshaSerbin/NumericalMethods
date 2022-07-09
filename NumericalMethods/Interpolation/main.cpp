#include"header.h"

int readSettings(const char* const fileName, size_t& n, char* const matrixName, char* const vectorName, char* const solveName) {
	std::ifstream ifile;
	ifile.open(fileName);
	if (!ifile.is_open()) {
		std::cerr << " Error : file with settings is not open !\n";
		return 1;
	}
	ifile >> n;
	std::string str;
	getline(ifile, str);
	ifile >> matrixName;
	getline(ifile, str);
	ifile >> vectorName;
	getline(ifile, str);
	ifile >> solveName;
	ifile.close();
	return 0;

}

T f1(T x) {
	return x*x;
}


T fConst(T x) {
	return 1;
}

T fLinear(T x) {
	return x;
}

T fTest17(T x) {
	

	return pow(4 * pow(x, 3) + 2 * pow(x, 2) - 4 * x + 2, sqrt(2)) +asin(-1 / (pow(x, 2) + x + 5)) - 5;
}

T fTest22(T x) {
	T R1 = log((-sqrt(3) * pow(x, 4) - pow(x, 2) + 5 * x + 1)) / log(10);
	T R2 = tanh((-pow(x, 5) - 2 * pow(x, 4) - pow(x, 3) +
	3 * pow(x, 2) + 6 * x + 3 - sqrt(5)) / (pow(x, 2) + 2 * x + 1));
	return R1 + R2 - 1.0;

}

T fRunge(T x){

	return 1 / (1 + 25 * x * x);
}

T fSinpx(T x) {
	return sin(PI * x);
}


int main() {

#ifdef CHECK_MEMORY
	_CrtMemState _ms;
	_CrtMemCheckpoint(&_ms);
#endif

	setlocale(LC_ALL, "RUS");


	MESH_TYPE mt;

	//for (int i = 0; i < numberOfTable; i++) {
	//	std::cout << CosT[i][1] << "\n";
	//}

	std::cout << "\n                                     ЗАДАНИЕ 1\n";
	T a = -1; T b = 1;
	T(*func)(T) = fTest17;
	std::cout << "ФУНКЦИЯ:  f(x) = из теста 17\n";
	std::cout << "ОТРЕЗОК:  [ " << a << "; " << b << "]\n\n";
	showTable1(func, a, b);

	std::cout << "                                       ЗАДАНИЕ 2 (a)\n";
	a = -1; b = 1;
	func = fConst;
	std::cout << "ФУНКЦИЯ:  f(x) = const = 1\n";
	std::cout << "ОТРЕЗОК:  [ " << a << "; " << b << "]\n\n";
	showTable1(func, a, b);

	std::cout << "                                       ЗАДАНИЕ 2 (b)\n";
	a = -1; b = 1;
	func = fLinear;
	std::cout << "ФУНКЦИЯ:  f(x) = x\n";
	std::cout << "ОТРЕЗОК:  [ " << a << "; " << b << "]\n\n";
	showTable1(func, a, b);

	std::cout << "                                       ЗАДАНИЕ 3 (функция Рунге) --- в Wolfram Mathematica\n\n";
	a = -1; b = 1;
	func = fRunge;

	for (size_t num = 4; num <= 32; num *= 2) {

		T* mesh = new T[num + 1];//возможно заменить mesh на просто x
		T* interpolant = new T[numOfPointsForGraphics + 1];

		mt = MT_UNIFORM;
		makeMeshGrid(num, a, b, mt, mesh);
		T* fval = new T[num + 1];

		for (size_t i = 0; i <= num; i++) {
			fval[i] = func(mesh[i]);

		}

		LagrangeInterpolation(mesh, fval,a,b, num, interpolant, numOfPointsForGraphics);

		T* x = new T[numOfPointsForGraphics + 1];
		makeMeshGrid(numOfPointsForGraphics, a, b, MT_UNIFORM, x);

		switch (num) {
		case 4:
			output(numOfPointsForGraphics, x, interpolant, "D:\\out_task3_uniform_num4.txt");
			break;
		case 8:
			output(numOfPointsForGraphics, x, interpolant, "D:\\out_task3_uniform_num8.txt");
			break;
		case 16:
			output(numOfPointsForGraphics, x, interpolant, "D:\\out_task3_uniform_num16.txt");
			break;
		case 32:
			output(numOfPointsForGraphics, x, interpolant, "D:\\out_task3_uniform_num32.txt");
			break;
		}


		mt = MT_CHEBYSHEV;
		makeMeshGrid(num, a, b, mt, mesh);
		for (size_t i = 0; i <= num; i++) {
			fval[i] = func(mesh[i]);
		}
		LagrangeInterpolation(mesh, fval,a,b, num, interpolant, numOfPointsForGraphics);
		makeMeshGrid(numOfPointsForGraphics, a, b, MT_UNIFORM, x);
		switch (num) {
		case 4:
			output(numOfPointsForGraphics, x, interpolant, "D:\\out_task3_chebyshev_num4.txt");
			break;
		case 8:
			output(numOfPointsForGraphics, x, interpolant, "D:\\out_task3_chebyshev_num8.txt");
			break;
		case 16:
			output(numOfPointsForGraphics, x, interpolant, "D:\\out_task3_chebyshev_num16.txt");
			break;
		case 32:
			output(numOfPointsForGraphics, x, interpolant, "D:\\out_task3_chebyshev_num32.txt");
			break;
		}
		delete[] interpolant;
		delete[] fval;
		delete[] mesh;
		delete[] x;

	}

	std::cout << "                                       ЗАДАНИЕ 4 (Сплайн интерполяция)\n";
	a = -1; b = 1;
	func = fLinear;
	size_t num = 128;
	T* mesh = new T[num + 1];//возможно заменить mesh на просто x
	T* interpolant = new T[numOfPointsForGraphics + 1];
	T* fval = new T[num + 1];
	T* x = new T[numOfPointsForGraphics + 1];
	{




		mt = MT_UNIFORM;
		makeMeshGrid(num, a, b, mt, mesh);


		for (size_t i = 0; i <= num; i++) {
			fval[i] = func(mesh[i]);

		}
		splineInterpolation(mesh, fval,a,b, num, interpolant, numOfPointsForGraphics);

		makeMeshGrid(numOfPointsForGraphics, a, b, MT_UNIFORM, x);
		std::cout << "Норма ошибки f(x)=x:  " << errorNorm(func, a, b, interpolant, numOfPointsForGraphics) << "\n";

		output(numOfPointsForGraphics, x, interpolant, "D:\\out_task4_1.txt");
	}

	func = f1;
	a = -1; b = 1;
	{
		makeMeshGrid(num, a, b, mt, mesh);

		for (size_t i = 0; i <= num; i++) {
			fval[i] = func(mesh[i]);

		}
		splineInterpolation(mesh, fval,a,b, num, interpolant, numOfPointsForGraphics);
		makeMeshGrid(numOfPointsForGraphics, a, b, MT_UNIFORM, x);
		std::cout << "Норма ошибки f(x) = x^2:  " << errorNorm(func, a, b, interpolant, numOfPointsForGraphics) << "\n";

		output(numOfPointsForGraphics, x, interpolant, "D:\\out_task4_2.txt");
	}


	//func = fTest17;
	//mt = MT_CHEBYSHEV;
	//makeMeshGrid(num, a, b, mt, mesh);


	//for (size_t i = 0; i <= num; i++) {
	//	fval[i] = func(mesh[i]);

	//}
	//splineInterpolation(mesh, fval, num, interpolant, numOfPointsForGraphics);

	//makeMeshGrid(numOfPointsForGraphics, a, b, MT_UNIFORM, x);

	//output(numOfPointsForGraphics, x, interpolant, "D:\\out_task0.txt");




	

	

	std::cout << "\n                                       ЗАДАНИЕ 5 (Скорость сходимости)\n";
	func = fSinpx;
	size_t N = 2;
	size_t q = 2;
	std::cout << "Начальное количество разбиений: " << N << "\n              q=" << q << "\n\n";
	a = -1; b = 1;
	
	std::cout << "\n\nЛагранж на равномерной сетке:\n";
	showTable2(func, a, b, N, q, LagrangeInterpolation);
	std::cout << "\n\nЛагранж на чебышевской сетке:\n";
	showTable2(func, a, b, N, q, LagrangeInterpolation, MT_CHEBYSHEV);
	std::cout << "\n\nСплайн-интерполяция на отрезке [-1;1]:\n";
	showTable2(func, a, b, N, q, splineInterpolation);
	std::cout << "\n\nСплайн-интерполяция на отрезке [-1.25;1.25]:\n";
	showTable2(func, -1.25, 1.25, N, q, splineInterpolation);




	makeGraphic(fTest17, -1, 1, 64, splineInterpolation);








	delete[] interpolant;
	delete[] fval;
	delete[] mesh;
	delete[] x;
	//T* a1 = new T[10];
	//T* b1 = new T[10];
	//T* c1 = new T[10];
	//T* d1 = new T[10];
	//T* x1 = new T[10];
	//for (size_t i = 0; i < 10; i++) {
	//	a1[i] = 0;
	//	b1[i] = 5;
	//	c1[i] = 0;
	//	d1[i] = 1;
	//}
	//solve_3d(10, a1, b1, c1, d1, x1);
	//std::cout << "\n";
	//for (size_t i = 0; i < 10; i++) {
	//	std::cout << x1[i] << "\t";
	//}




	//T* _a = new T[10];
	//T* _c = new T[10];
	//T* _b = new T[10];
	//T* _d = new T[10];
	//T* _x = new T[10];
	//T* trueSol3d = new T[10];
	//for (size_t i = 0; i < 10; i++) {
	//	trueSol3d[i] = 2 - (i + 1) % 2;
	//}
	//for (size_t i = 0; i < 10; i++) {
	//	_a[i] = 1;
	//	_b[i] = 4;
	//	_c[i] = 1;
	//	_d[i] = 10 - 2 * ((i + 1) % 2);
	//}
	//_a[0] = 0; _c[10 - 1] = 0;
	//_d[0] = 6; _d[10 - 1] = 9 - 3 * (10 % 2);

	//std::cout << "\n\n                          Трехдиагональная матрица:\n";
	//solve_3d(10, _a, _b, _c, _d, x1);



	//std::cout << "\nРешение системы для трехдиагональной матрицы:\n";
	//for (size_t i = 0; i < 10; i++) {
	//	std::cout << std::setw(15) << x1[i] << '\t';
	//}

#ifdef CHECK_MEMORY
	_CrtMemDumpAllObjectsSince(&_ms);
#endif

	system("pause");
	return 0;
}
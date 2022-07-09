#include"header.h"

//func = y; mesh = x;
void LagrangeInterpolation(T* mesh, T* func, T A, T B, size_t num, T* interpolant, size_t m) {

	//T h = (mesh[num] - mesh[0]) / (1.*m);
	T h = (B - A) / (1. * m);
	for (size_t i = 0; i <= m; i++) {

		T x = A + i * h;
		T Ln = 0;
		T* c = new T[num+1];
		for (size_t k = 0; k <= num; k++) {
			c[k] = 1;
			for (size_t j = 0; j <= num; j++) {
				if (j != k) {
					c[k] *= (x - mesh[j]) / (mesh[k] - mesh[j]);
				}
			}

		}
		for (size_t k = 0; k <= num; k++) {
			Ln += c[k] * func[k];
		}

		interpolant[i] = Ln;


		delete[] c;


	}
	//if ( num == 128 && interpolant[0] = 0) {
	//	std::cout << "FUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUCK";
	//}
	//else {
	//	std::cout << ")))))))))))))";
	//}
}

void solve_3d(size_t n, T* a, T* b, T* c, T* d, T* x) {

	T* alpha = new T[n+1];
	T* beta = new T[n+1];
	alpha[1] = -c[0]*1.0 / b[0];
	beta[1] = d[0]*1.0 / b[0];
	for (size_t i = 1; i < n - 1; i++) {
		T z = 1.*(b[i] + a[i] * alpha[i]);
		alpha[i + 1] = -c[i] / z;
		beta[i + 1] = (d[i] - a[i] * beta[i]) / z;
	}

	x[n - 1] = (d[n-1] - a[n-1] * beta[n-1]) / (b[n-1] + a[n-1] * alpha[n-1]);
	for(size_t i = n - 2; i > 0; i--) {
		x[i] = alpha[i + 1] * x[i + 1] + beta[i + 1];
	}
	x[0] = alpha[1] * x[1] + beta[1];
	delete[] alpha;
	delete[] beta;
	//std::cout << "коснтанты:\n";
	//for (size_t i = 0; i < n; i++) {
	//	std::cout << a[i] << "\t";
	//}
	//std::cout << "\n";
}



T errorNorm(T(*func)(T ), T a, T b, T* interpolant, size_t num) {

	T h = (b - a) / (1.*num);
	T* fval = new T[num+1];
	for (size_t i = 0; i <= num; i++) {
		fval[i] = func(a+i*h);
	}
	//std::cout << "shit" << func(a);

	T maxError = 0;
//	size_t ind;
	for (size_t i = 0; i <= num; i++) {
		T er = abs(fval[i] - interpolant[i]);
		if (er > maxError) {
			maxError = er; //ind = i;
		}
	}
	//std::cout << "dsdsds" << fval[ind];
	//std::cout << "ssdsds" << interpolant[num];
	delete[] fval;

	return maxError;
}

void showTable1(T(*func)(T), T a, T b) {

	MESH_TYPE mt;

	std::cout << "  Количество узлов   ||   Норма ошибки интерполяции   ||   Норма ошибки интерполяции    \n";
	std::cout << "                     ||      на равномерной сетке     ||     на чебышевской сетке    \n\n";
	for (size_t num = 2; num <= 128; num *= 2) {

		std::cout << std::setw(12) << num << " & " << '\t';

		T* mesh = new T[num + 1];//возможно заменить mesh на просто x
		T* interpolant = new T[numOfPointsForGraphics + 1];

		mt = MT_UNIFORM;
		makeMeshGrid(num, a, b, mt, mesh);
		T* fval = new T[num + 1];
		//	std::cout << "\n";
		for (size_t i = 0; i <= num; i++) {
			fval[i] = func(mesh[i]);
			//		std::cout << mesh[i] << " ";
		}
		//	std::cout << "\n";

		LagrangeInterpolation(mesh, fval, a, b, num, interpolant, numOfPointsForGraphics);

		T* x = new T[numOfPointsForGraphics + 1];
		makeMeshGrid(numOfPointsForGraphics, a, b, MT_UNIFORM, x);
		if (num == 4) {
			output(numOfPointsForGraphics, x, interpolant, "D:\\out_task1_uniform_num4.txt");
		}
		else {
			if (num == 128) {
				output(numOfPointsForGraphics, x, interpolant, "D:\\out_task1_uniform_num128.txt");
			}
		}

		std::cout << std::setw(23) << errorNorm(func, a, b, interpolant, numOfPointsForGraphics) << " & " << "\t\t";

		//for (size_t i = 0; i <= numOfPointsForGraphics; i++) {
		//	std::cout << "[" << x[i] << ";" << interpolant[i]<<"]\n";
		//}


		mt = MT_CHEBYSHEV;
		makeMeshGrid(num, a, b, mt, mesh);
		//	std::cout << "\n";
		for (size_t i = 0; i <= num; i++) {
			//	std::cout << mesh[i] << " ";
			fval[i] = func(mesh[i]);
		}
		//	std::cout << "\n";
			//numOfPointsForGraphics звучит не оч тк мы юзаем его не только для графиков
		LagrangeInterpolation(mesh, fval, a, b, num, interpolant, numOfPointsForGraphics);
		std::cout << std::setw(23) << errorNorm(func, a, b, interpolant, numOfPointsForGraphics) << "\\\\\ \\hline\n";
		makeMeshGrid(numOfPointsForGraphics, a, b, MT_UNIFORM, x);
		if (num == 4) {
			output(numOfPointsForGraphics, x, interpolant, "D:\\out_task1_chebyshev_num4.txt");
		}
		else {
			if (num == 128) {
				output(numOfPointsForGraphics, x, interpolant, "D:\\out_task1_chebyshev_num128.txt");
			}
		}
		
		delete[] interpolant;
		delete[] fval;
		delete[] mesh;
		delete[] x;

	}
	std::cout << "----------------------------------------------------------------------------------------\n\n";

}


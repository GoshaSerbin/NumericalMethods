#include"Header.h"

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

//Ввод данных из файлов
int input(const size_t& n, T** A, T* b, T* trueSol, char* const matrixName, char* const vectorName, char* const solveName) {

	std::ifstream ifile;
	ifile.open(matrixName);
	if (!ifile.is_open()) {
		std::cerr << " Error : file with matrix is not open !\n";
		return 1;
	}
	ifile >> A[0][0];
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < n; j++) {
			ifile >> A[i][j];
		}
	}
	ifile.close();

	ifile.open(vectorName);
	if (!ifile.is_open()) {
		std::cerr << " Error : file with vector is not open !\n";
		return 1;
	}
	ifile >> b[0];
	for (size_t i = 0; i < n; i++) {
		ifile >> b[i];
	}
	ifile.close();

	ifile.open(solveName);
	if (!ifile.is_open()) {
		std::cerr << " Error : file with vector is not open !\n";
		return 1;
	}
	ifile >> trueSol[0];
	for (size_t i = 0; i < n; i++) {
		ifile >> trueSol[i];
	}
	ifile.close();
	return 0;
}

int main() {

#ifdef CHECK_MEMORY
	_CrtMemState _ms;
	_CrtMemCheckpoint(&_ms);
#endif

	setlocale(LC_ALL, "RUS");

	size_t n;
	char matrixName[maxLengthOfFileName];
	char vectorName[maxLengthOfFileName];
	char solveName[maxLengthOfFileName];

	std::cout.precision(4);
	std::cout << "Считываем данные...\n";
	readSettings("D:\\settingsFile.txt", n, matrixName, vectorName, solveName);

	T** A = new T * [n];
	for (size_t i = 0; i < n; i++)	A[i] = new T[n];

	T* b = new T[n];

	T* trueSol = new T[n];

	STOP_CRITERION  sc = CRITERION_NORM;

	input(n, A, b, trueSol, matrixName, vectorName, solveName);

	std::cout << "\nИсходная СЛАУ:\n\n";
	showLinearSystem(n, A, b);

	T* x = new T[n];
	T* dx = new T[n];
	int er= simpleIterationMethod(n, A, b, x, 0.01, sc, trueSol);
	if (er == -1) {
		std::cout << "\nМетод расходится, т.к.  ||C||>1\n";
	}
	else {
		std::cout << "\nРешение системы:\n";
		for (size_t i = 0; i < n; i++) {
			std::cout << std::setw(15) << x[i] << '\t';
		}

		std::cout << "\nНевязка решения: " << residual(n, A, b, x);

		for (size_t i = 0; i < n; i++) { dx[i] = x[i] - trueSol[i]; }
		std::cout << "\nНорма вектора ошибки: " << norm1(n, dx);

	}




	er = JacobiMethod(n, A, b, x, sc, trueSol);
	if (er == -1) {
		std::cout << "\nМетод расходится, т.к.  ||C||>1\n";
	}
	else {
		std::cout << "\nРешение системы:\n";
		for (size_t i = 0; i < n; i++) {
			std::cout << std::setw(15) << x[i] << '\t';
		}
		std::cout << "\nНевязка решения: " << residual(n, A, b, x);

		for (size_t i = 0; i < n; i++) { dx[i] = x[i] - trueSol[i]; }
		std::cout << "\nНорма вектора ошибки: " << norm1(n, dx);
	}



	er = SeidelsMethod(n, A, b, x, sc, trueSol);
	if (er == -1) {
		std::cout << "\nМетод расходится, т.к.  ||C||>1\n";
	}
	else {
		std::cout << "\nРешение системы:\n";
		for (size_t i = 0; i < n; i++) {
			std::cout << std::setw(15) << x[i] << '\t';
		}
		std::cout << "\nНевязка решения: " << residual(n, A, b, x);

		for (size_t i = 0; i < n; i++) { dx[i] = x[i] - trueSol[i]; }
		std::cout << "\nНорма вектора ошибки: " << norm1(n, dx);
	}


	
	T* _a = new T[sizeOf3dMatrix];
	T* _c = new T[sizeOf3dMatrix];
	T* _b = new T[sizeOf3dMatrix];
	T* _d = new T[sizeOf3dMatrix];
	T* _x = new T[sizeOf3dMatrix];
	T* trueSol3d = new T[sizeOf3dMatrix];
	for (size_t i = 0; i < sizeOf3dMatrix; i++) {
		trueSol3d[i] = 2 - (i + 1) % 2;
	}
	for (size_t i = 0; i < sizeOf3dMatrix; i++) {
		_a[i] =1;
		_b[i] =4;
		_c[i] =1;
		_d[i] =10 - 2*((i+1)%2);
	}
	_a[0] = 0; _c[sizeOf3dMatrix - 1] = 0;
	_d[0] = 6; _d[sizeOf3dMatrix - 1] = 9 - 3 * (sizeOf3dMatrix % 2);
	
	std::cout << "\n\n                          Трехдиагональная матрица:\n";
	relaxationMethod_3d(sizeOf3dMatrix, _a, _b, _c, _d, _x, 1, sc, trueSol3d);



	//std::cout << "\nРешение системы для трехдиагональной матрицы:\n";
	//for (size_t i = 0; i < sizeOf3dMatrix; i++) {
	//	std::cout << std::setw(15) << _x[i] << '\t';
	//}
	std::cout << "\n Для матрицы 6 8 1 :\n";
	for (size_t i = 0; i < sizeOf3dMatrix; i++) {
		_a[i] = 6;
		_b[i] = 8;
		_c[i] = 1;
		_d[i] = i;
	}
	_a[0] = 0; _c[sizeOf3dMatrix - 1] = 0;
	relaxationMethod_3d(sizeOf3dMatrix, _a, _b, _c, _d, _x, 1, sc, trueSol3d);

	std::cout << "\n Для матрицы 1 8 6 :\n";
	for (size_t i = 0; i < sizeOf3dMatrix; i++) {
		_a[i] = 1;
		_b[i] = 8;
		_c[i] = 6;
		_d[i] = i;
	}
	_a[0] = 0; _c[sizeOf3dMatrix - 1] = 0;
	relaxationMethod_3d(sizeOf3dMatrix, _a, _b, _c, _d, _x, 1, sc, trueSol3d);

	delete[] _a;
	delete[] _b;
	delete[] _c;
	delete[] _d;
	delete[] _x;
	delete[] trueSol3d;
	
	for (size_t i = 0; i < n; i++) {
		delete[] A[i];
	}
	delete[] A;

	delete[] dx;
	delete[] b;
	delete[] x;
	delete[] trueSol;

#ifdef CHECK_MEMORY
	_CrtMemDumpAllObjectsSince(&_ms);
#endif

	system("pause");
	return 0;
}
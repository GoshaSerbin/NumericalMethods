#include"header.h"

int readSettings(const char* const fileName, size_t& n, char* const matrixName, char* const vectorName, char* const outputName) {
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
	ifile >> outputName;
	ifile.close();
	return 0;

}

//Ввод данных из файлов
int input(const size_t& n, T** A, T* b, char* const matrixName, char* const vectorName) {

	std::ifstream ifile;
	ifile.open(matrixName);
	ifile >> A[0][0];
	if (!ifile.is_open()) {
		std::cerr << " Error : file with matrix is not open !\n";
		return 1;
	}
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < n; j++) {
			ifile >> A[i][j];
		}
	}
	ifile.close();

	ifile.open(vectorName);
	ifile >> b[0];
	if (!ifile.is_open()) {
		std::cerr << " Error : file with vector is not open !\n";
		return 1;
	}
	
	for (size_t i = 0; i < n; i++) {
		ifile >> b[i];
	}

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
	char outputName[maxLengthOfFileName];

	std::cout.precision(4);
	std::cout << "Считываем данные...\n";
	readSettings("D:\\settingsFile.txt", n, matrixName, vectorName, outputName);

	T** A = new T * [n];
	for (size_t i = 0; i < n; i++)	A[i] = new T[n];

	T* b = new T[n];

	input(n, A, b, matrixName, vectorName);

	std::cout << "\nИсходная СЛАУ:\n\n";
	showLinearSystem(n, A, b);

	T* x = new T[n];
	int er = gauss_method(n, A, b, x);
	if (er == -1) {
		std::cout << "detA = 0\nИсходная матрица вырождена!\n";
		return -1;
	}

	std::cout << "\nРешение системы методом Гаусса:\n";
	for (size_t i = 0; i < n; i++) {
		std::cout << std::setw(15) << x[i] << '\t';
	}
	std::cout << "\nНевязка решения: " << residual(n, A, b, x) << "\n\n";


	er =QR_method(n, A, b, x);
	if (er == -1) {
		std::cout << "detA = 0\nИсходная матрица вырождена!\n";
		return -1;
	}
	std::cout << "\nРешение системы методом QR:\n";
	for (size_t i = 0; i < n; i++) {
		std::cout << std::setw(15) << x[i] << '\t';
	}
	std::cout << "\nНевязка решения: " << residual(n, A, b, x) << "\n\n";





	T** Q = new T * [n];
	for (size_t i = 0; i < n; i++)	Q[i] = new T[n];

	T** R = new T * [n];
	for (size_t i = 0; i < n; i++) 	R[i] = new T[n];

	QR_decomposition(n, A, Q, R);

	std::cout << "\n\nОртогональная матрица вращений Q=T^-1=T^T:\n\n";
	showMatrix(n, Q);

	std::cout << "\nВерхне треугольная матрица R:\n\n";
	showMatrix(n, R);

	std::cout << "\nПроизведение матриц Q * R:\n\n";
	showMatrixMultiplication(n, Q, R);

	T** invA = inverse(n, A);
	std::cout << "\nОбратная матрица A^-1:\n\n";
	showMatrix(n, invA);


	std::cout << "\nПроизведение A * A^-1:\n\n";
	showMatrixMultiplication(n, A, invA);

	T norm_x = norm1(n, x);
	T norm_b = norm1(n, b);
	T conditionEstimate = 0;
	T* _b = new T[n];
	T* _x = new T[n];
	T* delta_x = new T[n];
	T* delta_b = new T[n];
	for (size_t i = 0; i < numOfDisturbances; i++) {

	

		for (size_t j = 0; j < n; j++) _b[j] = b[j];
		disturbance(n, _b);
		QR_solve(n, Q, R, _b, _x);
		for (size_t j = 0; j < n; j++) {
			delta_x[j] = x[j] - _x[j];
			delta_b[j] = b[j] - _b[j];
		}
		T cond = norm_b * norm1(n, delta_x) / (norm_x * norm1(n, delta_b));
		if (cond > conditionEstimate) {
			conditionEstimate = cond;
		}

	}
	delete[] _b;
	delete[] delta_x;
	delete[] delta_b;
	delete[] _x;
	
	std::cout << "\n\nОценка снизу числа обусловленности cond A = " << conditionEstimate << "\n";

	std::cout << "\n\nЧисло обусловленности для нормы ||*||_1 по определению cond A = ||A|| * ||A^-1|| = " << norm1(n, A) * norm1(n, invA)<<"\n\n";
	std::cout << "Число обусловленности для нормы ||*||_inf по определению cond A = ||A|| * ||A^-1|| = " << normInf(n, A) * normInf(n, invA) << "\n\n";


	for (size_t i = 0; i < n; i++) {
		delete[] A[i];
	}
	delete[] A;

	for (size_t i = 0; i < n; i++) {
		delete[] invA[i];
	}
	delete[] invA;

	delete[] b;
	delete[] x;

	for (size_t i = 0; i < n; i++)		delete[] Q[i];
	delete[] Q;

	for (size_t i = 0; i < n; i++) 		delete[] R[i];
	delete[] R;

#ifdef CHECK_MEMORY
	_CrtMemDumpAllObjectsSince(&_ms);
#endif

	system("pause");
	return 0;
}
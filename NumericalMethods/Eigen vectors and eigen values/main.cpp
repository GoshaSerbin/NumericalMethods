#include"header.h"


size_t iterationCounter;
size_t multiOperCounter;

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

//���� ������ �� ������
int input(const size_t& n, T** A,  char* const matrixName, char* const vectorName, char* const solveName) {

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
	std::cout << "��������� ������...\n";
	readSettings("D:\\settingsFile.txt", n, matrixName, vectorName, solveName);

	T** A = new T * [n];
	for (size_t i = 0; i < n; i++)	A[i] = new T[n];


	input(n, A, matrixName, vectorName, solveName);

	std::cout << "\n                           �������� �������:\n\n";
	showMatrix(n, A);

	T* eigenvalues = new T[n];

	T** eigenvectors = new T * [n];
	for (size_t i = 0; i < n; i++)	eigenvectors[i] = new T[n];

	

	std::cout << "\n\n                                 ����� QR-����������:\n";

	std::cout << "����������� ����� ������� ��� �������:\n";
	iterationCounter = 0;
	multiOperCounter = 0;
	eigenvaluesSearch(n, A, false,false, eigenvalues);
	for (size_t i = 0; i < n; i++) {
		std::cout << std::setw(15) << eigenvalues[i] << '\t';
	}
	std::cout << "\n���������� ��������: "<<iterationCounter<<"\n";
	std::cout << "\n���������� ���������: " << multiOperCounter << "\n";


	std::cout << "\n����������� ����� ������� c �������������� �������:\n";
	iterationCounter = 0;
	multiOperCounter = 0;
	eigenvaluesSearch(n, A, true,false, eigenvalues);
	for (size_t i = 0; i < n; i++) {
		std::cout << std::setw(15) << eigenvalues[i] << '\t';
	}
	std::cout << "\n���������� ��������: " << iterationCounter << "\n";
	std::cout << "\n���������� ���������: " << multiOperCounter << "\n";
	 
		std::cout << "\n\n                        ����� QR-���������� � ����������� � ������� �����������:\n";
		iterationCounter = 0;
		multiOperCounter = 0;
		std::cout << "����������� ����� ������� ��� �������:\n";
		eigenvaluesSearch(n, A, false,true, eigenvalues);
		for (size_t i = 0; i < n; i++) {
			std::cout << std::setw(15) << eigenvalues[i] << '\t';
		}
		std::cout << "\n���������� ��������: " << iterationCounter << "\n";
		std::cout << "\n���������� ���������: " << multiOperCounter << "\n";

		std::cout << "\n����������� ����� ������� c �������������� �������:\n";
		iterationCounter = 0;
		multiOperCounter = 0;
		eigenvaluesSearch(n, A, true,true, eigenvalues);
		for (size_t i = 0; i < n; i++) {
			std::cout << std::setw(15) << eigenvalues[i] << '\t';
		}
		std::cout << "\n���������� ��������: " << iterationCounter << "\n";
		std::cout << "\n���������� ���������: " << multiOperCounter << "\n";

		std::cout << "\n\n                        ����� �������� �������� ��� ���������� ����������� �.�.:\n";
		eigenvaluesSearch(n, A, true, false, eigenvalues);
		std::cout << "����������� ������� (�� ��������):\n";
		eigenvectorsSearch(n, A, eigenvalues, eigenvectors);
		showMatrix(n, eigenvectors);
		std::cout << "������:\n";
		T* dx = new T[n];
		for (size_t s = 0; s < n; s++) {
			T sum = 0;
			for (size_t i = 0; i < n; i++) {
				for (size_t j = 0; j < n; j++) {
					if (i == j)	sum += (A[i][j] - eigenvalues[s]) * eigenvectors[s][j];
					else sum += A[i][j] * eigenvectors[s][j];
				}
				dx[i] = sum;
			}
			std::cout << norm1(n, dx) << "\n";
		}
		//delete[] dx;

		std::cout << "\n\n                        ����������� ������ �������� ��������:\n";
		iterationCounter = 0;
		eigenSearch(n,A,eigenvalues,eigenvectors);
		std::cout << "����������� ����� �������:\n";
		for (size_t i = 0; i < n; i++) {
			std::cout << std::setw(15) << eigenvalues[i] << '\t';
		}
		std::cout << "\n����������� ������� (�� ��������):\n";
		showMatrix(n, eigenvectors);
		std::cout << "������:\n";
		for (size_t s = 0; s < n; s++) {
			T sum = 0;
			for (size_t i = 0; i < n; i++) {
				for (size_t j = 0; j < n; j++) {
					if (i == j)	sum += (A[i][j] - eigenvalues[s]) * eigenvectors[s][j];
					else sum += A[i][j] * eigenvectors[s][j];
				}
				dx[i] = sum;
			}
			std::cout << norm1(n, dx) << "\n";
		}
		delete[] dx;


	for (size_t i = 0; i < n; i++) {
		delete[] A[i];
	}
	delete[] A;



#ifdef CHECK_MEMORY
	_CrtMemDumpAllObjectsSince(&_ms);
#endif

	system("pause");
	return 0;
}
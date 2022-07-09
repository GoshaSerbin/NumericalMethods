#include"header.h"

T f1(T x1, T x2) {
//	return 2 - 2 * x2;
	return -4;

}


T G1(T x1) { //x2=0
//	return 1 + x1 * x1;
	return x1*x1;
	//return 1 - x1;
}

T G2(T x1) { //x2=L2
//	return -1 - x1 * x1;
	return 1+x1*x1;
}

T G3(T x2) { //x1=0
	return 0;
	return x2 * x2;
	return 1-x2;
}

T G4(T x2) { //x1=L1
	return 2;
	return 1 + x2 * x2;
	return 2-2*x2;
}

int main() {

#ifdef CHECK_MEMORY
	_CrtMemState _ms;
	_CrtMemCheckpoint(&_ms);
#endif

	setlocale(LC_ALL, "RUS");
	std::cout << "\n                                            Лабораторная работа №4 \n";
	T tmax = 2;	
	T L1 = 1;	
	T L2 = 1;	


	std::cout << "\n             ТЕСТ 1 \n";
	//DSolvePoissonEquation(L1, L2, f1, BC_1, G1, BC_1, G2, BC_2, G3, BC_2, G4);
	std::cout << "\n             ТЕСТ 2 \n";
	showTable1(L1, L2, f1, BC_1, G1, BC_1, G2, BC_2, G3, BC_2, G4);



#ifdef CHECK_MEMORY
	_CrtMemDumpAllObjectsSince(&_ms);
#endif

	system("pause");
	return 0;
}
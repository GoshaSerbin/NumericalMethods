#include"header.h"


void makeMeshGrid(size_t n, T a, T b, MESH_TYPE mt, T* mesh) {
	T h, c1, c2, c3;
	switch (mt) {
	case MT_UNIFORM:

		h = (b - a) / (1.*n);
		for (size_t i = 0; i <= n; i++) {
			mesh[i] = a + i * h;
		}
		break;
	case MT_CHEBYSHEV:

		c1 = (a + b) / 2.;
		c2 = (b - a) / 2.;
		c3 = PI / (2. * (n + 1.));

		for (size_t i = 0; i <= n; i++) {
			mesh[i] = c1 + c2 * cos((2. * i + 1) * c3);
		}

		for (size_t i = 0; i <= (n+1)/2; i++) {
			T mesh_i = mesh[i];
			mesh[i] = mesh[n - i];
			mesh[n - i] = mesh_i;
		}


		break;
	}


}

size_t makeMeshGrid(T h, T a, T b, T* mesh) {
	T x = a;
	size_t i;
	for (i = 0; x <= b; i++) {
		mesh[i] = x;
		x = x + i * h;
	}
	return i;

}

void showVector(size_t num, T* mesh) {
	for (size_t i = 0; i < num; i++) {
		std::cout << mesh[i] << " ";
	}
	std::cout << "\n";
}

int output(size_t n, T* x, T* y, char const* const fileName) {

	std::ofstream out;          // поток для записи
	out.open(fileName); // окрываем файл для записи
	if (out.is_open())
	{
		for (size_t i = 0; i < n; i++) {
			out << x[i] << " " << y[i] << "\n";
		}
	}
	out.close();
	return 0;
}

#pragma once
#include "myvector.h"
#include "matrix.h"

namespace slae
{
	struct SLAE
	{
		int n; //размерность СЛАУ
		matrix::Matrix A;
		myvector::MyVector b;//вектор правой части

		void initialize(int size, int unzero_matrix_elements);

		void reinitialize();
	};
}


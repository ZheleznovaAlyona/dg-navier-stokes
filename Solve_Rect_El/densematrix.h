#pragma once
#include <vector>
#include "myvector.h"

namespace densematrix
{
	class DenseMatrix
	{
	public:

		int n_lines, n_columns;
		std::vector <myvector::MyVector> ar;
	
		DenseMatrix();
		DenseMatrix(int size1, int size2);

		myvector::MyVector& operator[](int j);
		//умножение на вектор
		myvector::MyVector operator*(myvector::MyVector a);

		~DenseMatrix();

		void initialize(int size1, int size2);
	};
}
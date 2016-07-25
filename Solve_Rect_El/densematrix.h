#pragma once
#include <fstream>
#include <vector>
using namespace std;

#include "myvector.h"
using namespace myvector;


namespace densematrix
{
	class DenseMatrix
	{
	public:

		int n_lines, n_columns;
		vector <MyVector> ar;
	
		DenseMatrix();
		DenseMatrix(int size1, int size2);

		MyVector& operator[](int j);
		//умножение на вектор
		MyVector operator*(MyVector a);

		~DenseMatrix();

		void initialize(int size1, int size2);
	};
}
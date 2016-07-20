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
	
		DenseMatrix(){};

		DenseMatrix(int size1, int size2)
		{
			initialize(size1, size2);
		}

		MyVector& operator[](int j) 
		{
			return ar[j];
		}

		//умножение на вектор
		MyVector operator*(MyVector a) 
		{
			MyVector new_vector(n_lines);

			assert(a.ar.size() == n_columns);
			for(int j = 0; j < n_columns; j++)
				for(int i = 0; i < n_lines; i++)
					new_vector[i] += ar[j][i] * a[j];

			return new_vector;
		}

		~DenseMatrix(){};

		void initialize(int size1, int size2)
		{
			n_lines = size1; n_columns = size2;

			ar.reserve(n_columns);
			for(int i = 0; i < n_columns; i++)
				ar.push_back(MyVector (n_lines));
		}
	};
}
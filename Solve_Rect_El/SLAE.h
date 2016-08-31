#pragma once
#include "myvector.h"
#include "matrix.h"

namespace slae
{
	class SLAE
	{
	public:

		int n; //����������� ����
		matrix::Matrix A; 
		myvector::MyVector b;//������ ������ �����

		void initialize(int size, int unzero_matrix_elements);
		void reinitialize();
	};
}


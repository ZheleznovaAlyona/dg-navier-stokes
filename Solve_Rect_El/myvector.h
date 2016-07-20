#pragma once

#include <fstream>
#include <vector>
#include <assert.h>
#include <iomanip>
using namespace std;

namespace myvector
{

	class MyVector
	{
	public:

		vector <double> ar;

		MyVector();
		MyVector(int size);
		~MyVector();

		double& operator[](int j);	
		MyVector operator+(MyVector a);
		MyVector operator-(MyVector a);
		MyVector operator*(double a);
		MyVector operator/(double a);
		void initialize(int size);
		void make_zero();
		double norm();
	};

	std::ofstream& operator<<(std::ofstream& os, MyVector& vec);

}
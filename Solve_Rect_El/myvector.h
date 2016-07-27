#pragma once
#include <fstream>
#include <vector>

namespace myvector
{

	class MyVector
	{
	public:

		std::vector <double> ar;

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
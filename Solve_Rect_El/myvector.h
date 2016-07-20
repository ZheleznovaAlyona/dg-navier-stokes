#pragma once
#include <fstream>

namespace myvector
{

	class MyVector
	{
	public:

		vector <double> ar;

		MyVector(){};

		MyVector(int size)
		{
			ar.resize(size);
			memset(&ar[0], 0, size * sizeof(double)); //обнуляем
		}

		~MyVector(){};

		double& operator[](int j) 
		{
			return ar[j];
		}

	
		MyVector operator+(MyVector a) 
		{
			MyVector new_vector = MyVector(ar.size());
			assert(a.ar.size() == ar.size());
			for(int i = 0; i < ar.size(); i++)
				 new_vector.ar[i] = ar[i] + a[i];
			return new_vector;
		}

		MyVector operator-(MyVector a) 
		{
			MyVector new_vector = MyVector(ar.size());
			assert(a.ar.size() == ar.size());
			for(int i = 0; i < ar.size(); i++)
				 new_vector.ar[i] = ar[i] - a[i];
			return new_vector;
		}

		MyVector operator*(double a) 
		{
			MyVector new_vector = MyVector(ar.size());
			for(int i = 0; i < ar.size(); i++)
				 new_vector.ar[i] = ar[i] * a;
			return new_vector;
		}

		MyVector operator/(double a) 
		{
			MyVector new_vector = MyVector(ar.size());
			for(int i = 0; i < ar.size(); i++)
				 new_vector.ar[i] = ar[i] / a;
			return new_vector;
		}

		void initialize(int size)
		{
			ar.resize(size);
			memset(&ar[0], 0, size * sizeof(double)); //обнуляем
		}

		void make_zero()
		{
			memset(&ar[0], 0, ar.size() * sizeof(double)); //обнуляем
		}

		double norm()
		{
			double sum = 0;
			int size = ar.size();
			for(int i = 0; i < size; i++)
				sum += ar[i] * ar[i];

			return sqrt(sum);
		}
	};

	std::ofstream& operator<<(std::ofstream& os, MyVector& vec)
	{
		int size = vec.ar.size();

		for(int i = 0; i < size; i++)
			os << setprecision(20) << vec[i] << endl;

		return os;
	}

}
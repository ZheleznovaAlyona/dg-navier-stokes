#pragma once
#include <vector>
#include "myvector.h"
#include "partition.h"

namespace matrix
{
	class Matrix
	{
	public:
		//разреженный строчный формат
		int n;//размерность матрицы
		int size;//размерность массивов,где хранятся недиагональные элементы 

		std::vector <double> ggl;//массив с нижнетреугольными недиагональными элементами
		std::vector <double> ggu;//массив с верхнетреугольными недиагональными элементами
		std::vector <double> di;//диагональ
		std::vector <double> LU_ggu; //верхнетреугольные недиагональные элементы U
		std::vector <double> LU_ggl; //нижнетреугольные недиагональные элементы L
		std::vector <double> LU_di; //диагональные элементы L
		std::vector <int> ig;//указатели начала строк(столбцов)
		std::vector <int> jg;//номера столбцов(строк) внедиагональных элементов
		myvector::MyVector yl; //решение системы Lyl=F
		myvector::MyVector yu; //решение системы Uyu=F
	
		Matrix();

		Matrix(int size1, int size2);

		void initialize(int size1, int size2);
		void reinitialize();

		//умножение на вектор
		myvector::MyVector operator*(myvector::MyVector a);
		myvector::MyVector operator/(myvector::MyVector a);

		~Matrix();

		myvector::MyVector Uv(myvector::MyVector v);

		void LU();
		void LYF(myvector::MyVector b);
		void LYFt(myvector::MyVector b);
		void UXY(myvector::MyVector b);
		void UXYt(myvector::MyVector b);
		myvector::MyVector Uv_(myvector::MyVector v);

		void add_element(int i, int j, double element);
		void put_element(int i, int j, double element);
		void create_portret(partition::Partition& p);
	};


}
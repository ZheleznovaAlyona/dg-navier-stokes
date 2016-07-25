#pragma once
#include <fstream>
#include <vector>
#include <assert.h>
#include "myvector.h"
#include "myfunctions.h"
#include "testing_parameters.h"

using namespace myvector;
using namespace std;
using namespace testingparameters;

namespace matrix
{
	class Matrix
	{
	public:
		//разреженный строчный формат
		int n;//размерность матрицы
		int size;//размерность массивов,где хранятся недиагональные элементы 

		vector <double> ggl;//массив с нижнетреугольными недиагональными элементами
		vector <double> ggu;//массив с верхнетреугольными недиагональными элементами
		vector <double> di;//диагональ
		vector <double> LU_ggu; //верхнетреугольные недиагональные элементы U
		vector <double> LU_ggl; //нижнетреугольные недиагональные элементы L
		vector <double> LU_di; //диагональные элементы L
		MyVector b;//вектор правой части
		vector <int> ig;//указатели начала строк(столбцов)
		vector <int> jg;//номера столбцов(строк) внедиагональных элементов
		MyVector yl; //решение системы Lyl=F
		MyVector yu; //решение системы Uyu=F
	
		Matrix();

		Matrix(int size1, int size2);

		void initialize(int size1, int size2);
		void reinitialize();

		//умножение на вектор
		MyVector operator*(MyVector a);
		MyVector operator/(MyVector a);

		~Matrix();

		MyVector Uv(MyVector v);

		void LU();
		void LYF(MyVector b);
		void LYFt(MyVector b);
		void UXY(MyVector b);
		void UXYt(MyVector b);
		MyVector Uv_(MyVector v);
	};
}
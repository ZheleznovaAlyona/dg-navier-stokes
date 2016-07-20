#pragma once
#include <fstream>

#include "functions.cpp"

using namespace myvector;

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
		//vector <double> LU_ggu; //верхнетреугольные недиагональные элементы U
		//vector <double> LU_ggl; //нижнетреугольные недиагональные элементы L
		//vector <double> LU_di; //диагональные элементы L
		MyVector b;//вектор правой части
		vector <int> ig;//указатели начала строк(столбцов)
		vector <int> jg;//номера столбцов(строк) внедиагональных элементов

	
		Matrix(){};

		Matrix(int size1, int size2)
		{
			initialize(size1, size2);
		}

		void initialize(int size1, int size2)
		{
			n = size1; size = size2;

			initialize_vector(ggl, size);
			initialize_vector(ggu, size);
			initialize_vector(di, n);
			b.initialize(n);
			initialize_vector(ig, n + 1);
			initialize_vector(jg, size);
		};
		void reinitialize()
		{
			memset(&ggl[0], 0, size * sizeof(double)); //обнуляем
			memset(&ggu[0], 0, size * sizeof(double)); //обнуляем
			memset(&di[0], 0, n * sizeof(double)); //обнуляем
			b.make_zero();
		};

		//умножение на вектор
		MyVector operator*(MyVector a) 
		{
			int i, j, k, kol;
			int iend;
			MyVector new_vector = MyVector(a.ar.size());

			assert(a.ar.size() == n);
			for(i = 0; i < n; i++)
			{
				kol = ig[i + 1] - ig[i];//количество ненулевых элементов строки (столбца) от первого
									  //ненулевого элемента до диагонального элемента (не включая его)
				iend = ig[i + 1];
				k = ig[i]; // адрес первого занятого элемента строки (столбца) 

				new_vector[i] = di[i] * a[i];//от главной диагонали

				for(; k < iend; k++)//проходим по всем элементам i строки (столбца)
				{
					j = jg[k];
					new_vector[i] += ggl[k] * a[j];//от нижнего треугольника
					new_vector[j] += ggu[k] * a[i];//от верхнего треугольника
				}
			}

			return new_vector;
		}

		MyVector operator/(MyVector a) 
		{
			int i, j, k, kol;
			int iend;
			MyVector new_vector = MyVector(a.ar.size());

			assert(a.ar.size() == n);
			for(i = 0; i < n; i++)
			{
				kol = ig[i + 1] - ig[i];//количество ненулевых элементов строки (столбца) от первого
									  //ненулевого элемента до диагонального элемента (не включая его)
				iend = ig[i + 1];
				k = ig[i]; // адрес первого занятого элемента строки (столбца) 

				new_vector[i] = di[i] * a[i];//от главной диагонали

				for(; k < iend; k++)//проходим по всем элементам i строки (столбца)
				{
					j = jg[k];
					new_vector[i] += ggu[k] * a[j];//от нижнего треугольника
					new_vector[j] += ggl[k] * a[i];//от верхнего треугольника
				}
			}

			return new_vector;
		}

		~Matrix(){};

		MyVector Uv(MyVector v)
		{
			int i, j, k, kol;
			int iend;
			MyVector new_vector = MyVector(v.ar.size());

			assert(v.ar.size() == n);
			for(i = 0; i < n; i++)
			{
				kol = ig[i+1] - ig[i];//количество ненулевых элементов столбца от первого
										//ненулевого элемента до диагонального элемента (не включая его)
				iend = ig[i+1];
				k = ig[i]; // адрес первого занятого элемента столбца

				new_vector[i] = v[i];//от главной диагонали (у U на диагонали 1)

				for(; k < iend; k++)//проходим по всем элементам i столбца
				{
					j = jg[k];
					new_vector[j] += ggu[k] * v[i];//от верхнего треугольника
				}
			}

			return new_vector;
		}

	};
}
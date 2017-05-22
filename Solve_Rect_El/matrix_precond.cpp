#include "matrix.h"
#include <assert.h>
#include <algorithm>
#include "myfunctions.h"
#include "testing_parameters.h"
#include "boundaries.h"
#include "boundary_conditions.h"
#include <iostream>

using namespace myvector;
using namespace std;
using namespace testingparameters;
using namespace boundary_conditions;
using namespace partition;
using namespace boundaries;
using namespace logger;

namespace matrix
{
	void Matrix::LU()
	{
		int i;
		int i0, j0;
		int iend;
		int num;
		int ki, kj;
		double suml, sumu, sumdg;
		int size2 = size;

		for (i = 0; i < size2; i++)
		{
			LU_ggu[i] = ggu[i];
			LU_ggl[i] = ggl[i];
		}

		for (i = 0; i < n; i++)
			LU_di[i] = di[i];

		for (i = 0; i < n; i++)
		{
			i0 = ig[i];
			iend = ig[i + 1];
			for (num = i0, sumdg = 0; num < iend; num++)
			{
				j0 = ig[jg[num]]; //в зависимости от номера фиксируем столбец,какой столбец l,такого столбца  ищем начальный эл у u 
				int jend = ig[jg[num] + 1];
				ki = i0;
				kj = j0;
				for (suml = 0, sumu = 0, ki = i0; ki < num; ki++) //для num учитываются все предыдущие элементы
				{
					for (int m = kj; m < jend; m++)
						if (jg[ki] == jg[m]) //ищем соответствующие ненулевые элементы для умножения
						{
							suml += LU_ggl[ki] * LU_ggu[m];
							sumu += LU_ggl[m] * LU_ggu[ki];//для симметричного элемента из U
						}
				}
				LU_ggl[num] -= suml;
				LU_ggu[num] = (LU_ggu[num] - sumu) / LU_di[jg[num]];
				sumdg += LU_ggl[num] * LU_ggu[num];//умножаются симметричные элементы	
			}
			LU_di[i] -= sumdg;
		}
	}
	void Matrix::LYF(MyVector& b)
	{
		int i, k;
		int i0;//адрес начала строки
		int iend;//адрес конца строки
		double sum;

		yl.make_zero();

		if (Testing_parameters::use_LU)
		{
			for (i = 0; i < n; i++)
			{
				i0 = ig[i]; iend = ig[i + 1];

				for (k = i0, sum = 0; k < iend; k++)
					sum += yl[jg[k]] * LU_ggl[k];

				yl[i] = (b[i] - sum) / LU_di[i];
			}
		}
		else
			yl = b;
	}
	void Matrix::LYFt(MyVector& b)
	{
		int i, k;
		int i0;//адрес начала строки
		int iend;//адрес конца строки
		double sum;

		yl.make_zero();
		if (Testing_parameters::use_LU)
		{
			MyVector bb(n);
			bb = b;
			for (i = n - 1; i >= 0; i--)
			{
				i0 = ig[i]; iend = ig[i + 1];
				yl[i] = bb[i] / LU_di[i];
				for (k = i0, sum = 0; k < iend; k++)
					bb[jg[k]] -= yl[i] * LU_ggl[k];
			}
		}
		else
			yl = b;
	}
	void Matrix::UXY(MyVector& b)
	{
		int i, k;
		int i0;
		int iend;

		yu.make_zero();
		if (Testing_parameters::use_LU)
		{
			for (i = n - 1; i >= 0; i--)//проход по столбцам с конца
			{
				yu[i] += b[i];
				i0 = ig[i]; iend = ig[i + 1];

				for (k = iend - 1; k >= i0; k--)//идём по столбцу с конца
					yu[jg[k]] -= yu[i] * LU_ggu[k];
			}
		}
		else
			yu = b;
	}
	void Matrix::UXYt(MyVector& b)
	{
		int i, k;
		int i0;
		int iend;

		yu.make_zero();
		if (Testing_parameters::use_LU)
		{
			for (i = n - 1; i >= 0; i--)//проход по столбцам с конца
			{
				yu[i] += b[i];
				i0 = ig[i]; iend = ig[i + 1];

				for (k = iend - 1; k >= i0; k--)//идём по столбцу с конца
					yu[i] -= yu[jg[k]] * LU_ggu[k];
			}
		}
		else
			yu = b;
	}

	MyVector Matrix::Uv_(MyVector& v)
	{
		int i, j, k, kol;
		int iend;
		if (Testing_parameters::use_LU)
		{ 
			MyVector new_vector = MyVector(v.ar.size());

			assert(v.ar.size() == n);
			return v;
			for (i = 0; i < n; i++)
			{
				kol = ig[i + 1] - ig[i];//количество ненулевых элементов столбца от первого
										//ненулевого элемента до диагонального элемента (не включая его)
				iend = ig[i + 1];
				k = ig[i]; // адрес первого занятого элемента столбца

				new_vector[i] = v[i];//от главной диагонали (у U на диагонали 1)

				for (; k < iend; k++)//проходим по всем элементам i столбца
				{
					j = jg[k];
					new_vector[j] += LU_ggu[k] * v[i];//от верхнего треугольника
				}
			}

			return new_vector;
		}

		return v;
	}
}
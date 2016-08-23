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

namespace matrix
{	
	Matrix::Matrix(){};

	Matrix::Matrix(int size1, int size2)
	{
		initialize(size1, size2);
	}

	void Matrix::initialize(int size1, int size2)
	{		
		n = size1; size = size2;

		initialize_vector(ggl, size);
		initialize_vector(ggu, size);
		initialize_vector(di, n);
		initialize_vector(ig, n + 1);
		initialize_vector(jg, size);
		yl.initialize(n);
		yu.initialize(n);

		if (Testing_parameters::use_LU)
		{
			initialize_vector(LU_ggl, size);
			initialize_vector(LU_ggu, size);
			initialize_vector(LU_di, n);
		}
	};

	void Matrix::reinitialize()
	{
		memset(&ggl[0], 0, size * sizeof(double)); //обнуляем
		memset(&ggu[0], 0, size * sizeof(double)); //обнуляем
		memset(&di[0], 0, n * sizeof(double)); //обнуляем

		if (Testing_parameters::use_LU)
		{
			memset(&LU_ggl[0], 0, size * sizeof(double));
			memset(&LU_ggu[0], 0, size * sizeof(double));
			memset(&LU_di[0], 0, n * sizeof(double));
		}
	};

	MyVector Matrix::operator*(MyVector a) 
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

	MyVector Matrix::operator/(MyVector a) 
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

	Matrix::~Matrix(){};

	MyVector Matrix::Uv(MyVector v)
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

	void Matrix::LU()
	{
		int i;
		int i0,j0;
		int iend;
		int num;
		int ki,kj;
		double suml,sumu,sumdg;
		int size2 = size;

		for(i = 0; i < size2; i++) 
		{
			LU_ggu[i] = ggu[i];
			LU_ggl[i] = ggl[i];
		}

		for(i = 0; i < n; i++) 
			LU_di[i] = di[i];

		for(i = 0; i < n; i++)
		{
			i0 = ig[i];
			iend = ig[i + 1];
			for(num = i0,sumdg = 0; num < iend; num++)
			{
				j0 = ig[jg[num]]; //в зависимости от номера фиксируем столбец,какой столбец l,такого столбца  ищем начальный эл у u 
				int jend=ig[jg[num] + 1];
				ki=i0;
				kj=j0;
				for(suml = 0, sumu = 0, ki = i0; ki < num; ki++) //для num учитываются все предыдущие элементы
				{
					for(int m = kj; m < jend; m++)
					if(jg[ki] == jg[m]) //ищем соответствующие ненулевые элементы для умножения
					{
						suml += LU_ggl[ki] * LU_ggu[m];
						sumu += LU_ggl[m] * LU_ggu[ki];//для симметричного элемента из U
					}
				}
				LU_ggl[num] -= suml;	
				LU_ggu[num] = (LU_ggu[num] - sumu) / LU_di[jg[num]];
			sumdg += LU_ggl[num] * LU_ggu[num];//умножаются симметричные элементы	
			}
			LU_di[i]-=sumdg;
		}
	}
	void Matrix::LYF(MyVector b)
	{
		int i, k;
		int i0;//адрес начала строки
		int iend;//адрес конца строки
		double sum;

		yl.make_zero();

		if(Testing_parameters::use_LU)
		{
			for(i = 0; i < n; i++)
			{
				i0 = ig[i]; iend = ig[i+1];

				for(k = i0, sum = 0; k < iend; k++)
					sum += yl[jg[k]] * LU_ggl[k];

				yl[i] = (b[i] - sum) / LU_di[i];
			}
		}
		else
			yl = b;
	}
	void Matrix::LYFt(MyVector b)
	{
		int i, k;
		int i0;//адрес начала строки
		int iend;//адрес конца строки
		double sum;

		yl.make_zero();
		if(Testing_parameters::use_LU)
		{
			MyVector bb(n);
			bb = b;
			for(i = n - 1; i >= 0; i--)
			{
				i0 = ig[i]; iend = ig[i + 1];
				yl[i] = bb[i] / LU_di[i];
				for(k = i0, sum = 0; k < iend; k++)
					bb[jg[k]] -= yl[i] * LU_ggl[k];
			}
		}
		else
			yl = b;
	}
	void Matrix::UXY(MyVector b)
	{
		int i, k;
		int i0;
		int iend;

		yu.make_zero();
		if(Testing_parameters::use_LU)
		{
			for(i = n - 1; i >= 0; i--)//проход по столбцам с конца
			{
				yu[i] += b[i];
				i0 = ig[i]; iend = ig[i + 1];

				for(k = iend - 1; k >= i0; k--)//идём по столбцу с конца
					yu[jg[k]] -= yu[i] * LU_ggu[k];
			}
		}
		else
			yu = b;
	}
	void Matrix::UXYt(MyVector b)
	{
		int i, k;
		int i0;
		int iend;

		yu.make_zero();
		if(Testing_parameters::use_LU)
		{
			for(i = n - 1; i >= 0; i--)//проход по столбцам с конца
			{
				yu[i] += b[i];
				i0 = ig[i]; iend = ig[i + 1];

				for(k = iend - 1; k >= i0; k--)//идём по столбцу с конца
					yu[i] -= yu[jg[k]] * LU_ggu[k];
			}
		}
		else
			yu = b;
	}

	MyVector Matrix::Uv_(MyVector v)
	{
		int i, j, k, kol;
		int iend;
		MyVector new_vector = MyVector(v.ar.size());

		assert(v.ar.size() == n);
		return v;
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
				new_vector[j] += LU_ggu[k] * v[i];//от верхнего треугольника
			}
		}

		return new_vector;
	}

	void Matrix::add_element(int i, int j, double element)
	{
		int id;
		bool flag;

		if(i == j)
			di[i] += element;
		else
		{
			if(i < j)
			{	
				flag = false;
				for(id = ig[j]; !flag && id <= ig[j + 1] - 1; id++)
					if(jg[id] == i) flag = true;
				 if(flag) ggu[id - 1] += element;
			}
			else
			{
				flag = false;
				for(id = ig[i]; !flag && id <= ig[i + 1] - 1; id++)
					if(jg[id] == j) flag = true;
				if(flag) ggl[id - 1] += element;
			}
		}
	}

	void Matrix::put_element(int i, int j, double element)
	{
		int id;
		bool flag;

		if(i == j)
			di[i] = element;
		else
		{
			if(i < j)
			{	
				flag = false;
				for(id = ig[j]; !flag && id <= ig[j + 1] - 1; id++)
					if(jg[id] == i) flag = true;
				 if(flag) ggu[id - 1] = element;
			}
			else
			{
				flag = false;
				for(id = ig[i]; !flag && id <= ig[i + 1] - 1; id++)
					if(jg[id] == j) flag = true;
				if(flag) ggl[id - 1] = element;
			}
		}
	}

	void Matrix::create_portret(Partition& p)
	{
		cout << "Creating matrix-portret..." << endl;
		vector <int> unzero_elements_list;
		vector <int> *lists;	
		int unzero_elements_lists_size; 
		int current_number;
		int n_edges = p.elements.size() * 4;

		lists = new vector <int>[n];
		int count_elements = p.elements.size();

		for(int i = 0; i < count_elements; i++)
		{
			//общий принцип сборки портрета для uu, pp, pu
			//--------------------------------------------
			//1. собираем в список:
			//*для к.э. глобальные номера dof по i и затем
			//*ненулевые для к.э. глобальные номера dof по j;
			//2. идём по первым элементам списка (dof по i)
			//и выбираем для каждого номера (для p номер=индекс цикла + число рёбер),
			//меньшие его, т.е. те, которые будут располагаться левее соответствующей 
			//диагонали, потому что портрет строим по строкам
			//а затем кладём в соответствующий список
			//2.a. сортируем список по возрастанию

			//структура СЛАУ:
			/*
			 ---------
			| UU | UP |
			-----|----
			| PU | PP |
			 ---------
			*/

			//блок UU
			//1
			unzero_elements_lists_size = p.create_unzero_elements_list(i, 
																	 unzero_elements_list, 
																	 4, 
																	 4, 
																	 p.elements[i].edges, 
																	 p.elements[i].edges, 
																	 true);
			//2
			for(int j = 0; j < 4; j++)
			{
				current_number = unzero_elements_list[j];
				for(int k = 4; k < unzero_elements_lists_size; k++)
					if(unzero_elements_list[k] < current_number)
						lists[current_number].push_back(unzero_elements_list[k]);
					// 2.a
					sort(lists[current_number].begin(), lists[current_number].end());
			}
			unzero_elements_list.clear();

			//блок PP
			//1
			unzero_elements_lists_size = p.create_unzero_elements_list(i, 
																	 unzero_elements_list, 
																	 4, 
																	 4, 
																	 p.elements[i].nodes, 
																	 p.elements[i].nodes, 
																	 false);
			//2
			for(int j = 0; j < 4; j++)
			{
				current_number = unzero_elements_list[j];
				for(int k = 4; k < unzero_elements_lists_size; k++)
					if(unzero_elements_list[k] < current_number)
						lists[current_number + n_edges].push_back(unzero_elements_list[k] + n_edges);	
					//2.a 
					//можно не сортировать, потому что потом всё равно добавятся ещё элементы из PU
			}
			unzero_elements_list.clear();

			//блок PU
			//1
			unzero_elements_lists_size = p.create_unzero_elements_list(i,
																	 unzero_elements_list, 
																	 4, 
																	 4, 
																	 p.elements[i].nodes, 
																	 p.elements[i].edges, 
																	 true);
			//2
			for(int j = 0; j < 4; j++)
			{
				current_number = unzero_elements_list[j];
				for(int k = 4; k < unzero_elements_lists_size; k++)
					if(unzero_elements_list[k] < current_number + n_edges)
						lists[current_number + n_edges].push_back(unzero_elements_list[k]);	
					// 2.a
					sort(lists[current_number + n_edges].begin(), 
						lists[current_number + n_edges].end());
			}
			unzero_elements_list.clear();
		}

		ig[0] = 0;

		for(int i = 0; i < n; i++)
		{
			if(!lists[i].empty())
				ig[i + 1] = ig[i] + lists[i].size();
			else ig[i + 1] = ig[i];
		}

		int k = 0;
		for(int i = 0; i < n; i++)
		{
			if(!lists[i].empty())
			{
				for(unsigned int j = 0; j < lists[i].size(); j++)
				{
					jg[k] = lists[i][j];
					k++;
				}
				lists[i].clear();
			}
		}

		delete[] lists;
	}
}
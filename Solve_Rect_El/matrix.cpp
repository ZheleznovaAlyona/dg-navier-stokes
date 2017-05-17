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

	MyVector Matrix::operator*(MyVector& a) 
	{
		int i, j, k, kol;
		int iend;
		MyVector new_vector(a.ar.size());

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

	MyVector Matrix::operator/(MyVector& a)
	{
		int i, j, k, kol;
		int iend;
		MyVector new_vector(a.ar.size());

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

	MyVector Matrix::Uv(MyVector& v)
	{
		int i, j, k, kol;
		int iend;
		MyVector new_vector(v.ar.size());

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

	void Matrix::create_portret(Partition& p, Logger& log)
	{
		log.send_message_create_portret();
		vector <int> unzero_elements_list;
		vector<vector <int>> lists; lists.resize(n);
		int unzero_elements_lists_size; 
		int current_number;
		int n_func_u_ = p.elements[0].ndof_u;
		int n_func_u_all = p.elements.size() * n_func_u_;

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
																	   p.elements[i].ndof_u, 
																	   p.elements[i].ndof_u, 
																	   p.elements[i].dof_u, 
																	   p.elements[i].dof_u, 
																	   true);
			//2
			for(int j = 0; j < p.elements[i].ndof_u; j++)
			{
				current_number = unzero_elements_list[j];
				for(int k = p.elements[i].ndof_u; k < unzero_elements_lists_size; k++)
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
																	 p.elements[i].ndof_p, 
																	 p.elements[i].ndof_p, 
																	 p.elements[i].dof_p, 
																	 p.elements[i].dof_p, 
																	 false);
			//2
			for(int j = 0; j < p.elements[i].ndof_p; j++)
			{
				current_number = unzero_elements_list[j];
				for(int k = p.elements[i].ndof_p; k < unzero_elements_lists_size; k++)
					if(unzero_elements_list[k] < current_number)
						lists[current_number].push_back(unzero_elements_list[k]);	
					//2.a 
					//можно не сортировать, потому что потом всё равно добавятся ещё элементы из PU
			}
			unzero_elements_list.clear();

			//блок PU
			//1
			unzero_elements_lists_size = p.create_unzero_elements_list(i,
																	   unzero_elements_list, 
																	   p.elements[i].ndof_p, 
																	   p.elements[i].ndof_u, 
																	   p.elements[i].dof_p, 
																	   p.elements[i].dof_u, 
																	   true);
			//2
			for(int j = 0; j < p.elements[i].ndof_p; j++)
			{
				current_number = unzero_elements_list[j];
				for(int k = p.elements[i].ndof_p; k < unzero_elements_lists_size; k++)
					if(unzero_elements_list[k] < current_number)
						lists[current_number].push_back(unzero_elements_list[k]);	
					// 2.a
					sort(lists[current_number].begin(), 
						lists[current_number].end());
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
	}
}
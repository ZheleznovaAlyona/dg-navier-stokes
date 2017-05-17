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
		memset(&ggl[0], 0, size * sizeof(double)); //��������
		memset(&ggu[0], 0, size * sizeof(double)); //��������
		memset(&di[0], 0, n * sizeof(double)); //��������

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
			kol = ig[i + 1] - ig[i];//���������� ��������� ��������� ������ (�������) �� �������
									//���������� �������� �� ������������� �������� (�� ������� ���)
			iend = ig[i + 1];
			k = ig[i]; // ����� ������� �������� �������� ������ (�������) 

			new_vector[i] = di[i] * a[i];//�� ������� ���������

			for(; k < iend; k++)//�������� �� ���� ��������� i ������ (�������)
			{
				j = jg[k];
				new_vector[i] += ggl[k] * a[j];//�� ������� ������������
				new_vector[j] += ggu[k] * a[i];//�� �������� ������������
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
			kol = ig[i + 1] - ig[i];//���������� ��������� ��������� ������ (�������) �� �������
									//���������� �������� �� ������������� �������� (�� ������� ���)
			iend = ig[i + 1];
			k = ig[i]; // ����� ������� �������� �������� ������ (�������) 

			new_vector[i] = di[i] * a[i];//�� ������� ���������

			for(; k < iend; k++)//�������� �� ���� ��������� i ������ (�������)
			{
				j = jg[k];
				new_vector[i] += ggu[k] * a[j];//�� ������� ������������
				new_vector[j] += ggl[k] * a[i];//�� �������� ������������
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
			kol = ig[i+1] - ig[i];//���������� ��������� ��������� ������� �� �������
									//���������� �������� �� ������������� �������� (�� ������� ���)
			iend = ig[i+1];
			k = ig[i]; // ����� ������� �������� �������� �������

			new_vector[i] = v[i];//�� ������� ��������� (� U �� ��������� 1)

			for(; k < iend; k++)//�������� �� ���� ��������� i �������
			{
				j = jg[k];
				new_vector[j] += ggu[k] * v[i];//�� �������� ������������
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
			//����� ������� ������ �������� ��� uu, pp, pu
			//--------------------------------------------
			//1. �������� � ������:
			//*��� �.�. ���������� ������ dof �� i � �����
			//*��������� ��� �.�. ���������� ������ dof �� j;
			//2. ��� �� ������ ��������� ������ (dof �� i)
			//� �������� ��� ������� ������ (��� p �����=������ ����� + ����� ����),
			//������� ���, �.�. ��, ������� ����� ������������� ����� ��������������� 
			//���������, ������ ��� ������� ������ �� �������
			//� ����� ����� � ��������������� ������
			//2.a. ��������� ������ �� �����������

			//��������� ����:
			/*
			 ---------
			| UU | UP |
			-----|----
			| PU | PP |
			 ---------
			*/

			//���� UU
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

			//���� PP
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
					//����� �� �����������, ������ ��� ����� �� ����� ��������� ��� �������� �� PU
			}
			unzero_elements_list.clear();

			//���� PU
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
#pragma once
#include <fstream>

#include "functions.cpp"

using namespace myvector;

namespace matrix
{
	class Matrix
	{
	public:

		//����������� �������� ������
		int n;//����������� �������
		int size;//����������� ��������,��� �������� �������������� �������� 

		vector <double> ggl;//������ � ����������������� ��������������� ����������
		vector <double> ggu;//������ � ������������������ ��������������� ����������
		vector <double> di;//���������
		//vector <double> LU_ggu; //����������������� �������������� �������� U
		//vector <double> LU_ggl; //���������������� �������������� �������� L
		//vector <double> LU_di; //������������ �������� L
		MyVector b;//������ ������ �����
		vector <int> ig;//��������� ������ �����(��������)
		vector <int> jg;//������ ��������(�����) ��������������� ���������

	
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
			memset(&ggl[0], 0, size * sizeof(double)); //��������
			memset(&ggu[0], 0, size * sizeof(double)); //��������
			memset(&di[0], 0, n * sizeof(double)); //��������
			b.make_zero();
		};

		//��������� �� ������
		MyVector operator*(MyVector a) 
		{
			int i, j, k, kol;
			int iend;
			MyVector new_vector = MyVector(a.ar.size());

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

		MyVector operator/(MyVector a) 
		{
			int i, j, k, kol;
			int iend;
			MyVector new_vector = MyVector(a.ar.size());

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

		~Matrix(){};

		MyVector Uv(MyVector v)
		{
			int i, j, k, kol;
			int iend;
			MyVector new_vector = MyVector(v.ar.size());

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

	};
}
#pragma once
#include <fstream>
#include <vector>
#include "myvector.h"
#include "myfunctions.h"

using namespace myvector;
using namespace std;

bool use_LU;
int test = 3;
int solver = 2;


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
		vector <double> LU_ggu; //����������������� �������������� �������� U
		vector <double> LU_ggl; //���������������� �������������� �������� L
		vector <double> LU_di; //������������ �������� L
		MyVector b;//������ ������ �����
		vector <int> ig;//��������� ������ �����(��������)
		vector <int> jg;//������ ��������(�����) ��������������� ���������
		MyVector yl; //������� ������� Lyl=F
		MyVector yu; //������� ������� Uyu=F

	
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

		void LU()
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
					j0 = ig[jg[num]]; //� ����������� �� ������ ��������� �������,����� ������� l,������ �������  ���� ��������� �� � u 
					int jend=ig[jg[num] + 1];
					ki=i0;
					kj=j0;
					for(suml = 0, sumu = 0, ki = i0; ki < num; ki++) //��� num ����������� ��� ���������� ��������
					{
						for(int m = kj; m < jend; m++)
						if(jg[ki] == jg[m]) //���� ��������������� ��������� �������� ��� ���������
						{
							suml += LU_ggl[ki] * LU_ggu[m];
							sumu += LU_ggl[m] * LU_ggu[ki];//��� ������������� �������� �� U
						}
					}
					LU_ggl[num] -= suml;	
					LU_ggu[num] = (LU_ggu[num] - sumu) / LU_di[jg[num]];
				sumdg += LU_ggl[num] * LU_ggu[num];//���������� ������������ ��������	
				}
				LU_di[i]-=sumdg;
			}
		}
		void LYF(MyVector b)
		{
			int i, k;
			int i0;//����� ������ ������
			int iend;//����� ����� ������
			double sum;

			yl.make_zero();

			if(use_LU)
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
		void LYFt(MyVector b)
		{
			int i, k;
			int i0;//����� ������ ������
			int iend;//����� ����� ������
			double sum;

			yl.make_zero();
			if(use_LU)
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
		void UXY(MyVector b)
		{
			int i, k;
			int i0;
			int iend;

			yu.make_zero();
			if(use_LU)
			{
				for(i = n - 1; i >= 0; i--)//������ �� �������� � �����
				{
					yu[i] += b[i];
					i0 = ig[i]; iend = ig[i + 1];

					for(k = iend - 1; k >= i0; k--)//��� �� ������� � �����
						yu[jg[k]] -= yu[i] * LU_ggu[k];
				}
			}
			else
				yu = b;
		}
		void UXYt(MyVector b)
		{
			int i, k;
			int i0;
			int iend;

			yu.make_zero();
			if(use_LU)
			{
				for(i = n - 1; i >= 0; i--)//������ �� �������� � �����
				{
					yu[i] += b[i];
					i0 = ig[i]; iend = ig[i + 1];

					for(k = iend - 1; k >= i0; k--)//��� �� ������� � �����
						yu[i] -= yu[jg[k]] * LU_ggu[k];
				}
			}
			else
				yu = b;
		}
		MyVector Uv_(MyVector v)
		{
			int i, j, k, kol;
			int iend;
			MyVector new_vector = MyVector(v.ar.size());

			assert(v.ar.size() == n);
			return v;
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
					new_vector[j] += LU_ggu[k] * v[i];//�� �������� ������������
				}
			}

			return new_vector;
		}
	};
}
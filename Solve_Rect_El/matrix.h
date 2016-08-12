#pragma once
#include <vector>
#include "myvector.h"
#include "partition.h"

namespace matrix
{
	class Matrix
	{
	public:
		//����������� �������� ������
		int n;//����������� �������
		int size;//����������� ��������,��� �������� �������������� �������� 

		std::vector <double> ggl;//������ � ����������������� ��������������� ����������
		std::vector <double> ggu;//������ � ������������������ ��������������� ����������
		std::vector <double> di;//���������
		std::vector <double> LU_ggu; //����������������� �������������� �������� U
		std::vector <double> LU_ggl; //���������������� �������������� �������� L
		std::vector <double> LU_di; //������������ �������� L
		std::vector <int> ig;//��������� ������ �����(��������)
		std::vector <int> jg;//������ ��������(�����) ��������������� ���������
		myvector::MyVector yl; //������� ������� Lyl=F
		myvector::MyVector yu; //������� ������� Uyu=F
	
		Matrix();

		Matrix(int size1, int size2);

		void initialize(int size1, int size2);
		void reinitialize();

		//��������� �� ������
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
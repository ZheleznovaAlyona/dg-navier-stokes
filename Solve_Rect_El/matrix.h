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
	
		Matrix();

		Matrix(int size1, int size2);

		void initialize(int size1, int size2);
		void reinitialize();

		//��������� �� ������
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
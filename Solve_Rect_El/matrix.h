#pragma once
#include <vector>
#include "myvector.h"
#include "partition.h"
#include "boundaries.h"
#include "boundary_conditions.h"

namespace matrix
{
	class MatrixSupport
	{
		
	};

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
		myvector::MyVector b;//������ ������ �����
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

		int count_unzero_matrix_elements(partition::Partition& p);
		int create_unzero_elements_list(int element_number, 
										vector <int> &list, 
										int dof_num_i, 
										int dof_num_j, 
										int *dof_i, 
										int *dof_j,
										bool dof_j_edge,
										partition::Partition& p);
		void create_portret(partition::Partition& p);

		void calculate_global_matrix(MyVector q_calc, 
									 partition::Partition& p,
									 boundaries::InternalBoundaries& internal_bs,
									 boundaries::OuterBoundaries& outer_bs,
									 boundary_conditions::BoundaryConditionsSupport&
										b_conditions);

	};


}
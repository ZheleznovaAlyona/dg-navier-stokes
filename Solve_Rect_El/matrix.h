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
		//разреженный строчный формат
		int n;//размерность матрицы
		int size;//размерность массивов,где хранятся недиагональные элементы 

		std::vector <double> ggl;//массив с нижнетреугольными недиагональными элементами
		std::vector <double> ggu;//массив с верхнетреугольными недиагональными элементами
		std::vector <double> di;//диагональ
		std::vector <double> LU_ggu; //верхнетреугольные недиагональные элементы U
		std::vector <double> LU_ggl; //нижнетреугольные недиагональные элементы L
		std::vector <double> LU_di; //диагональные элементы L
		myvector::MyVector b;//вектор правой части
		std::vector <int> ig;//указатели начала строк(столбцов)
		std::vector <int> jg;//номера столбцов(строк) внедиагональных элементов
		myvector::MyVector yl; //решение системы Lyl=F
		myvector::MyVector yu; //решение системы Uyu=F
	
		Matrix();

		Matrix(int size1, int size2);

		void initialize(int size1, int size2);
		void reinitialize();

		//умножение на вектор
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
#pragma once
#include <fstream>
#include <vector>
#include "partition.h"
#include "parameters.h"
#include "myvector.h"
#include "integration.h"
#include "basis.h"

namespace boundary_conditions
{
	class BoundaryCondition
	{
	public:

		int elem;
		int edges[4]; //�����,������,������, �������: 1 - ����, 0 - ���
		int formula_number;
	};

	class BoundaryConditionsSupport : virtual public partition::Partition,
									  virtual public parameters::Parameters,
									  virtual public integration::Gauss_integration,
									  virtual public basis::Basis
	{
	public:

		std::vector <BoundaryCondition> boundaries1; //������ ������� �������
		std::vector <BoundaryCondition> boundaries2; //������ -//-
		std::vector <BoundaryCondition> boundaries3; //������ -//-

		double mu1;

		void initialize_penalty_parameters();

		void calculate_all_boundaries1(myvector::MyVector& b);
		void calculate_boundaries1(int number, myvector::MyVector& b);

		void calculate_boundaries1_left(int number, myvector::MyVector& b);
		void calculate_boundaries1_right(int number, myvector::MyVector& b);
		void calculate_boundaries1_low(int number, myvector::MyVector& b);
		void calculate_boundaries1_up(int number, myvector::MyVector& b);
	};

	std::ifstream& operator>>(std::ifstream& is, std::vector <BoundaryCondition>& boundaries);
}
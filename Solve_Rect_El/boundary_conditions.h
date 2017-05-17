#pragma once
#include <fstream>
#include <vector>
#include "partition.h"
#include "parameters.h"
#include "myvector.h"
#include "integration.h"
#include "basis.h"
#include "EdgeLocation.h"
#include <string>

namespace boundary_conditions
{
	class BoundaryCondition
	{
	public:

		int elem;
		int edges[4]; //�����,������,������, �������: 1 - ����, 0 - ���
		int formula_number[4];
	};

	class BoundaryConditionsSupport : virtual public partition::Partition,
									  virtual public parameters::Parameters,
									  virtual public integration::Gauss_integration,
									  virtual public basis::Basis
	{
		double gamma;

		void calculate_boundaries1_for_side(int number, myvector::MyVector& b, EdgeSide edgeSide);
		void calculate_boundaries1(int number, myvector::MyVector& b);

	protected:

		std::vector <BoundaryCondition> boundaries1; //������ ������� �������
		std::vector <BoundaryCondition> boundaries2; //������ -//-
		std::vector <BoundaryCondition> boundaries3; //������ -//-

		void initialize_penalty_parameters(std::string fileName);
		void calculate_all_boundaries1(myvector::MyVector& b);
	};

	std::ifstream& operator>>(std::ifstream& is, std::vector <BoundaryCondition>& boundaries);
}
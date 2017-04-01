#pragma once
#include <fstream>
#include <vector>
#include "partition.h"
#include "parameters.h"
#include "myvector.h"
#include "integration.h"
#include "basis.h"
#include "EdgeLocation.h"

namespace boundary_conditions
{
	//enum Side {left_edge, right_edge, low_edge, up_edge};

	class BoundaryCondition
	{
	public:

		int elem;
		int edges[4]; //левое,правое,нижнее, верхнее: 1 - есть, 0 - нет
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

		std::vector <BoundaryCondition> boundaries1; //первые краевые условия
		std::vector <BoundaryCondition> boundaries2; //вторые -//-
		std::vector <BoundaryCondition> boundaries3; //третьи -//-

		void initialize_penalty_parameters();
		void calculate_all_boundaries1(myvector::MyVector& b);
	};

	std::ifstream& operator>>(std::ifstream& is, std::vector <BoundaryCondition>& boundaries);
}
#pragma once
#include <fstream>
#include <vector>
#include "partition.h"

namespace boundary_conditions
{
	class BoundaryCondition
	{
	public:

		int elem;
		int edges[4]; //�����,������,������, �������: 1 - ����, 0 - ���
		int formula_number;
	};

	class BoundaryConditionsSupport : public partition::Partition
	{
	public:

		std::vector <BoundaryCondition> boundaries1; //������ ������� �������
		std::vector <BoundaryCondition> boundaries2; //������ -//-
		std::vector <BoundaryCondition> boundaries3; //������ -//-

		void calculate_all_boundaries1();
		void calculate_boundaries1(int number);

		void calculate_boundaries1_left(int number);
		void calculate_boundaries1_right(int number);
		void calculate_boundaries1_low(int number);
		void calculate_boundaries1_up(int number);
	};

	std::ifstream& operator>>(std::ifstream& is, std::vector <BoundaryCondition>& boundaries);
}
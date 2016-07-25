#pragma once
#include <fstream>
#include <vector>
#include "partition.h"

using namespace std;

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

		vector <BoundaryCondition> boundaries1; //������ ������� �������
		vector <BoundaryCondition> boundaries2; //������ -//-
		vector <BoundaryCondition> boundaries3; //������ -//-

		void input_boundaries1(FILE* f_in);
		void input_boundaries2(FILE* f_in);
		void input_boundaries3(FILE* f_in);

		void calculate_all_boundaries1();
		void calculate_boundaries1(int number);

		void calculate_boundaries1_left(int number);
		void calculate_boundaries1_right(int number);
		void calculate_boundaries1_low(int number);
		void calculate_boundaries1_up(int number);
	};

	std::ifstream& operator>>(std::ifstream& is, vector <BoundaryCondition>& boundaries);
}
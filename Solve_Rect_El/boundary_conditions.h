#pragma once
#include <fstream>

namespace boundary_conditions
{
	class BoundaryCondition
	{
	public:

		int elem;
		int edges[4]; //левое,правое,нижнее, верхнее: 1 - есть, 0 - нет
		int formula_number;
	};

	std::ifstream& operator>>(std::ifstream& is, vector <BoundaryCondition>& boundaries)
	{
		int count;
		BoundaryCondition tmp;

		is >> count;
		boundaries.reserve(count);

		for(int i = 1; i <= count; i++)
		{
			is >> tmp.elem;
			is >> tmp.formula_number;
			is >> tmp.edges[0];
			is >> tmp.edges[1];
			is >> tmp.edges[2];
			is >> tmp.edges[3];
			boundaries.push_back(tmp);
		}

		return is;
	}
}
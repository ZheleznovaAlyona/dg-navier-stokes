#pragma once
#include <fstream>
#include <vector>

namespace element
{
	const int n_func_u_elem = 4;
	const int n_func_p_elem = 4;

	class Element
	{
	public:
		int nodes[4];
		int edges[4];
		int dof_u[n_func_u_elem];
		int ndof_u;
		int dof_p[n_func_p_elem];
		int ndof_p;
		int number_of_area;
		int neighbors[4]; //левый, правый, нижний, верхний

		Element& operator=(Element element);
	};

	std::ifstream& operator>>(std::ifstream& is, std::vector <Element>& elements);
}
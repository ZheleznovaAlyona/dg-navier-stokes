#include "element.h"
#include "basis.h"

using namespace std;
using namespace basis;

namespace element
{
	Element& Element::operator=(Element element)
	{
		ndof_u = element.ndof_u;
		ndof_p = element.ndof_p;
		for(int i = 0; i < 4; i++)
			nodes[i] = element.nodes[i];
		for(int i = 0; i < 4; i++)
			edges[i] = element.edges[i];
		for (int i = 0; i < ndof_u; i++)
			dof_u[i] = element.dof_u[i];
		for (int i = 0; i < ndof_p; i++)
			dof_p[i] = element.dof_p[i];
		number_of_area = element.number_of_area;
		for(int i = 0; i < 4; i++)
			neighbors[i] = element.neighbors[i];

		return *this;
	}

	ifstream& operator>>(ifstream& is, vector <Element>& elements)
	{
		int tmp;
		Element element_tmp;

		is >> tmp;
		elements.reserve(tmp);

		for(int i = 0; i < tmp; i++)
		{
			is >> element_tmp.number_of_area;
			for(int j = 0; j < 4; j++)
				is >> element_tmp.nodes[j];
			for(int j = 0; j < 4; j++)
				is >> element_tmp.neighbors[j];
			for(int j = 0; j < 4; j++)
				is >> element_tmp.edges[j];

			element_tmp.ndof_u = n_func_u;
			element_tmp.ndof_p = n_func_p;

			element_tmp.dof_u[0] = element_tmp.ndof_u * i;
			for (int j = 1; j < element_tmp.ndof_u; j++)
				element_tmp.dof_u[j] = element_tmp.dof_u[0] + j;

			element_tmp.dof_p[0] = tmp * n_func_u + element_tmp.ndof_p * i;
			for (int j = 1; j < element_tmp.ndof_p; j++)
				element_tmp.dof_p[j] = element_tmp.dof_p[0] + j;

			elements.push_back(element_tmp);
		}

		return is;
	}
}
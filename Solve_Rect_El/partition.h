#pragma once
#include <fstream>
#include <vector>
#include "element.h"
#include "point.h"


namespace partition
{
	class Partition
	{
	public:

		std::vector <element::Element> elements;
		std::vector <point::Point> nodes;
		void input(std::ifstream& grid_f_in, std::ifstream& elements_f_in);
		double get_hx(int element_number);
		double get_hy(int element_number);
		int count_unzero_matrix_elements();
		int create_unzero_elements_list(int element_number, 
										std::vector <int> &list, 
										int dof_num_i, 
										int dof_num_j, 
										int *dof_i, 
										int *dof_j,
										bool dof_j_edge);
		int search_element(double x, double y);
	};
}
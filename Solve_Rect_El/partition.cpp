#include "partition.h"

using namespace std;
using namespace element;
using namespace point;

namespace partition
{
	void Partition::input(ifstream& grid_f_in, ifstream& elements_f_in)
	{
		grid_f_in >> nodes;
		elements_f_in >> elements;
	};

	double Partition::get_hx(int element_number)
	{
		Element element = elements[element_number];
		return nodes[element.nodes[1]].x - nodes[element.nodes[0]].x;
	}

	double Partition::get_hy(int element_number)
	{
		Element element = elements[element_number];
		return nodes[element.nodes[2]].y - nodes[element.nodes[0]].y;
	}

	int Partition::count_unzero_matrix_elements()
	{
		int count_uu = 0;
		for(unsigned int i = 0; i < elements.size(); i++)
		{
			//с собой
			count_uu += 4;
			//с соседями
			for(int j = 0; j < 4; j++)
				if(elements[i].neighbors[j] != -1)
					count_uu += 4;
		}
		count_uu *= 4;
		count_uu -= elements.size() * 4;
		//так как нужно для одного треугольника ввиду симметричности портрета,
		//то необходимо полученное количество поделить на 2
		count_uu /= 2;

		int count_pp = 0;
		for(unsigned int i = 0; i < elements.size(); i++)
		{
			//с собой
			count_pp += 4;
			//с соседями
			for(int j = 0; j < 4; j++)
				if(elements[i].neighbors[j] != -1)
					count_pp += 4;
		}
		count_pp *= 4;
		count_pp -= nodes.size();
		//так как нужно для одного треугольника ввиду симметричности портрета,
		//то необходимо полученное количество поделить на 2
		count_pp /= 2;

		int count_up = 0;
		for(unsigned int i = 0; i < elements.size(); i++)
		{
			//с собой
			count_up += 4;
			//с соседями
			for(int j = 0; j < 4; j++)
				if(elements[i].neighbors[j] != -1)
					count_up += 4;
		}
		count_up *= 4;

		return count_uu + count_pp + count_up;
	}

	int Partition::create_unzero_elements_list(int element_number, 
									vector <int> &list, 
									int dof_num_i, 
									int dof_num_j, 
									int *dof_i, 
									int *dof_j,
									bool dof_j_edge)
	{
		int neighbor;

		//свои по строкам
		for(int i = 0; i < dof_num_i; i++)
			list.push_back(dof_i[i]);

		//свои
		for(int i = 0; i < dof_num_j; i++)
			list.push_back(dof_j[i]);

		//соседей
		for(int j = 0; j < 4; j++)
		{
			neighbor = elements[element_number].neighbors[j];
			if(neighbor != -1)
			for(int i = 0; i < dof_num_j; i++)
			{
				if(dof_j_edge) list.push_back(elements[neighbor].edges[i]);
				else list.push_back(elements[neighbor].nodes[i]);
			}			
		}

		return list.size();
	}

	int Partition::search_element(double x, double y)
	{
		double x_left, x_right, y_low, y_up;
		int size = elements.size();

		for(int i = 0; i < size; i++)
		{
			x_left = nodes[elements[i].nodes[0]].x;
			x_right = nodes[elements[i].nodes[1]].x;
			y_low = nodes[elements[i].nodes[0]].y;
			y_up = nodes[elements[i].nodes[3]].y;
			if(x_left <= x && x <= x_right && y_low <= y && y <= y_up)
				return i;
		}
	}

}
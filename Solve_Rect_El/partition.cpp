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
}
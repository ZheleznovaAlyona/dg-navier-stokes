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
	};
}
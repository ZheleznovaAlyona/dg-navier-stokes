#pragma once
#include <fstream>
#include "element.h"
#include "point.h"


using namespace std;
using namespace element;
using namespace point;

namespace partition
{
	class Partition
	{
	public:

		vector <Element> elements;
		vector <Point> nodes;
		void input(ifstream& grid_f_in, ifstream& elements_f_in)
		{
			grid_f_in >> nodes;
			elements_f_in >> elements;
		};
	};

}
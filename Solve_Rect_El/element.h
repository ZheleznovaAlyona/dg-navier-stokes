#pragma once
#include <fstream>
#include <vector>

using namespace std;

namespace element
{
	class Element
	{
	public:
		int nodes[4];
		int edges[4];
		int number_of_area;
		int neighbors[4]; //левый, правый, нижний, верхний

		Element& operator=(Element element);
	};

	std::ifstream& operator>>(std::ifstream& is, vector <Element>& elements);
}
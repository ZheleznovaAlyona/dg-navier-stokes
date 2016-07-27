#include "element.h"

using namespace std;

namespace element
{
	Element& Element::operator=(Element element)
	{
		for(int i = 0; i < 4; i++)
			nodes[i] = element.nodes[i];
		for(int i = 0; i < 4; i++)
			edges[i] = element.edges[i];
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

			elements.push_back(element_tmp);
		}

		return is;
	}
}
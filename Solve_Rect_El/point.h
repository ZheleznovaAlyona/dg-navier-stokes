#pragma once
#include <fstream>
#include <vector>

namespace point
{
	class Point
	{
	public:

		double x;
		double y;
		bool operator==(Point point);
	};

	std::ifstream& operator>>(std::ifstream& is, std::vector <Point>& points);
}
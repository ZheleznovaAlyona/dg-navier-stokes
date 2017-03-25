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
		Point();
		Point(double xx, double yy);
		~Point();
		bool operator==(Point point);
		Point operator*(double a);
		Point operator+(Point pt);
		double norm();
	};

	std::ifstream& operator>>(std::ifstream& is, std::vector <Point>& points);
	double pt_scal(Point p1, Point p2);
}
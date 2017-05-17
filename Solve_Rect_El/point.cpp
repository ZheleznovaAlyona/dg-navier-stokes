#include "point.h"

using namespace std;

namespace point
{
	Point::Point()
	{}

	Point::Point(double xx, double yy)
	{
		x = xx;
		y = yy;
	}

	Point::~Point()
	{}

	bool Point::operator==(Point& point)
	{
		if (point.x == x && point.y == y)
			return true;
		else
			return false;
	}

	Point Point::operator*(double a)
	{
		return Point(a * x, a * y);
	}

	Point Point::operator+(Point& pt)
	{
		return Point(x + pt.x, y + pt.y);
	}

	double Point::norm()
	{
		return sqrt(x * x + y * y);
	}

	ifstream& operator>>(ifstream& is, vector <Point>& points)
	{
		int tmp;
		Point point_tmp;

		is >> tmp;

		points.reserve(tmp);
		for (int i = 0; i < tmp; i++)
		{
			is >> point_tmp.x >> point_tmp.y;
			points.push_back(point_tmp);
		}

		return is;
	}

	double pt_scal(Point& p1, Point& p2)
	{
		return p1.x * p2.x + p1.y * p2.y;
	}
}
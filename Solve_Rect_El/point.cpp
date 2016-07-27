#include "point.h"

using namespace std;

namespace point
{
	bool Point::operator==(Point point)
	{
		if(point.x == x && point.y == y)
			return true;
		else
			return false;
	}

	ifstream& operator>>(ifstream& is, vector <Point>& points)
	{
		int tmp;
		Point point_tmp;

		is >> tmp;

		points.reserve(tmp);
		for(int i = 0; i < tmp; i++)
		{
			is >> point_tmp.x >> point_tmp.y;
			points.push_back(point_tmp);
		}

		return is;
	}
}
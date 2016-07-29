#include "integration.h"
#include <math.h>

namespace integration
{
	void Gauss_integration::initialize()
	{
		double tmp_gauss_points[2][9] =
		{
			{0.0, 0.0, 0.0, sqrt(3.0 / 5.0), -sqrt(3.0 / 5.0), sqrt(3.0 / 5.0), sqrt(3.0 / 5.0), -sqrt(3.0 / 5.0), -sqrt(3.0 / 5.0)},
			{0.0, sqrt(3.0 / 5.0), -sqrt(3.0 / 5.0), 0.0, 0.0, sqrt(3.0 / 5.0), -sqrt(3.0 / 5.0), sqrt(3.0 / 5.0), -sqrt(3.0 / 5.0)}
		};
		for(int i = 0; i < 2; i++)
			for(int j = 0; j < 9; j++)
				gauss_points[i][j] = tmp_gauss_points[i][j];

		gauss_weights[0] = 64.0 / 81.0;
		gauss_weights[1] = gauss_weights[2] = gauss_weights[3] = gauss_weights[4] = 40.0 / 81.0;
		gauss_weights[5] = gauss_weights[6] = gauss_weights[7] = gauss_weights[8] = 25.0 / 81.0;

		gauss_points_1[0] = -sqrt(3.0 / 5.0); gauss_points_1[1] = 0.0; gauss_points_1[2] = sqrt(3.0 / 5.0);
		gauss_weights_1[0] = 5.0 / 9.0; gauss_weights_1[1] = 8.0 / 9.0; gauss_weights_1[2] =  5.0 / 9.0;
	}
}
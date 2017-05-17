#include "integration.h"
#include <math.h>

namespace integration
{
	void Gauss_integration::initialize()
	{
		if (n_ip1D == 3)
		{
			gauss_points_1[0] = -sqrt(3.0 / 5.0);
			gauss_points_1[1] = 0.0;
			gauss_points_1[2] = -gauss_points_1[0];

			gauss_weights_1[0] = 5.0 / 9.0;
			gauss_weights_1[1] = 8.0 / 9.0;
			gauss_weights_1[2] = gauss_weights_1[0];
		}

		if (n_ip1D == 5)
		{
			gauss_points_1[0] = 0;
			gauss_points_1[1] = (1.0 / 3.0) * sqrt(5.0 - 2.0 * sqrt(10.0 / 7.0));
			gauss_points_1[2] = -gauss_points_1[1];
			gauss_points_1[3] = (1.0 / 3.0) * sqrt(5.0 + 2.0 * sqrt(10.0 / 7.0));
			gauss_points_1[4] = -gauss_points_1[3];

			gauss_weights_1[0] = 128.0 / 225.0;
			gauss_weights_1[1] = (322.0 + 13.0 * sqrt(70.0)) / 900.0;
			gauss_weights_1[2] = gauss_weights_1[1];
			gauss_weights_1[3] = (322.0 - 13.0 * sqrt(70.0)) / 900.0;
			gauss_weights_1[4] = gauss_weights_1[3];
		}

		for (size_t i = 0; i < n_ip1D; i++)
		{
			for (size_t j = 0; j < n_ip1D; j++)
			{
				gauss_points[0][i * n_ip1D + j] = gauss_points_1[j];
				gauss_points[1][i * n_ip1D + j] = gauss_points_1[i];
				gauss_weights[i * n_ip1D + j] = gauss_weights_1[i] * gauss_weights_1[j];
			}
		}
	}
}
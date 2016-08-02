#pragma once
namespace integration
{
	class Gauss_integration
	{
	public:
		double gauss_points[2][9];//точки гаусса
		double gauss_weights[9];// веса гаусса
		double gauss_points_1[3];//точки гаусса
		double gauss_weights_1[3];// веса гаусса

		void initialize();
	};

}
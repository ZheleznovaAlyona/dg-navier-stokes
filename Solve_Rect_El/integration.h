#pragma once
namespace integration
{
	//const int n_ip = 9;
	//const int n_ip1D = 3;
	const int n_ip = 25;
	const int n_ip1D = 5;

	class Gauss_integration
	{
	public:
		double gauss_points[2][n_ip];//точки гаусса
		double gauss_weights[n_ip];// веса гаусса
		double gauss_points_1[n_ip1D];//точки гаусса
		double gauss_weights_1[n_ip1D];// веса гаусса

		void initialize();
	};

}
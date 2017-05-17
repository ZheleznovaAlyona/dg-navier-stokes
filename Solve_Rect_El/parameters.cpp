#include "parameters.h"
#include "testing_parameters.h"
#include "myfunctions.h"

using namespace testingparameters;

namespace parameters
{
	double Parameters::gx(int formula_number, double x, double y)
	{
		if(Testing_parameters::test == 1)
			switch(formula_number)
			{
				case 0: return x;
				case 1:	return x;
				default: return 1.0;
			}

		if(Testing_parameters::test == 3)
			switch(formula_number)
			{
				case 0: return 20 * x * y * y * y;
				case 1:	return 20 * x * y * y * y;
				default: return 1.0;
			}

		return 1.0;
	}

	double Parameters::gy(int formula_number, double x, double y)
	{
		if(Testing_parameters::test == 1)
			switch(formula_number)
			{
				case 0: return -y;
				case 1:	return -y;
				default: return 1.0;
			}

		if(Testing_parameters::test == 3)
			switch(formula_number)
			{
				case 0: return 5 * x * x * x * x - 5 * y * y * y * y;
				case 1:	return 5 * x * x * x * x - 5 * y * y * y * y;
				default: return 1.0;
			}

		return 1.0;
	}

	double Parameters::calculate_fx(int area_number, double x, double y)
	{
		if(Testing_parameters::test == 1)
			switch(area_number)
			{
				case 0: return 1.0 / calculate_rho(area_number) * 0;
				case 1:	return 1.0 / calculate_rho(area_number) * 0;
				default: return 1.0;
			}

		if(Testing_parameters::test == 3)
			switch(area_number)
			{
			case 0: return -calculate_lambda(area_number) * 120 * x * y + 
							1.0 / calculate_rho(area_number) * 120 * x * y + 
							400 * x * pow_i(6, y) + 300 * (pow_i(5, x) * y * y
							- x * pow_i(6, y));
			case 1:	return -calculate_lambda(area_number) * 120 * x * y + 
							1.0 / calculate_rho(area_number) * 120 * x * y + 
							400 * x * pow_i(6, y) + 300 * (pow_i(5, x) * y * y
							- x * pow_i(6, y));
			default: return 1.0;
			}

		return 1.0;
	}

	double Parameters::calculate_fy(int area_number, double x, double y)
	{
		if(Testing_parameters::test == 1)
			switch(area_number)
			{
				case 0: return 1.0 / calculate_rho(area_number) * 0;
				case 1:	return 1.0 / calculate_rho(area_number) * 0;
				default: return 1.0;
			}

		if(Testing_parameters::test == 3)
			switch(area_number)
			{
			case 0: return -calculate_lambda(area_number) * (60 * x * x - 60 * y * y) + 
							1.0 / calculate_rho(area_number) * (60 * x * x - 60 * y * y) + 
							400 * pow_i(4, x)  * pow_i(3, y) + 100 * (pow_i(7, y) 
							- pow_i(4, x)  * pow_i(3, y));
			case 1:	return -calculate_lambda(area_number) * (60 * x * x - 60 * y * y) + 
							1.0 / calculate_rho(area_number) * (60 * x * x - 60 * y * y) + 
							400 * pow_i(4, x)  * pow_i(3, y) + 100 * (pow_i(7, y)
							- pow_i(4, x)  * pow_i(3, y));
			default: return 1.0;
			}

		return 1.0;
	}

	double Parameters::calculate_ux_analytic(int area_number, double x, double y)
	{
		if(Testing_parameters::test == 1)
			switch(area_number)
			{
				case 0: return x;
				case 1:	return x;
				default: return 1.0;
			}

		if(Testing_parameters::test == 3)
			switch(area_number)
			{
				case 0: return 20 * x * y * y * y;
				case 1:	return 20 * x * y * y * y;
				default: return 1.0;
			}

		return 1.0;
	}

	double Parameters::calculate_uy_analytic(int area_number, double x, double y)
	{
		if(Testing_parameters::test == 1)
			switch(area_number)
			{
				case 0: return -y;
				case 1:	return -y;
				default: return 1.0;
			}

		if(Testing_parameters::test == 3)
			switch(area_number)
			{
				case 0: return 5 * x * x * x * x - 5 * y * y * y * y;
				case 1:	return 5 * x * x * x * x - 5 * y * y * y * y;
				default: return 1.0;
			}

		return 1.0;
	}

	double Parameters::calculate_uxdx_analytic(int area_number, double x, double y)
	{
		if(Testing_parameters::test == 1)
			switch(area_number)
			{
				case 0: return 1;
				case 1:	return 1;
				default: return 1.0;
			}

		if(Testing_parameters::test == 3)
			switch(area_number)
			{
				case 0: return 20.0 * y * y * y;
				case 1:	return 20.0 * y * y * y;
				default: return 1.0;
			}

		return 1.0;
	}

	double Parameters::calculate_uydy_analytic(int area_number, double x, double y)
	{
		if(Testing_parameters::test == 1)
			switch(area_number)
			{
				case 0: return -1;
				case 1:	return -1;
				default: return 1.0;
			}
		if(Testing_parameters::test == 3)
			switch(area_number)
			{
				case 0: return -20.0 * y * y * y;
				case 1:	return -20.0 * y * y * y;
				default: return 1.0;
			}

		return 1.0;
	}

	double Parameters::calculate_p_analytic(int area_number, double x, double y)
	{
		if(Testing_parameters::test == 1)
			switch(area_number)
			{
				case 0: return 0;
				case 1:	return 0;
				default: return 1.0;
			}
		if(Testing_parameters::test == 3)
			switch(area_number)
			{
				case 0: return 60 * x * x * y - 20 * y * y * y - 5;
				case 1:	return 60 * x * x * y - 20 * y * y * y - 5;
				default: return 1.0;
			}

		return 1.0;
	}

	double Parameters::calculate_lambda(int area_number)
	{
		if(Testing_parameters::test == 1)
			switch(area_number)
			{
				case 0: return 1.0;
				case 1:	return 1.0;
				default: return 1.0;
			}

		if(Testing_parameters::test == 3)
			switch(area_number)
			{
				case 0: return 1.0;
				case 1:	return 1.0;
				default: return 1.0;
			}

		return 1.0;
	}

	double Parameters::calculate_rho(int area_number)
	{
		if(Testing_parameters::test == 1)
			switch(area_number)
			{
				case 0: return 1.0;
				case 1:	return 1.0;
				default: return 1.0;
			}
		if(Testing_parameters::test == 3)
			switch(area_number)
			{
				case 0: return 1.0;
				case 1:	return 1.0;
				default: return 1.0;
			}

		return 1.0;
	}
}
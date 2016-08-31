#include "basis.h"

using namespace std;
using namespace partition;

namespace basis
{
	int Basis::get_n_func_u()
	{
		return n_func_u;
	}

	int Basis::get_n_func_p()
	{
		return n_func_p;
	}

	void Basis::initialize()
	{
		n_func_p = 4;
		n_func_u = 4;

		phix[0] = [](double ksi, double etta) { return 0.5 * (1 - ksi); };
		phix[1] = [](double ksi, double etta) { return 0.5 * (1 + ksi); };
		phix[2] = [](double ksi, double etta) { return 0.0; };
		phix[3] = [](double ksi, double etta) { return 0.0; };

		phiy[0] = [](double ksi, double etta) { return 0.0; };
		phiy[1] = [](double ksi, double etta) { return 0.0; };
		phiy[2] = [](double ksi, double etta) { return 0.5 * (1 - etta); };
		phiy[3] = [](double ksi, double etta) { return 0.5 * (1 + etta); };

		dphixksi[0] = [](double ksi, double etta) { return -0.5; };
		dphixksi[1] = [](double ksi, double etta) { return 0.5; };
		dphixksi[2] = [](double ksi, double etta) { return 0.0; };
		dphixksi[3] = [](double ksi, double etta) { return  0.0; };

		dphiyksi[0] = [](double ksi, double etta) { return 0.0; };
		dphiyksi[1] = [](double ksi, double etta) { return 0.0; };
		dphiyksi[2] = [](double ksi, double etta) { return 0.0; };
		dphiyksi[3] = [](double ksi, double etta) { return 0.0; };

		dphixetta[0] = [](double ksi, double etta) { return 0.0; };
		dphixetta[1] = [](double ksi, double etta) { return 0.0; };
		dphixetta[2] = [](double ksi, double etta) { return 0.0; };
		dphixetta[3] = [](double ksi, double etta) { return 0.0; };

		dphiyetta[0] = [](double ksi, double etta) { return 0.0; };
		dphiyetta[1] = [](double ksi, double etta) { return 0.0; };
		dphiyetta[2] = [](double ksi, double etta) { return -0.5; };
		dphiyetta[3] = [](double ksi, double etta) { return 0.5; };

		psi[0] = [](double ksi, double etta) { return 0.25 * (1 - ksi) * (1 - etta); };
		psi[1] = [](double ksi, double etta) { return 0.25 * (1 + ksi) * (1 - etta); };
		psi[2] = [](double ksi, double etta) { return 0.25 * (1 - ksi) * (1 + etta); };
		psi[3] = [](double ksi, double etta) { return 0.25 * (1 + ksi) * (1 + etta); };

		dpsiksi[0] = [](double ksi, double etta) { return -0.25 * (1 - etta); };
		dpsiksi[1] = [](double ksi, double etta) { return 0.25 * (1 - etta); };
		dpsiksi[2] = [](double ksi, double etta) { return -0.25 * (1 + etta); };
		dpsiksi[3] = [](double ksi, double etta) { return 0.25 * (1 + etta); };

		dpsietta[0] = [](double ksi, double etta) { return -0.25 * (1 - ksi); };
		dpsietta[1] = [](double ksi, double etta) { return -0.25 * (1 + ksi); };
		dpsietta[2] = [](double ksi, double etta) { return 0.25 * (1 - ksi); };
		dpsietta[3] = [](double ksi, double etta) { return 0.25 * (1 + ksi); };
	}

	double Basis::phix_i(int i, double x, double y, int element_number, Partition& p)
	{
		double x_left = p.nodes[p.elements[element_number].nodes[0]].x;
		double x_right = p.nodes[p.elements[element_number].nodes[1]].x;
		double y_low = p.nodes[p.elements[element_number].nodes[0]].y;
		double y_up = p.nodes[p.elements[element_number].nodes[3]].y;
		double hx = x_right - x_left, hy = y_up - y_low;
		double ksi = 2 * (x - (x_left + x_right) / 2) / hx, 
			etta = 2 * (y - (y_low + y_up) / 2) / hy;
		return phix[i](ksi, etta);
	}

	double Basis::phiy_i(int i, double x, double y, int element_number, Partition& p)
	{
		double x_left = p.nodes[p.elements[element_number].nodes[0]].x;
		double x_right = p.nodes[p.elements[element_number].nodes[1]].x;
		double y_low = p.nodes[p.elements[element_number].nodes[0]].y;
		double y_up = p.nodes[p.elements[element_number].nodes[3]].y;
		double hx = x_right - x_left, hy = y_up - y_low;
		double ksi = 2 * (x - (x_left + x_right) / 2) / hx, 
			etta = 2 * (y - (y_low + y_up) / 2) / hy;
		return phiy[i](ksi, etta);
	}

	double Basis::phixdx_i(int i, double x, double y, int element_number, Partition& p)
	{
		double x_left = p.nodes[p.elements[element_number].nodes[0]].x;
		double x_right = p.nodes[p.elements[element_number].nodes[1]].x;
		double y_low = p.nodes[p.elements[element_number].nodes[0]].y;
		double y_up = p.nodes[p.elements[element_number].nodes[3]].y;
		double hx = x_right - x_left, hy = y_up - y_low;
		double ksi = 2 * (x - (x_left + x_right) / 2) / hx, 
			etta = 2 * (y - (y_low + y_up) / 2) / hy;
		return dphixksi[i](ksi, etta) / hx;
	}

	double Basis::phiydy_i(int i, double x, double y, int element_number, Partition& p)
	{
		double x_left = p.nodes[p.elements[element_number].nodes[0]].x;
		double x_right = p.nodes[p.elements[element_number].nodes[1]].x;
		double y_low = p.nodes[p.elements[element_number].nodes[0]].y;
		double y_up = p.nodes[p.elements[element_number].nodes[3]].y;
		double hx = x_right - x_left, hy = y_up - y_low;
		double ksi = 2 * (x - (x_left + x_right) / 2) / hx, 
			etta = 2 * (y - (y_low + y_up) / 2) / hy;
		return dphiyetta[i](ksi, etta) / hy;
	}

	double Basis::psi_i(int i, double x, double y, int element_number, Partition& p)
	{
		double x_left = p.nodes[p.elements[element_number].nodes[0]].x;
		double x_right = p.nodes[p.elements[element_number].nodes[1]].x;
		double y_low = p.nodes[p.elements[element_number].nodes[0]].y;
		double y_up = p.nodes[p.elements[element_number].nodes[3]].y;
		double hx = x_right - x_left, hy = y_up - y_low;
		double ksi = 2 * (x - (x_left + x_right) / 2) / hx, 
			etta = 2 * (y - (y_low + y_up) / 2) / hy;
		return psi[i](ksi, etta);
	}
}
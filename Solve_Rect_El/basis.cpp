#include "basis.h"

using namespace std;
using namespace partition;

namespace basis
{
	void Basis::initialize()
	{
		phix[0] = [](double ksi, double etta) { return 0.5 * (1 - ksi); };
		phix[1] = [](double ksi, double etta) { return 0.5 * (1 + ksi); };
		phix[2] = [](double ksi, double etta) { return 0.0; };
		phix[3] = [](double ksi, double etta) { return 0.0; };
		phix[4] = [](double ksi, double etta) { return 0.5 * (1 - ksi) * etta; };
		phix[5] = [](double ksi, double etta) { return 0.5 * (1 + ksi) * etta; };
		phix[6] = [](double ksi, double etta) { return 0.0; };
		phix[7] = [](double ksi, double etta) { return 0.0; };
		phix[8] = [](double ksi, double etta) { return 1 - ksi * ksi; };
		phix[9] = [](double ksi, double etta) { return (1 - ksi * ksi) * etta; };
		phix[10] = [](double ksi, double etta) { return 0.0; };
		phix[11] = [](double ksi, double etta) { return 0.0; };

		phiy[0] = [](double ksi, double etta) { return 0.0; };
		phiy[1] = [](double ksi, double etta) { return 0.0; };
		phiy[2] = [](double ksi, double etta) { return 0.5 * (1 - etta); };
		phiy[3] = [](double ksi, double etta) { return 0.5 * (1 + etta); };
		phiy[4] = [](double ksi, double etta) { return 0.0; };
		phiy[5] = [](double ksi, double etta) { return 0.0; };
		phiy[6] = [](double ksi, double etta) { return 0.5 * (1 - etta) * ksi; };
		phiy[7] = [](double ksi, double etta) { return 0.5 * (1 + etta) * ksi; };
		phiy[8] = [](double ksi, double etta) { return 0.0; };
		phiy[9] = [](double ksi, double etta) { return 0.0; };
		phiy[10] = [](double ksi, double etta) { return 1 - etta * etta; };
		phiy[11] = [](double ksi, double etta) { return (1 - etta * etta) * ksi; };

		dphixksi[0] = [](double ksi, double etta) { return -0.5; };
		dphixksi[1] = [](double ksi, double etta) { return 0.5; };
		dphixksi[2] = [](double ksi, double etta) { return 0.0; };
		dphixksi[3] = [](double ksi, double etta) { return  0.0; };
		dphixksi[4] = [](double ksi, double etta) { return -0.5 * etta; };
		dphixksi[5] = [](double ksi, double etta) { return 0.5 * etta; };
		dphixksi[6] = [](double ksi, double etta) { return 0.0; };
		dphixksi[7] = [](double ksi, double etta) { return  0.0; };
		dphixksi[8] = [](double ksi, double etta) { return -2 * ksi; };
		dphixksi[9] = [](double ksi, double etta) { return -2 * ksi * etta; };
		dphixksi[10] = [](double ksi, double etta) { return 0.0; };
		dphixksi[11] = [](double ksi, double etta) { return  0.0; };

		dphiyksi[0] = [](double ksi, double etta) { return 0.0; };
		dphiyksi[1] = [](double ksi, double etta) { return 0.0; };
		dphiyksi[2] = [](double ksi, double etta) { return 0.0; };
		dphiyksi[3] = [](double ksi, double etta) { return 0.0; };
		dphiyksi[4] = [](double ksi, double etta) { return 0.0; };
		dphiyksi[5] = [](double ksi, double etta) { return 0.0; };
		dphiyksi[6] = [](double ksi, double etta) { return 0.5 * (1 - etta); };
		dphiyksi[7] = [](double ksi, double etta) { return  0.5 * (1 + etta); };
		dphiyksi[8] = [](double ksi, double etta) { return 0.0; };
		dphiyksi[9] = [](double ksi, double etta) { return 0.0; };
		dphiyksi[10] = [](double ksi, double etta) { return 0.0; };
		dphiyksi[11] = [](double ksi, double etta) { return  1 - etta * etta; };

		dphixetta[0] = [](double ksi, double etta) { return 0.0; };
		dphixetta[1] = [](double ksi, double etta) { return 0.0; };
		dphixetta[2] = [](double ksi, double etta) { return 0.0; };
		dphixetta[3] = [](double ksi, double etta) { return 0.0; };
		dphixetta[4] = [](double ksi, double etta) { return 0.5 * (1 - ksi); };
		dphixetta[5] = [](double ksi, double etta) { return 0.5 * (1 + ksi); };
		dphixetta[6] = [](double ksi, double etta) { return 0.0; };
		dphixetta[7] = [](double ksi, double etta) { return  0.0; };
		dphixetta[8] = [](double ksi, double etta) { return 0.0; };
		dphixetta[9] = [](double ksi, double etta) { return 1 - ksi * ksi; };
		dphixetta[10] = [](double ksi, double etta) { return 0.0; };
		dphixetta[11] = [](double ksi, double etta) { return  0.0; };

		dphiyetta[0] = [](double ksi, double etta) { return 0.0; };
		dphiyetta[1] = [](double ksi, double etta) { return 0.0; };
		dphiyetta[2] = [](double ksi, double etta) { return -0.5; };
		dphiyetta[3] = [](double ksi, double etta) { return 0.5; };
		dphiyetta[4] = [](double ksi, double etta) { return 0.0; };
		dphiyetta[5] = [](double ksi, double etta) { return 0.0; };
		dphiyetta[6] = [](double ksi, double etta) { return -0.5 * ksi; };
		dphiyetta[7] = [](double ksi, double etta) { return  0.5 * ksi; };
		dphiyetta[8] = [](double ksi, double etta) { return 0.0; };
		dphiyetta[9] = [](double ksi, double etta) { return 0.0; };
		dphiyetta[10] = [](double ksi, double etta) { return - 2 * etta; };
		dphiyetta[11] = [](double ksi, double etta) { return  - 2 * etta * ksi; };

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
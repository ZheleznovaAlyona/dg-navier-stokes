#include "basis.h"

using namespace std;

namespace basis
{
	void Basis::initialize()
	{
		n_func = 4;

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
}
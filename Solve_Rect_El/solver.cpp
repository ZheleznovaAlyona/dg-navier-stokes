#include "solver.h"

using namespace solverparameters;
using namespace std;
using namespace myvector;
using namespace slae;
using namespace densematrix;

namespace solver
{
	Solver::Solver(){}

	Solver::Solver(string& s_parameters_f)
	{
		s_parameters.initialize(s_parameters_f);
	}

	Solver::~Solver(){}
}
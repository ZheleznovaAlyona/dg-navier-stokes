#include "solver_parameters.h"

using namespace std;

#include <libconfig.h++>

using namespace libconfig;

namespace solverparameters
{
	void SolverParameters::initialize(std::string file_name)
	{
		Config cfg;
		cfg.readFile(file_name.c_str());
		Setting& root = cfg.getRoot();

		Setting& solver = root["solver"];

		if(!solver.lookupValue("epsilon", this->epsilon))
			this->epsilon = 1E-10;

		if(!solver["gmres"].lookupValue("m", this->gmres_m)) 
			this->gmres_m = 15;

		if(!solver.lookupValue("max_number_of_iterations", this->max_number_of_iterations))
			max_number_of_iterations = 1000;

		if(!solver.lookupValue("max_number_of_iterations_non_linear", this->max_number_of_iterations_non_lin))
			max_number_of_iterations_non_lin = 10;
	}
}
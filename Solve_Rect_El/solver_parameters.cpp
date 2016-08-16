#include "solver_parameters.h"

using namespace std;

// по необъяснимым причинам библиотека сама не линкуется
#if defined(_WIN32) && !defined(_WIN64)
#pragma comment(lib, "C:/Users/прол/Documents/проекты/NS/NS/DG/Solve_Rect_El/packages/libconfig_vc120.1.4.9.4/build/native/bin/libconfig-x86-v120-mt-1_4_9_4.imp.lib")
#else
#pragma comment(lib, "C:/Users/прол/Documents/проекты/NS/NS/DG/Solve_Rect_El/packages/libconfig_vc120.1.4.9.4/build/native/bin/libconfig-x64-v120-mt-1_4_9_4.imp.lib")
#endif

// https://habrahabr.ru/company/abbyy/blog/136909/
#pragma warning (push)
#pragma warning (disable: 4290)
#include <libconfig.h++>
#pragma warning (pop)

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
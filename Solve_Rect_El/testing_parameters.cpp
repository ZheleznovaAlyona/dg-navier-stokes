#include "testing_parameters.h"

namespace testingparameters
{

	bool Testing_parameters::use_LU = false;
	int Testing_parameters::test = 0;
	int Testing_parameters::solver = 0;

	Testing_parameters::Testing_parameters() {}

	Testing_parameters::Testing_parameters(bool use_LU_in, int test_in, int solver_in)
	{
		use_LU = use_LU_in;
		test = test_in;
		solver = solver_in;
	}

	Testing_parameters::~Testing_parameters(){}

	const Testing_parameters& Testing_parameters::instance()
    {  
            static Testing_parameters the_single_instance(false, 3, 2);
            return the_single_instance;
    }

	Testing_parameters::Testing_parameters(Testing_parameters& parameters){}

}
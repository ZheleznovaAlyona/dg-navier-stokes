#pragma once
#include <fstream>
#include <string>

namespace solverparameters
{
	class SolverParameters
	{
	public:		
		int max_number_of_iterations;
		int max_number_of_iterations_non_lin; 
		double epsilon;
		int gmres_m;

		void initialize(std::string file_name);
	};
}
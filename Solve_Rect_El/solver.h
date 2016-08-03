#pragma once
#include "myvector.h"
#include "densematrix.h"
#include <fstream>
#include "SLAE.h"

namespace solver
{
	class Solver
	{
	public:
		virtual myvector::MyVector Solve(myvector::MyVector U_begin,
										 double &normL2u,
										 double &normL2p,
										 slae::SLAE& slae_in) = 0;
		void si_print(std::ofstream& log_f, 
					  int iteration_number,
					  double &normL2u,
					  double &normL2p,
					  slae::SLAE& slae_in);
		double find_relaxation_parameter(myvector::MyVector q_current, 
										 myvector::MyVector q_previous, 
										 double &residual_previous,
										 slae::SLAE& slae_in);
		void simple_iterations(slae::SLAE& slae_in);
		void run(std::ofstream& solution_f_out,
				 std::ofstream& info_f_out,
				 slae::SLAE& slae_in);
	};

	class BCG : public Solver
	{
		myvector::MyVector Solve(myvector::MyVector U_begin, double &normL2u, double &normL2p, myvector::MyVector& solution);
	};

	class BiCGStab : public Solver
	{
		myvector::MyVector Solve(myvector::MyVector U_begin, double &normL2u, double &normL2p, myvector::MyVector& solution);
	};

	class GMRES : public Solver
	{
		myvector::MyVector Solve(myvector::MyVector U_begin, double &normL2u, double &normL2p, myvector::MyVector& solution);
		void solve_min_sqr_problem(myvector::MyVector d, densematrix::DenseMatrix H, myvector::MyVector &result);
	};
}
#pragma once
#include "myvector.h"
#include "densematrix.h"
#include <fstream>
#include "SLAE.h"
#include "solver_parameters.h"
#include "logger.h"

namespace solver
{
	class Solver
	{

		public:
			solverparameters::SolverParameters s_parameters;

			Solver();
			Solver(std::string& s_parameters_f);
			~Solver();

			virtual myvector::MyVector solve(myvector::MyVector U_begin,
											 double &normL2u,
											 double &normL2p,
											 slae::SLAE& slae_in,
											 logger::Logger& my_logger) = 0;
	};

	class BCG : virtual public Solver
	{
	public:
		myvector::MyVector solve(myvector::MyVector U_begin, double &normL2u, double &normL2p, slae::SLAE& slae_in, logger::Logger& my_logger);
	};

	class BiCGStab : public Solver
	{
	public:
		myvector::MyVector solve(myvector::MyVector U_begin, double &normL2u, double &normL2p, slae::SLAE& slae_in, logger::Logger& my_logger) final;
	};

	class GMRES : virtual public Solver
	{
	public:
		myvector::MyVector solve(myvector::MyVector U_begin, double &normL2u, double &normL2p, slae::SLAE& slae_in, logger::Logger& my_logger);
		void solve_min_sqr_problem(myvector::MyVector d, densematrix::DenseMatrix H, myvector::MyVector &result);
	};

	class BCGandGMRESSolver : public GMRES, public BCG
	{
		myvector::MyVector solve(myvector::MyVector U_begin, double &normL2u, double &normL2p, slae::SLAE& slae_in, logger::Logger& my_logger);
	};

}
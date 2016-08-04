#include "solver.h"

using namespace solverparameters;
using namespace std;
using namespace myvector;
using namespace slae;
using namespace densematrix;

namespace solver
{
	Solver::Solver(){}

	Solver::Solver(ifstream& s_parameters_f)
	{
		s_parameters_f >> s_parameters;
	}

	Solver::~Solver(){}

	void Solver::si_print(ofstream& log_f, 
						  int iteration_number,
						  double &normL2u,
						  double &normL2p,
						  SLAE& slae_in)
	{
		string f_name_s, f_name_i;

		log_f << "---" << iteration_number << "---" << endl;
		f_name_s = string("s_") + to_string(iteration_number) + ".txt";
		f_name_i = string("i_") + to_string(iteration_number) + ".txt";
		ofstream solution_f_out(f_name_s), info_f_out(f_name_i);

		output(solution_f_out, info_f_out, normL2u, normL2p);
		solution_f_out.close();
		info_f_out.close();

	}

	double Solver::find_relaxation_parameter(MyVector q_current, 
											 MyVector q_previous, 
											 double &residual_previous,
											 SLAE& slae_in)
	{
		MyVector qNonL(n);
		double w = 1, residual;

		do
		{
			qNonL = q_current * w + q_previous * (1 - w);
			slae_in.reinitialize();
			slae_in.calculate_global_matrix(qNonL);
			residual = (slae_in.A.b - slae_in.A * qNonL).norm() / slae_in.A.b.norm();
			printf("w = %.3lf -> residual = %.4le\n", w, residual);
			if(w > 0.1) w -= 0.1;
			else w -= 0.01;		
		}
		while(residual > residual_previous && w > 0.01);

		residual_previous = residual;
		return w;
	}

	void Solver::simple_iterations(SLAE& slae_in)
	{
		//FILE *in_f;
		//in_f = fopen("tmp.txt", "w");
		int k_it = 0;
		double w, residual_previous;
		MyVector q0(n), q1(n);
		double normL2u, normL2p;
		MyVector sol_0;
		sol_0.initialize(n);

		slae_in.A.create_portret();

		printf("\n%d iteration is in process\n", k_it);

		slae_in.A.calculate_global_matrix(sol_0);
		//for(int i = 0; i < A.ggl.size(); i++)
			//fprintf(in_f, "%.20lf\n", A.ggl[i]);
	
		Solve(sol_0, normL2u, normL2p);
		q0 = q_prev;

		get_vector_solution_in_nodes_ux(q0, Ux_numerical);
		get_vector_solution_in_nodes_uy(q0, Uy_numerical);
		get_vector_solution_in_nodes_p(q0, P_numerical);

		normL2u = diff_normL2_u(q0);
		normL2p = diff_normL2_p(q0);

		printf("%le\n", normL2u);
		printf("%le\n", normL2p);

		si_print(logger.log_f, k_it, normL2u, normL2p);

		do
		{
			k_it++;

			printf("\n%d iteration is in process\n", k_it);

			reinitialize();
			calculate_global_matrix(q0);
			Solve(q0, normL2u, normL2p);
			q1 = q_prev;

			if(k_it != 1 && k_it != 2)
			{
				w = find_relaxation_parameter(q1, q0, residual_previous);
				q1 = q1 * w + q0 * (1 - w);
			}
			else
			{
				reinitialize();
				calculate_global_matrix(q1);
				residual_previous = (A.b - A * q1).norm() / A.b.norm();
				///q1.output(in_f);
				printf("\n%.20lf residual\n", residual_previous);
			}

			q0 = q1;
			q_prev = q1;

			get_vector_solution_in_nodes_ux(q0, Ux_numerical);
			get_vector_solution_in_nodes_uy(q0, Uy_numerical);
			get_vector_solution_in_nodes_p(q0, P_numerical);

			normL2u = diff_normL2_u(q0);
			normL2p = diff_normL2_p(q0);

			printf("%le\n", normL2u);
			printf("%le\n", normL2p);

			si_print(logger.log_f, k_it, normL2u, normL2p);
		} while(k_it <= max_iter_nonlinear && residual_previous > 1e-6);
	}

	void Solver::run(ofstream& solution_f_out,
					 ofstream& info_f_out,
					 SLAE& slae_in)
	{
		bool one_research = true;
		double normL2u, normL2p;
		MyVector sol_0;
		sol_0.initialize(n);

		create_portret();

		calculate_global_matrix(sol_0);
		Solve(sol_0, normL2u, normL2p);
		output(solution_f_out, info_f_out, normL2u, normL2p);
	}
}
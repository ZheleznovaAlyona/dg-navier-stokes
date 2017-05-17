#include "main_solver.h"
#include <iostream>
#include <iomanip>
#include <math.h>
using namespace std;

#include "element.h"
#include "partition.h"
#include "testing_parameters.h"
#include "parameters.h"
#include "myfunctions.h"

using namespace boundary_conditions;
using namespace element;
using namespace partition;
using namespace myvector;
using namespace matrix;
using namespace boundaries;
using namespace logger;
using namespace parameters;
using namespace solver;
using namespace slae;
using namespace testingparameters;
using namespace integration;


namespace mainsolver
{
	MainSolver::MainSolver() {}

	MainSolver::MainSolver(std::ifstream& grid_f_in,
						   std::ifstream& elements_f_in,
						   std::string log_f,
						   std::ifstream& boundary1,
						   std::ifstream& boundary2,
						   std::ifstream& boundary3,
						   std::string fileInPenaltyParameters)
	{
		initialize(grid_f_in,
			elements_f_in,
			log_f,
			boundary1,
			boundary2,
			boundary3,
			fileInPenaltyParameters);
	}

	MainSolver::~MainSolver() {}

	void MainSolver::initialize(ifstream& grid_f_in,
								ifstream& elements_f_in,
								string log_f,
								ifstream& boundary1,
								ifstream& boundary2,
								ifstream& boundary3,
								std::string fileInPenaltyParameters)
	{
		boundary1 >> boundaries1;
		boundary2 >> boundaries2;
		boundary3 >> boundaries3;

		Partition::input(grid_f_in, elements_f_in);
		Basis::initialize();
		Gauss_integration::initialize();
		int slae_size = (basis::n_func_p + basis::n_func_u) * elements.size();

		InternalBoundaries::initialize_penalty_parameters(fileInPenaltyParameters);
		OuterBoundaries::initialize_penalty_parameters(fileInPenaltyParameters);
		BoundaryConditionsSupport::initialize_penalty_parameters(fileInPenaltyParameters);

		logger.open(log_f);

		int tmp = Partition::count_unzero_matrix_elements();
		my_slae.initialize(slae_size, tmp);

		Ux_numerical.initialize(nodes.size());
		Uy_numerical.initialize(nodes.size());
		P_numerical.initialize(nodes.size());
		q_prev.initialize(slae_size);
	}

	void MainSolver::reinitialize()
	{
		Ux_numerical.make_zero();
		Uy_numerical.make_zero();
		P_numerical.make_zero();
		my_slae.reinitialize();
	}

	void MainSolver::build_slae(MyVector q_calc)
	{
		logger.send_message_build_slae();
		int size = Partition::elements.size();

		//локаьные матрицы и вектор правой части
		for(int el_i = 0; el_i < size; el_i++)
			calculate_locals(el_i, q_calc);

		//матрицы межэлементных границ
		for(int el_i = 0; el_i < size; el_i++)
			InternalBoundaries::calculate_internal_boundaries(el_i, my_slae.A);

		//расчёт матриц границ области
			for(int el_i = 0; el_i < size; el_i++)
				OuterBoundaries::calculate_outer_boundaries(el_i, my_slae.A);

		//учёт первых краевых условий
		BoundaryConditionsSupport::calculate_all_boundaries1(my_slae.b);
	}

	double MainSolver::find_relaxation_parameter(MyVector q_current, 
												MyVector q_previous, 
												double &residual_previous)
	{
		MyVector qNonL(my_slae.n);
		double w = 1, residual;

		do
		{
			qNonL = q_current * w + q_previous * (1 - w);
			my_slae.reinitialize();
			build_slae(qNonL);
			residual = (my_slae.b - my_slae.A * qNonL).norm() / my_slae.b.norm();
			cout << "w = " << setprecision(3) << w << " -> residual = " << 
							  scientific << setprecision(4) << residual << endl;
			if(w > 0.1) w -= 0.1;
			else w -= 0.01;		
		}
		while(residual > residual_previous && w > 0.01);

		residual_previous = residual;
		return w;
	}


	void MainSolver::simple_iterations(Solver& s)
	{
		int k_it = 0;
		double w, residual_previous;
		MyVector q0(my_slae.n), q1(my_slae.n);
		double normL2u, normL2p;
		MyVector sol_0;
		sol_0.initialize(my_slae.n);

		static auto& part = static_cast<Partition>(*this);

		my_slae.A.create_portret(part, logger);

		logger.send_current_information_to_screen_si(k_it);

		build_slae(sol_0);

		logger.send_message_solution();
		q_prev = s.solve(sol_0, normL2u, normL2p, my_slae, logger);
		q0 = q_prev;

		get_vector_solution_in_nodes_ux(q0, Ux_numerical);
		get_vector_solution_in_nodes_uy(q0, Uy_numerical);
		get_vector_solution_in_nodes_p(q0, P_numerical);

		logger.send_message_norms();
		normL2u = diff_normL2_u(q0);
		logger.send_inf_UL2_norm(normL2u);
		normL2p = diff_normL2_p(q0);	
		logger.send_inf_PL2_norm(normL2p);

		logger.si_print(k_it, normL2u, normL2p, Ux_numerical, Uy_numerical, P_numerical, part);

		do
		{
			k_it++;

			logger.send_current_information_to_screen_si(k_it);

			reinitialize();

			build_slae(q0);

			logger.send_message_solution();
			q_prev = s.solve(q0, normL2u, normL2p, my_slae, logger);
			q1 = q_prev;

			if(k_it != 1 && k_it != 2)
			{
				w = find_relaxation_parameter(q1, q0, residual_previous);
				q1 = q1 * w + q0 * (1 - w);
			}
			else
			{
				reinitialize();

				build_slae(q1);

				residual_previous = (my_slae.b - my_slae.A * q1).norm() / my_slae.b.norm();
				cout << endl << "residual = " << setprecision(20) << residual_previous  << endl;
			}

			q0 = q1;
			q_prev = q1;

			get_vector_solution_in_nodes_ux(q0, Ux_numerical);
			get_vector_solution_in_nodes_uy(q0, Uy_numerical);
			get_vector_solution_in_nodes_p(q0, P_numerical);

			logger.send_message_norms();
			normL2u = diff_normL2_u(q0);
			logger.send_inf_UL2_norm(normL2u);
			normL2p = diff_normL2_p(q0);
			logger.send_inf_PL2_norm(normL2p);

			logger.si_print(k_it, normL2u, normL2p, Ux_numerical, Uy_numerical, P_numerical, part);
		} while(k_it <= s.s_parameters.max_number_of_iterations_non_lin && residual_previous > 1e-6);
	}

	void MainSolver::linear(ofstream& solution_f_out,
							ofstream& info_f_out,
							Solver& s)
	{
		double normL2u, normL2p;
		MyVector sol_0;
		sol_0.initialize(my_slae.n);

		my_slae.A.create_portret(static_cast<Partition>(*this), logger);

		build_slae(sol_0);

		logger.send_message_solution();
		q_prev = s.solve(sol_0, normL2u, normL2p, my_slae, logger);

		get_vector_solution_in_nodes_ux(q_prev, Ux_numerical);
		get_vector_solution_in_nodes_uy(q_prev, Uy_numerical);
		get_vector_solution_in_nodes_p(q_prev, P_numerical);

		logger.send_message_norms();
		normL2u = diff_normL2_u(q_prev); 
		logger.send_inf_UL2_norm(normL2u);
		normL2p = diff_normL2_p(q_prev);		
		logger.send_inf_PL2_norm(normL2p);

		logger.output(solution_f_out, info_f_out, normL2u, normL2p, Ux_numerical, Uy_numerical, P_numerical, static_cast<Partition>(*this));
		//TO DO: output in tecplot
	}


void MainSolver::solve()
{
	MyVector q(my_slae.n);
	Solver *s;

	switch(Testing_parameters::solver)
	{
	case 1:
		{
			s = new BiCGStab();
		}
		break;
	case 2:
		{
			s = new GMRES();
		}
		break;
	case 3:
		{
			s = new BCGandGMRESSolver();
		}
		break;
	case 4:
		{
			s = new BCG();
		}
		break;
	default:
		{
			s = new GMRES();
		}
	};

	s->s_parameters.initialize("solver.json");

	simple_iterations(*s);

	q_prev = q;
}

void MainSolver::solve(std::ofstream & solution_f_out, std::ofstream & info_f_out)
{
	MyVector q(my_slae.n);
	Solver *s;

	switch (Testing_parameters::solver)
	{
	case 1:
		{
			s = new BiCGStab();
		}
		break;
	case 2:
		{
			s = new GMRES();
		}
		break;
	case 3:
		{
			s = new BCGandGMRESSolver();
		}
	break;
	case 4:
		{
			s = new BCG();
		}
	break;
	default:
		{
			s = new GMRES();
		}
	};

	if(Testing_parameters::use_LU) my_slae.A.LU();

	s->s_parameters.initialize("solver.json");

	linear(solution_f_out, info_f_out, *s);
}

double MainSolver::diff_normL2_p(MyVector q_solution)
{
	double diff_local, diff;
	int size = elements.size();
	double p, function;
	double x0, y0, hx, hy;
	double jacobian;

	diff = 0;
	for(int i = 0; i < size; i++)
	{
		x0 = nodes[elements[i].nodes[0]].x;
		y0 = nodes[elements[i].nodes[0]].y;
		hx = get_hx(i);
		hy = get_hy(i);
		jacobian = hx * hy / 4.0;

		diff_local = 0;
		for(int k = 0; k < n_ip; k++)
		{
			double p_x = hx * gauss_points[0][k] + x0;
			double p_y = hy * gauss_points[1][k] + y0;
			p = get_solution_in_point_p(p_x, p_y, i, q_solution);
			function = calculate_p_analytic(elements[i].number_of_area, p_x, p_y);
			function -= p;
			diff_local += gauss_weights[k] * function * function;
		}
		diff_local *= jacobian;

		diff += diff_local;
	}

	return sqrt(diff / size);
}

double MainSolver::diff_normL2_u(MyVector q_solution)
{
	double diff_local, diff;
	int size = elements.size();
	double function;
	double x0, y0, hx, hy;
	double jacobian;
	double ux, uy, uxdx, uydy;
	double ux_an, uy_an, uxdx_an, uydy_an;

	diff = 0;
	for(int i = 0; i < size; i++)
	{
		x0 = nodes[elements[i].nodes[0]].x;
		y0 = nodes[elements[i].nodes[0]].y;
		hx = get_hx(i);
		hy = get_hy(i);
		jacobian = hx * hy / 4.0;

		diff_local = 0;
		for(int k = 0; k < n_ip; k++)
		{
			double p_x = hx * gauss_points[0][k] + x0;
			double p_y = hy * gauss_points[1][k] + y0;

			ux = get_solution_in_point_ux(p_x, p_y, i, q_solution);
			uy = get_solution_in_point_uy(p_x, p_y, i, q_solution);
			uxdx = get_solution_in_point_uxdx(p_x, p_y, i, q_solution);
			uydy = get_solution_in_point_uydy(p_x, p_y, i, q_solution);
			ux_an = calculate_ux_analytic(elements[i].number_of_area, p_x, p_y);
			uy_an = calculate_uy_analytic(elements[i].number_of_area, p_x, p_y);
			uxdx_an = calculate_uxdx_analytic(elements[i].number_of_area, p_x, p_y);
			uydy_an = calculate_uydy_analytic(elements[i].number_of_area, p_x, p_y);

			function = (ux - ux_an) * (ux - ux_an) + (uy - uy_an) * (uy - uy_an);
			//function += (uxdx - uxdx_an + uydy - uydy_an) * 
						//(uxdx - uxdx_an + uydy - uydy_an);

			diff_local += gauss_weights[k] * function;
		}
		diff_local *= jacobian;

		diff += diff_local;
	}

	return sqrt(diff / size);
}

}

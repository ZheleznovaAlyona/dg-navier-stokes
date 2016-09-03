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


namespace mainsolver
{
	MainSolver::MainSolver() {}

	MainSolver::MainSolver(std::ifstream& grid_f_in,
		std::ifstream& elements_f_in,
		std::string log_f,
		std::ifstream& boundary1,
		std::ifstream& boundary2,
		std::ifstream& boundary3)
	{
		initialize(grid_f_in,
			elements_f_in,
			log_f,
			boundary1,
			boundary2,
			boundary3);
	}

	MainSolver::~MainSolver() {}

	void MainSolver::initialize(ifstream& grid_f_in,
								ifstream& elements_f_in,
								string log_f,
								ifstream& boundary1,
								ifstream& boundary2,
								ifstream& boundary3)
	{
		boundary1 >> boundaries1;
		boundary2 >> boundaries2;
		boundary3 >> boundaries3;

		Partition::input(grid_f_in, elements_f_in);

		int slae_size = nodes.size() + elements.size() * 4; //узлы и рёбра

		InternalBoundaries::initialize_penalty_parameters();
		OuterBoundaries::initialize_penalty_parameters();
		BoundaryConditionsSupport::initialize_penalty_parameters();

		logger.open(log_f);

		int tmp = Partition::count_unzero_matrix_elements();
		my_slae.initialize(slae_size, tmp);

		Ux_numerical.initialize(nodes.size());
		Uy_numerical.initialize(nodes.size());
		P_numerical.initialize(nodes.size());
		q_prev.initialize(slae_size);

		Basis::initialize();
		Gauss_integration::initialize();
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

		my_slae.A.create_portret(part, logger, n_func_u, n_func_p);

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
		bool one_research = true;
		double normL2u, normL2p;
		MyVector sol_0;
		sol_0.initialize(my_slae.n);

		my_slae.A.create_portret(static_cast<Partition>(*this), logger, n_func_u, n_func_p);

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
	}


#pragma region локальные матрицы и векторы

	void MainSolver::calculate_locals(int element_number, MyVector q_calc)
	{
		Element element = Partition::elements[element_number];
	
		calculate_G(element_number);
		calculate_C(element_number, q_calc);
		calculate_P1(element_number);
		calculate_P2(element_number);
		calculate_F(element_number);
	}

	void MainSolver::calculate_G(int element_number)
	{
		double hx, hy, hx2, hy2, g1, g2;
		vector <vector<double>> G;
		G.resize(n_func_u);
		for (int i = 0; i < n_func_u; i++)
			initialize_vector(G[i], n_func_u);
		Element element = Partition::elements[element_number];
		double lambda = calculate_lambda(element.number_of_area);

		hx = get_hx(element_number);
		hy = get_hy(element_number);
		hx2 = hx * hx;
		hy2 = hy * hy;

		double jacobian = hx * hy / 4.0;

		for(int i = 0; i < n_func_u; i++)
		{
			for(int j = i; j < n_func_u; j++)
			{
				g1 = 0;
				g2 = 0;
				for(int k = 0; k < 9; k++)
				{
					double p_ksi = gauss_points[0][k], 
						   p_etta = gauss_points[1][k];
					g1 += gauss_weights[k] * 
						  (dphixksi[i](p_ksi, p_etta) * dphixksi[j](p_ksi, p_etta) +
						   dphiyksi[i](p_ksi, p_etta) * dphiyksi[j](p_ksi, p_etta));
					g2 += gauss_weights[k] * 
						  (dphixetta[i](p_ksi, p_etta) * dphixetta[j](p_ksi, p_etta) +
						  dphiyetta[i](p_ksi, p_etta) * dphiyetta[j](p_ksi, p_etta));
				}
				G[i][j] = (g1 * jacobian / hx2 + g2 * jacobian / hy2) * lambda;
			}
		}

		for(int i = 1; i < 4; i++)
			for(int j = 0; j < i; j++)
				G[i][j] = G[j][i];

		for (int i = 0; i < n_func_u; i++)
		{
			int id_i = element.edges[i];
			for (int j = 0; j < n_func_u; j++)
			{
				int id_j = element.edges[j];
				my_slae.A.add_element(id_i, id_j, G[i][j]);
			}
		}
	}
	void MainSolver::calculate_C(int element_number, MyVector q_calc)
	{
		double hx, hy, c1, c2;
		Element element = elements[element_number];

		vector <vector<double>> C;
		C.resize(n_func_u);
		for (int i = 0; i < n_func_u; i++)
			initialize_vector(C[i], n_func_u);

		hx = get_hx(element_number);
		hy = get_hy(element_number);

		double jacobian = hx * hy / 4.0;
		double u_x, u_y;
		double x0 = nodes[element.nodes[0]].x;
		double y0 = nodes[element.nodes[0]].y;

		for(int i = 0; i < n_func_u; i++)
		{
			for(int j = 0; j < n_func_u; j++)
			{
				C[i][j] = 0;
				for(int k = 0; k < 9; k++)
				{
					double p_ksi = gauss_points[0][k], 
						   p_etta = gauss_points[1][k];
					double p_x = p_ksi * hx + x0, p_y = p_etta * hy + y0;
					u_x = get_solution_in_point_ux(p_x, p_y, element_number, q_calc);
					u_y = get_solution_in_point_uy(p_x, p_y, element_number, q_calc);
					//c1 = u_x * (dphixksi[j](p_ksi, p_etta) / hx * phix[i](p_ksi, p_etta) +
					//		    dphixetta[j](p_ksi, p_etta) / hy * phiy[i](p_ksi, p_etta));
					//c2 = u_y * (dphiyksi[j](p_ksi, p_etta) / hx * phix[i](p_ksi, p_etta) +
					//		    dphiyetta[j](p_ksi, p_etta) / hy * phiy[i](p_ksi, p_etta));
					c1 = u_x * (dphixksi[j](p_ksi, p_etta) / hx * phix[i](p_ksi, p_etta) +
								dphiyksi[j](p_ksi, p_etta) / hx * phiy[i](p_ksi, p_etta)-
								dphixksi[i](p_ksi, p_etta) / hx * phix[j](p_ksi, p_etta) -
								dphiyksi[i](p_ksi, p_etta) / hx * phiy[j](p_ksi, p_etta));
					c2 = u_y * (dphixetta[j](p_ksi, p_etta) / hy * phix[i](p_ksi, p_etta) +
								dphiyetta[j](p_ksi, p_etta) / hy * phiy[i](p_ksi, p_etta)-
								dphixetta[i](p_ksi, p_etta) / hy * phix[j](p_ksi, p_etta) -
								dphiyetta[i](p_ksi, p_etta) / hy * phiy[j](p_ksi, p_etta));			
					C[i][j] += gauss_weights[k] * 0.5 * (c1 + c2);
				}
				C[i][j] *= jacobian;
			}
		}

		for (int i = 0; i < n_func_u; i++)
		{
			int id_i = element.edges[i];
			for (int j = 0; j < n_func_u; j++)
			{
				int id_j = element.edges[j];
				my_slae.A.add_element(id_i, id_j, C[i][j]);
			}
		}
	}
	void MainSolver::calculate_P1(int element_number)
	{
		int n_edges = Partition::elements.size() * 4;
		Element element = elements[element_number];

		vector <vector<double>> P1;
		P1.resize(n_func_u);
		for (int i = 0; i < n_func_u; i++)
			initialize_vector(P1[i], n_func_p);

		double hx = get_hx(element_number);
		double hy = get_hy(element_number);

		double jacobian = hx * hy / 4.0;

		double rho = calculate_rho(element.number_of_area);
		rho = - 1 / rho;

		for(int i = 0; i < n_func_u; i++)
		{
			for(int j = 0; j < n_func_p; j++)
			{
				P1[i][j] = 0;
				for(int k = 0; k < 9; k++)
				{
					double p_ksi = gauss_points[0][k], p_etta = gauss_points[1][k];
					P1[i][j] += gauss_weights[k] * psi[j](p_ksi, p_etta) *
						   		(dphixksi[i](p_ksi, p_etta) / hx + 
								dphiyetta[i](p_ksi, p_etta) / hy);
				}
				P1[i][j] *= jacobian * rho;
			}
		}

		for (int i = 0; i < n_func_u; i++)
		{
			int id_i = element.edges[i];
			for (int j = 0; j < n_func_p; j++)
			{
				int id_j = element.nodes[j] + n_edges;
				my_slae.A.add_element(id_i, id_j, P1[i][j]);
			}
		}
	}
	void MainSolver::calculate_P2(int element_number)
	{
		int n_edges = Partition::elements.size() * 4;
		Element element = elements[element_number];
		vector <vector<double>> P2;
		P2.resize(n_func_p);
		for (int i = 0; i < n_func_p; i++)
			initialize_vector(P2[i], n_func_u);

		double hx = get_hx(element_number);
		double hy = get_hy(element_number);

		double jacobian = hx * hy / 4.0;

		for(int i = 0; i < n_func_p; i++)
		{
			for(int j = 0; j < n_func_u; j++)
			{
				P2[i][j] = 0;
				for(int k = 0; k < 9; k++)
				{
					double p_ksi = gauss_points[0][k], p_etta = gauss_points[1][k];
					P2[i][j] += gauss_weights[k] * psi[i](p_ksi, p_etta) *
						   		(dphixksi[j](p_ksi, p_etta) / hx + 
								dphiyetta[j](p_ksi, p_etta) / hy);
				}
				P2[i][j] *= jacobian;
			}
		}

		for (int i = 0; i < n_func_p; i++)
		{
			int id_i = element.nodes[i] + n_edges;
			for (int j = 0; j < n_func_u; j++)
			{
				int id_j = element.edges[j];
				my_slae.A.add_element(id_i, id_j, P2[i][j]);
			}
		}
	}
	void MainSolver::calculate_F(int element_number)
	{
		double f_x, f_y;
		Element element = elements[element_number];
		int number_of_area = element.number_of_area;

		vector<double> F;
		initialize_vector(F, n_func_u);

		double x0 = nodes[element.nodes[0]].x;
		double y0 = nodes[element.nodes[0]].y;

		double hx = get_hx(element_number);
		double hy = get_hy(element_number);

		double jacobian = hx * hy / 4.0;

		for(int i = 0;  i < n_func_u; i++)
		{
			F[i] = 0;
			for(int j = 0; j < 9; j++)
			{
				double p_ksi = gauss_points[0][j], p_etta = gauss_points[1][j];
				double p_x = p_ksi * hx + x0, p_y = p_etta * hy + y0;
				f_x = calculate_fx(number_of_area, p_x, p_y);	
				f_y = calculate_fy(number_of_area, p_x, p_y);
				F[i] += gauss_weights[j] * (f_x * phix[i](p_ksi, p_etta) 
						+ f_y * phiy[i](p_ksi, p_etta));
			}
			 F[i] *= jacobian;
		}

		for (int i = 0; i < n_func_u; i++)
		{
			int id_i = element.edges[i];
			my_slae.b[id_i] += F[i];
		}
	}

#pragma endregion

#pragma region расчёт решений и производных

double MainSolver::get_solution_in_point_ux(double x, double y, int element_number, MyVector qi)
{
	vector <int> indexes;
	indexes.resize(n_func_u);

	vector <double> qi_local;
	qi_local.resize(n_func_u);

	double u_in_point;

	//собираем глобальные номера с элемента
	for(int j = 0; j < n_func_u; j++)
		indexes[j] = elements[element_number].edges[j];

	//собираем локальный набор весов
	for(int j = 0; j < n_func_u; j++)
		qi_local[j] = qi[indexes[j]];

	//вычисляем в решение в точке
	static auto& part = static_cast<Partition>(*this);
	u_in_point = 0;
	for(int j = 0; j < n_func_u; j++)
		u_in_point += qi_local[j] * phix_i(j, x, y, element_number, part);

	return u_in_point;
}

double MainSolver::get_solution_in_point_uy(double x, double y, int element_number, MyVector qi)
{
	vector <int> indexes;
	indexes.resize(n_func_u);

	vector <double> qi_local;
	qi_local.resize(n_func_u);

	double u_in_point;

	//собираем глобальные номера с элемента
	for(int j = 0; j < n_func_u; j++)
		indexes[j] = elements[element_number].edges[j];

	//собираем локальный набор весов
	for(int j = 0; j < n_func_u; j++)
		qi_local[j] = qi[indexes[j]];

	//вычисляем в решение в точке
	static auto& part = static_cast<Partition>(*this);
	u_in_point = 0;
	for(int j = 0; j < n_func_u; j++)
		u_in_point += qi_local[j] * phiy_i(j, x, y, element_number, part);

	return u_in_point;
}

double MainSolver::get_solution_in_point_uxdx(double x, double y, int element_number, MyVector qi)
{
	vector <int> indexes;
	indexes.resize(n_func_u);

	vector <double> qi_local;
	qi_local.resize(n_func_u);

	double du_in_point;

	//собираем глобальные номера с элемента
	for(int j = 0; j < n_func_u; j++)
		indexes[j] = elements[element_number].edges[j];

	//собираем локальный набор весов
	for(int j = 0; j < n_func_u; j++)
		qi_local[j] = qi[indexes[j]];

	//вычисляем в решение в точке
	static auto& part = static_cast<Partition>(*this);
	du_in_point = 0;
	for(int j = 0; j < n_func_u; j++)
		du_in_point += qi_local[j] * phixdx_i(j, x, y, element_number, part);

	return du_in_point;
}

double MainSolver::get_solution_in_point_uydy(double x, double y, int element_number, MyVector qi)
{
	vector <int> indexes;
	indexes.resize(n_func_u);

	vector <double> qi_local;
	qi_local.resize(n_func_u);

	double du_in_point;

	//собираем глобальные номера с элемента
	for(int j = 0; j < n_func_u; j++)
		indexes[j] = elements[element_number].edges[j];

	//собираем локальный набор весов
	for(int j = 0; j < n_func_u; j++)
		qi_local[j] = qi[indexes[j]];

	//вычисляем в решение в точке
	static auto& part = static_cast<Partition>(*this);
	du_in_point = 0;
	for(int j = 0; j < n_func_u; j++)
		du_in_point += qi_local[j] * phiydy_i(j, x, y, element_number, part);
	
	return du_in_point;
}

double MainSolver::get_solution_in_point_p(double x, double y, int element_number, MyVector qi)
{
	int n_edges = elements.size() * 4;

	vector <int> indexes;
	indexes.resize(n_func_p);

	vector <double> qi_local;
	qi_local.resize(n_func_p);

	double p_in_point;

	//собираем глобальные номера с элемента
	for(int j = 0; j < n_func_p; j++)
		indexes[j] = elements[element_number].nodes[j] + n_edges;

	//собираем локальный набор весов
	for(int j = 0; j < n_func_p; j++)
		qi_local[j] = qi[indexes[j]];

	//вычисляем в решение в точке
	static auto& part = static_cast<Partition>(*this);
	p_in_point = 0;
	for(int j = 0; j < n_func_p; j++)
		p_in_point += qi_local[j] * psi_i(j, x, y, element_number, part);

	return p_in_point;
}

double MainSolver::get_solution_in_point2_ux(double x, double y, MyVector qi)
{
	int element_number = search_element(x, y);
	return get_solution_in_point_ux(x, y, element_number, qi);
}

double MainSolver::get_solution_in_point2_uy(double x, double y, MyVector qi)
{
	int element_number = search_element(x, y);
	return get_solution_in_point_uy(x, y, element_number, qi);
}

double MainSolver::get_solution_in_point2_p(double x, double y, MyVector qi)
{
	int element_number = search_element(x, y);
	return get_solution_in_point_p(x, y, element_number, qi);
}

void MainSolver::get_vector_solution_in_nodes_ux(MyVector qi, MyVector &solution)
{
	logger.send_message_Ux();
	int indexes_nodes[4];
	int size = elements.size();
	double u_local[4];
	double x, y;

	vector <int> indexes;
	indexes.resize(n_func_u);

	vector <double> qi_local;
	qi_local.resize(n_func_u);

	static auto& part = static_cast<Partition>(*this);

	for(int i = 0; i < size; i++)
	{
		//собираем глобальные номера с элемента
		for(int j = 0; j < n_func_u; j++)
			indexes[j] = elements[i].edges[j];

		for(int j = 0; j < 4; j++)
			indexes_nodes[j] = elements[i].nodes[j];

		//собираем локальный набор весов
		for(int j = 0; j < n_func_u; j++)
			qi_local[j] = qi[indexes[j]];

		//вычисляем в узлах элемента решение
		for(int j = 0; j < 4; j++)
		{
			u_local[j] = 0;
			x = nodes[indexes_nodes[j]].x;
			y = nodes[indexes_nodes[j]].y;
			for(int k = 0; k < n_func_u; k++)
				u_local[j] += qi_local[k] * phix_i(k, x, y, i, part);
		}

		//кладём в результирующий вектор
		for(int j = 0; j < 4; j++)
			solution[indexes_nodes[j]] = u_local[j];
	}
}

void MainSolver::get_vector_solution_in_nodes_uy(MyVector qi, MyVector &solution)
{
	logger.send_message_Uy();
	int indexes_nodes[4];
	int size = elements.size();
	double u_local[4];
	double x, y;

	vector <int> indexes;
	indexes.resize(n_func_u);

	vector <double> qi_local;
	qi_local.resize(n_func_u);

	static auto& part = static_cast<Partition>(*this);

	for(int i = 0; i < size; i++)
	{
		//собираем глобальные номера с элемента
		for(int j = 0; j < n_func_u; j++)
			indexes[j] = elements[i].edges[j];

		//собираем локальный набор весов
		for(int j = 0; j < n_func_u; j++)
			qi_local[j] = qi[indexes[j]];

		for(int j = 0; j < 4; j++)
			indexes_nodes[j] = elements[i].nodes[j];

		//вычисляем в узлах элемента решение
		for(int j = 0; j < 4; j++)
		{
			u_local[j] = 0;
			x = nodes[indexes_nodes[j]].x;
			y = nodes[indexes_nodes[j]].y;
			for(int k = 0; k < n_func_u; k++)
				u_local[j] += qi_local[k] * phiy_i(k, x, y, i, part);
		}

		//кладём в результирующий вектор
		for(int j = 0; j < 4; j++)
			solution[indexes_nodes[j]] = u_local[j];
	}
}

void MainSolver::get_vector_solution_in_nodes_p(MyVector qi, MyVector &solution)
{
	logger.send_message_P();
	int indexes_nodes[4], n_edges = elements.size() * 4;
	int size = elements.size();
	double p_local[4];
	double x, y;

	vector <int> indexes;
	indexes.resize(n_func_p);

	vector <double> qi_local;
	qi_local.resize(n_func_p);

	static auto& part = static_cast<Partition>(*this);

	for(int i = 0; i < size; i++)
	{
		//собираем глобальные номера с элемента
		for(int j = 0; j < n_func_p; j++)
			indexes[j] = elements[i].nodes[j];

		//собираем локальный набор весов
		for(int j = 0; j < n_func_p; j++)
			qi_local[j] = qi[indexes[j]+ n_edges];

		for (int j = 0; j < 4; j++)
			indexes_nodes[j] = elements[i].nodes[j];

		//вычисляем в узлах элемента решение
		for(int j = 0; j < 4; j++)
		{
			p_local[j] = 0;
			x = nodes[indexes_nodes[j]].x;
			y = nodes[indexes_nodes[j]].y;
			for(int k = 0; k < n_func_p; k++)
				p_local[j] += qi_local[k] * psi_i(k, x, y, i, part);
		}

		//кладём в результирующий вектор
		for(int j = 0; j < 4; j++)
			solution[indexes_nodes[j]] = p_local[j];
	}
}

#pragma endregion

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
		for(int k = 0; k < 9; k++)
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

	return sqrt(diff / nodes.size());
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
		for(int k = 0; k < 9; k++)
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
			function += (uxdx - uxdx_an + uydy - uydy_an) * 
						(uxdx - uxdx_an + uydy - uydy_an);

			diff_local += gauss_weights[k] * function;
		}
		diff_local *= jacobian;

		diff += diff_local;
	}

	return sqrt(diff / nodes.size());
}

}

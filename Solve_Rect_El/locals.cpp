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
#include "point.h"

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
using namespace point;
using namespace integration;

namespace mainsolver
{

	void MainSolver::calculate_locals(int element_number, MyVector& q_calc)
	{
		Element element = Partition::elements[element_number];

		calculate_G(element_number);
		//calculate_C(element_number, q_calc);
		calculate_P1(element_number);
		calculate_P2(element_number);
		calculate_F(element_number);
	}

	void MainSolver::calculate_G(int element_number)
	{
		Element element = Partition::elements[element_number];
		double hx, hy, hx2, hy2, g1, g2;
		vector <vector<double>> G;
		G.resize(element.ndof_u);
		for (int i = 0; i < element.ndof_u; i++)
			initialize_vector(G[i], element.ndof_u);
		
		double lambda = calculate_lambda(element.number_of_area);

		hx = get_hx(element_number);
		hy = get_hy(element_number);
		hx2 = hx * hx;
		hy2 = hy * hy;

		double jacobian = hx * hy / 4.0;

		for (int i = 0; i < element.ndof_u; i++)
		{
			for (int j = i; j < element.ndof_u; j++)
			{
				g1 = 0;
				g2 = 0;
				for (int k = 0; k < n_ip; k++)
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

		for (int i = 1; i < element.ndof_u; i++)
			for (int j = 0; j < i; j++)
				G[i][j] = G[j][i];

		for (int i = 0; i < element.ndof_u; i++)
		{
			int id_i = element.dof_u[i];
			for (int j = 0; j < element.ndof_u; j++)
			{
				int id_j = element.dof_u[j];
				my_slae.A.add_element(id_i, id_j, G[i][j]);
			}
		}
	}
	void MainSolver::calculate_C(int element_number, MyVector& q_calc)
	{
		double hx, hy, c1, c2;
		Element element = elements[element_number];

		hx = get_hx(element_number);
		hy = get_hy(element_number);

		double jacobian = hx * hy / 4.0;
		double u_x, u_y;
		double x0 = nodes[element.nodes[0]].x;
		double y0 = nodes[element.nodes[0]].y;

		for (int i = 0; i < element.ndof_u; i++)
		{
			int id_i = element.dof_u[i];
			for (int j = 0; j < element.ndof_u; j++)
			{
				int id_j = element.dof_u[j];
				double cij = 0;
				for (int k = 0; k < n_ip; k++)
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
						        dphiyksi[j](p_ksi, p_etta) / hx * phiy[i](p_ksi, p_etta) -
						        dphixksi[i](p_ksi, p_etta) / hx * phix[j](p_ksi, p_etta) -
						        dphiyksi[i](p_ksi, p_etta) / hx * phiy[j](p_ksi, p_etta));
					c2 = u_y * (dphixetta[j](p_ksi, p_etta) / hy * phix[i](p_ksi, p_etta) +
								dphiyetta[j](p_ksi, p_etta) / hy * phiy[i](p_ksi, p_etta) -
								dphixetta[i](p_ksi, p_etta) / hy * phix[j](p_ksi, p_etta) -
								dphiyetta[i](p_ksi, p_etta) / hy * phiy[j](p_ksi, p_etta));
					cij += gauss_weights[k] * 0.5 * (c1 + c2);
				}
				cij *= jacobian;
				my_slae.A.add_element(id_i, id_j, cij);
			}
		}
	}
	void MainSolver::calculate_P1(int element_number)
	{
		Element element = elements[element_number];

		double hx = get_hx(element_number);
		double hy = get_hy(element_number);

		double jacobian = hx * hy / 4.0;

		double rho = calculate_rho(element.number_of_area);
		rho = -1 / rho;

		for (int i = 0; i < element.ndof_u; i++)
		{
			int id_i = element.dof_u[i];
			for (int j = 0; j < element.ndof_p; j++)
			{
				int id_j = element.dof_p[j];
				double p1ij = 0;
				for (int k = 0; k < n_ip; k++)
				{
					double p_ksi = gauss_points[0][k], p_etta = gauss_points[1][k];
					p1ij += gauss_weights[k] * psi[j](p_ksi, p_etta) *
								(dphixksi[i](p_ksi, p_etta) / hx +
								 dphiyetta[i](p_ksi, p_etta) / hy);
				}
				p1ij *= jacobian * rho;
				my_slae.A.add_element(id_i, id_j, p1ij);
			}
		}
	}
	void MainSolver::calculate_P2(int element_number)
	{
		Element element = elements[element_number];

		double hx = get_hx(element_number);
		double hy = get_hy(element_number);

		double jacobian = hx * hy / 4.0;

		for (int i = 0; i < element.ndof_p; i++)
		{
			int id_i = element.dof_p[i];
			for (int j = 0; j < element.ndof_u; j++)
			{
				int id_j = element.dof_u[j];
				double p2ij = 0;
				for (int k = 0; k < n_ip; k++)
				{
					double p_ksi = gauss_points[0][k], p_etta = gauss_points[1][k];
					p2ij += gauss_weights[k] * psi[i](p_ksi, p_etta) *
								(dphixksi[j](p_ksi, p_etta) / hx +
								 dphiyetta[j](p_ksi, p_etta) / hy);
				}
				p2ij *= jacobian;
				my_slae.A.add_element(id_i, id_j, p2ij);
			}
		}
	}
	void MainSolver::calculate_F(int element_number)
	{
		double f_x, f_y;
		Element element = elements[element_number];
		int number_of_area = element.number_of_area;

		double x0 = nodes[element.nodes[0]].x;
		double y0 = nodes[element.nodes[0]].y;

		double hx = get_hx(element_number);
		double hy = get_hy(element_number);

		double jacobian = hx * hy / 4.0;

		for (int i = 0; i < element.ndof_u; i++)
		{
			int id_i = element.dof_u[i];
			double fi = 0;
			for (int k = 0; k < n_ip; k++)
			{
				double p_ksi = gauss_points[0][k], p_etta = gauss_points[1][k];
				double p_x = p_ksi * hx + x0, p_y = p_etta * hy + y0;
				f_x = calculate_fx(number_of_area, p_x, p_y);
				f_y = calculate_fy(number_of_area, p_x, p_y);
				fi += gauss_weights[k] * (f_x * phix[i](p_ksi, p_etta)
										  + f_y * phiy[i](p_ksi, p_etta));
			}
			fi *= jacobian;
			my_slae.b[id_i] += fi;
		}
	}
}
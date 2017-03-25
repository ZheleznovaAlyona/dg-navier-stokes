#include "boundary_conditions.h"
#include "element.h"
#include "partition.h"
#include "myfunctions.h"

using namespace std;
using namespace element;
using namespace partition;
using namespace parameters;
using namespace myvector;
using namespace basis;

namespace boundary_conditions
{
	ifstream& operator>>(ifstream& is, vector <BoundaryCondition>& boundaries)
	{
		int count;
		BoundaryCondition tmp;

		is >> count;
		boundaries.reserve(count);

		for(int i = 1; i <= count; i++)
		{
			is >> tmp.elem;
			is >> tmp.edges[0];
			is >> tmp.edges[1];
			is >> tmp.edges[2];
			is >> tmp.edges[3];
			is >> tmp.formula_number[0];
			is >> tmp.formula_number[1];
			is >> tmp.formula_number[2];
			is >> tmp.formula_number[3];
			boundaries.push_back(tmp);
		}

		return is;
	}

	void BoundaryConditionsSupport::initialize_penalty_parameters()
	{
		mu1 = 1;
	}

	void BoundaryConditionsSupport::calculate_all_boundaries1(MyVector& b)
	{
		int size_b = boundaries1.size();
		for(int i = 0; i < size_b; i++)
			calculate_boundaries1(i, b);
	}

	void BoundaryConditionsSupport::calculate_boundaries1_horizontal(int number, myvector::MyVector & b, Side side)
	{
		Element element = elements[boundaries1[number].elem];
		double hy = get_hy(boundaries1[number].elem);
		double hx = get_hx(boundaries1[number].elem);

		double x0 = nodes[element.nodes[0]].x;
		double y0 = nodes[element.nodes[0]].y;

		double lambda = calculate_lambda(element.number_of_area);

		double n_vec[2], p_etta, p_y;
		int formula;

		if (side == low_edge)
		{
			n_vec[0] = 0; n_vec[1] = -1;
			p_etta = -1;
			p_y = y0;
			formula = boundaries1[number].formula_number[2];
		}
		else
		{
			n_vec[0] = 0; n_vec[1] = 1;
			p_etta = 1;
			p_y = y0 + hy;
			formula = boundaries1[number].formula_number[3];
		}

		double Ugi, Pgi;
		double g_x, g_y;
		double jacobian = hx / 2;
		int n_func_u_all = elements.size() * n_func_u;

		for (int i = 0; i < n_func_u; i++)
		{
			Ugi = 0;
			for (int j = 0; j < 3; j++)
			{
				double  p_ksi = gauss_points_1[j];
				double p_x = p_ksi * hx + x0;
				g_x = gx(formula, p_x, p_y);
				g_y = gy(formula, p_x, p_y);
				Ugi += lambda * gauss_weights_1[j] *
					(g_x * n_vec[0] * dphixksi[i](p_ksi, p_etta) / hx +
						g_x * n_vec[1] * dphixetta[i](p_ksi, p_etta) / hy +
						g_y * n_vec[0] * dphiyksi[i](p_ksi, p_etta) / hx +
						g_y * n_vec[1] * dphiyetta[i](p_ksi, p_etta) / hy) +
					gauss_weights_1[j] * mu1 * (g_x * phix[i](p_ksi, p_etta) +
						g_y * phiy[i](p_ksi, p_etta));
			}
			Ugi *= jacobian;
			b[element.dof_u[i]] += Ugi;
		}

		for (int i = 0; i < n_func_p; i++)
		{
			Pgi = 0;
			for (int j = 0; j < 3; j++)
			{
				double  p_ksi = gauss_points_1[j];
				double p_x = p_ksi * hx + x0;
				g_x = gx(formula, p_x, p_y);
				g_y = gy(formula, p_x, p_y);
				Pgi -= gauss_weights_1[j] * psi[i](p_ksi, p_etta) *
					(g_x * n_vec[0] + g_y * n_vec[1]);
			}
			Pgi *= jacobian;
			b[element.dof_p[i] + n_func_u_all] += Pgi;
		}
	}
	void BoundaryConditionsSupport::calculate_boundaries1_vertical(int number, myvector::MyVector & b, Side side)
	{
		Element element = elements[boundaries1[number].elem];
		double hy = get_hy(boundaries1[number].elem);
		double hx = get_hx(boundaries1[number].elem);

		double x0 = nodes[element.nodes[0]].x;
		double y0 = nodes[element.nodes[0]].y;

		double lambda = calculate_lambda(element.number_of_area);
		double n_vec[2], p_ksi, p_x;
		int formula;

		if(side == left_edge) 
		{
			n_vec[0] = -1; n_vec[1] = 0;
			p_ksi = -1;
			p_x = x0;
			formula = boundaries1[number].formula_number[0];
		}
		else
		{
			n_vec[0] = 1; n_vec[1] = 0;
			p_ksi = 1;
			p_x = x0 + hx;
			formula = boundaries1[number].formula_number[1];
		}
		

		double Ugi, Pgi;

		double g_x, g_y;
		double jacobian = hy / 2;
		int n_func_u_all = elements.size() * n_func_u;

		for (int i = 0; i < n_func_u; i++)
		{
			Ugi = 0;
			for (int j = 0; j < 3; j++)
			{
				double  p_etta = gauss_points_1[j];
				double p_y = p_etta * hy + y0;
				g_x = gx(formula, p_x, p_y);
				g_y = gy(formula, p_x, p_y);
				Ugi += lambda * gauss_weights_1[j] *
					(g_x * n_vec[0] * dphixksi[i](p_ksi, p_etta) / hx +
						g_x * n_vec[1] * dphixetta[i](p_ksi, p_etta) / hy +
						g_y * n_vec[0] * dphiyksi[i](p_ksi, p_etta) / hx +
						g_y * n_vec[1] * dphiyetta[i](p_ksi, p_etta) / hy) +
					gauss_weights_1[j] * mu1 * (g_x * phix[i](p_ksi, p_etta) +
						g_y * phiy[i](p_ksi, p_etta));

			}
			Ugi *= jacobian;
			b[element.dof_u[i]] += Ugi;
		}

		for (int i = 0; i < n_func_p; i++)
		{
			Pgi = 0;
			for (int j = 0; j < 3; j++)
			{
				double  p_etta = gauss_points_1[j];
				double p_y = p_etta * hy + y0;
				g_x = gx(formula, p_x, p_y);
				g_y = gy(formula, p_x, p_y);
				Pgi -= gauss_weights_1[j] * psi[i](p_ksi, p_etta) *
					(g_x * n_vec[0] + g_y * n_vec[1]);

			}
			Pgi *= jacobian;
			b[element.dof_p[i] + n_func_u_all] += Pgi;
		}
	}

	void BoundaryConditionsSupport::calculate_boundaries1(int number, MyVector& b)
	{
		if(boundaries1[number].edges[0] == 1) calculate_boundaries1_vertical(number, b, left_edge);
		if(boundaries1[number].edges[1] == 1) calculate_boundaries1_vertical(number, b, right_edge);
		if(boundaries1[number].edges[2] == 1) calculate_boundaries1_horizontal(number, b, low_edge);
		if(boundaries1[number].edges[3] == 1) calculate_boundaries1_horizontal(number, b, up_edge);
	}

}

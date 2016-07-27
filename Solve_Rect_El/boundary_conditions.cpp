#include "boundary_conditions.h"
#include "element.h"
#include "partition.h"

using namespace std;
using namespace element;
using namespace partition;

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
			is >> tmp.formula_number;
			is >> tmp.edges[0];
			is >> tmp.edges[1];
			is >> tmp.edges[2];
			is >> tmp.edges[3];
			boundaries.push_back(tmp);
		}

		return is;
	}

	void BoundaryConditionsSupport::calculate_all_boundaries1()
	{
		int size_b = boundaries1.size();
		for(int i = 0; i < size_b; i++)
			calculate_boundaries1(i);
	}


	void BoundaryConditionsSupport::calculate_boundaries1_left(int number)
	{
		Element element = elements[boundaries1[number].elem];
		double hy = get_hy(boundaries1[number].elem);
		double hx = get_hx(boundaries1[number].elem);

		double x0 = nodes[element.nodes[0]].x;
		double y0 = nodes[element.nodes[0]].y;

		double lambda = calculate_lambda(element.number_of_area);
		double n_vec[2] = {-1, 0};

		double Ug_vector[4];
		double Pg_vector[4];
		double g_x, g_y;
		double jacobian = hy / 2;
		int n_edges = elements.size() * 4;

		for(int i = 0;  i < 4; i++)
		{
			Ug_vector[i] = 0;
			Pg_vector[i] = 0;
			for(int j = 0; j < 3; j++)
			{
				double  p_etta = gauss_points_1[j];
				double p_y = p_etta * hy + y0;
				g_x = gx(boundaries1[number].formula_number, x0, p_y);
				g_y = gy(boundaries1[number].formula_number, x0, p_y);
				Ug_vector[i] += lambda * gauss_weights_1[j] * 
					(g_x * n_vec[0] * dphixksi[i](-1, p_etta) / hx +
					g_x * n_vec[1] * dphixetta[i](-1, p_etta) / hy +
					g_y * n_vec[0] * dphiyksi[i](-1, p_etta) / hx +
					g_y * n_vec[1] * dphiyetta[i](-1, p_etta) / hy) +
					gauss_weights_1[j] * mu1 * (g_x * phix[i](-1, p_etta) +
					g_y * phiy[i](-1, p_etta));
				Pg_vector[i] -= gauss_weights_1[j] * psi[i](-1, p_etta) *
					(g_x * n_vec[0] + g_y * n_vec[1]);

			}
			Ug_vector[i] *= jacobian;
			Pg_vector[i] *= jacobian;
			A.b[element.edges[i]] += Ug_vector[i];
			A.b[element.nodes[i] + n_edges] += Pg_vector[i];
		}
	}

	void BoundaryConditionsSupport::calculate_boundaries1_right(int number)
	{
		Element element = elements[boundaries1[number].elem];
		double hy = get_hy(boundaries1[number].elem);
		double hx = get_hx(boundaries1[number].elem);

		double x0 = nodes[element.nodes[0]].x;
		double y0 = nodes[element.nodes[0]].y;

		double lambda = calculate_lambda(element.number_of_area);
		double n_vec[2] = {1, 0};

		double Ug_vector[4];
		double Pg_vector[4];
		double g_x, g_y;
		double jacobian = hy / 2;
		int n_edges = elements.size() * 4;

		for(int i = 0;  i < 4; i++)
		{
			Ug_vector[i] = 0;
			Pg_vector[i] = 0;
			for(int j = 0; j < 3; j++)
			{
				double  p_etta = gauss_points_1[j];
				double p_y = p_etta * hy + y0;
				g_x = gx(boundaries1[number].formula_number, x0 + hx, p_y);
				g_y = gy(boundaries1[number].formula_number, x0 + hx, p_y);
				Ug_vector[i] += lambda * gauss_weights_1[j] * 
					(g_x * n_vec[0] * dphixksi[i](1, p_etta) / hx +
					g_x * n_vec[1] * dphixetta[i](1, p_etta) / hy +
					g_y * n_vec[0] * dphiyksi[i](1, p_etta) / hx +
					g_y * n_vec[1] * dphiyetta[i](1, p_etta) / hy) +
					gauss_weights_1[j] * mu1 * (g_x * phix[i](1, p_etta) +
					g_y * phiy[i](1, p_etta));
				Pg_vector[i] -= gauss_weights_1[j] * psi[i](1, p_etta) *
					(g_x * n_vec[0] + g_y * n_vec[1]);			
			}
			Ug_vector[i] *= jacobian;
			Pg_vector[i] *= jacobian;
			A.b[element.edges[i]] += Ug_vector[i];
			A.b[element.nodes[i] + n_edges] += Pg_vector[i];
		}
	}

	void BoundaryConditionsSupport::calculate_boundaries1_low(int number)
	{
		Element element = elements[boundaries1[number].elem];
		double hy = get_hy(boundaries1[number].elem);
		double hx = get_hx(boundaries1[number].elem);

		double x0 = nodes[element.nodes[0]].x;
		double y0 = nodes[element.nodes[0]].y;

		double lambda = calculate_lambda(element.number_of_area);
		double n_vec[2] = {0, -1};

		double Ug_vector[4];
		double Pg_vector[4];
		double g_x, g_y;
		double jacobian = hx / 2;
		int n_edges = elements.size() * 4;

		for(int i = 0;  i < 4; i++)
		{
			Ug_vector[i] = 0;
			Pg_vector[i] = 0;
			for(int j = 0; j < 3; j++)
			{
				double  p_ksi = gauss_points_1[j];
				double p_x = p_ksi * hx + x0;
				g_x = gx(boundaries1[number].formula_number, p_x, y0);
				g_y = gy(boundaries1[number].formula_number, p_x, y0);
				Ug_vector[i] += lambda * gauss_weights_1[j] * 
					(g_x * n_vec[0] * dphixksi[i](p_ksi, -1) / hx +
					g_x * n_vec[1] * dphixetta[i](p_ksi, -1) / hy +
					g_y * n_vec[0] * dphiyksi[i](p_ksi, -1) / hx +
					g_y * n_vec[1] * dphiyetta[i](p_ksi, -1) / hy) +
					gauss_weights_1[j] * mu1 * (g_x * phix[i](p_ksi, -1) +
					g_y * phiy[i](p_ksi, -1));
				Pg_vector[i] -= gauss_weights_1[j] * psi[i](p_ksi, -1) *
					(g_x * n_vec[0] + g_y * n_vec[1]);			
			}
			Ug_vector[i] *= jacobian;
			Pg_vector[i] *= jacobian;
			A.b[element.edges[i]] += Ug_vector[i];
			A.b[element.nodes[i] + n_edges] += Pg_vector[i];
		}
	}

	void BoundaryConditionsSupport::calculate_boundaries1_up(int number)
	{
		Element element = elements[boundaries1[number].elem];
		double hy = get_hy(boundaries1[number].elem);
		double hx = get_hx(boundaries1[number].elem);

		double x0 = nodes[element.nodes[0]].x;
		double y0 = nodes[element.nodes[0]].y;

		double lambda = calculate_lambda(element.number_of_area);
		double n_vec[2] = {0, 1};

		double Ug_vector[4];
		double Pg_vector[4];
		double g_x, g_y;
		double jacobian = hx / 2;
		int n_edges = elements.size() * 4;

		for(int i = 0;  i < 4; i++)
		{
			Ug_vector[i] = 0;
			Pg_vector[i] = 0;
			for(int j = 0; j < 3; j++)
			{
				double  p_ksi = gauss_points_1[j];
				double p_x = p_ksi * hx + x0;
				g_x = gx(boundaries1[number].formula_number, p_x, y0 + hy);
				g_y = gy(boundaries1[number].formula_number, p_x, y0 + hy);
				Ug_vector[i] += lambda * gauss_weights_1[j] * 
					(g_x * n_vec[0] * dphixksi[i](p_ksi, 1) / hx +
					g_x * n_vec[1] * dphixetta[i](p_ksi, 1) / hy +
					g_y * n_vec[0] * dphiyksi[i](p_ksi, 1) / hx +
					g_y * n_vec[1] * dphiyetta[i](p_ksi, 1) / hy) +
					gauss_weights_1[j] * mu1 * (g_x * phix[i](p_ksi, 1) +
					g_y * phiy[i](p_ksi, 1));
				Pg_vector[i] -= gauss_weights_1[j] * psi[i](p_ksi, 1) *
					(g_x * n_vec[0] + g_y * n_vec[1]);			
			}
			Ug_vector[i] *= jacobian;
			Pg_vector[i] *= jacobian;
			A.b[element.edges[i]] += Ug_vector[i];
			A.b[element.nodes[i] + n_edges] += Pg_vector[i];
		}
	}


	void BoundaryConditionsSupport::calculate_boundaries1(int number)
	{
		if(boundaries1[number].edges[0] == 1) calculate_boundaries1_left(number);
		if(boundaries1[number].edges[1] == 1) calculate_boundaries1_right(number);
		if(boundaries1[number].edges[2] == 1) calculate_boundaries1_low(number);
		if(boundaries1[number].edges[3] == 1) calculate_boundaries1_up(number);
	}

}

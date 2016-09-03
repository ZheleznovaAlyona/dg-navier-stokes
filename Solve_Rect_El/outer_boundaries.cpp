#include "boundaries.h"
#include "element.h"
#include "myfunctions.h"

using namespace std;
using namespace element;
using namespace partition;
using namespace parameters;
using namespace basis;
using namespace matrix;

namespace boundaries
{

	void OuterBoundaries::initialize_penalty_parameters()
	{
		sigma = 1;
		mu2 = 1;
	}

	void OuterBoundaries::calculate_outer_boundaries(int element_number, matrix::Matrix& A)
	{
		Element element = elements[element_number];

		int last_node = nodes.size() - 1;
		int left_low_corner_node = element.nodes[0];
		int right_up_corner_node = element.nodes[3];

		if(nodes[left_low_corner_node].x == nodes[0].x)
		{
			calculate_ES_out_left(element_number, A);
			calculate_P_1_out_left(element_number, A);
			calculate_P_2_out_left(element_number, A);
			calculate_SP_out_left(element_number, A);
		}//вертикальная левая граница
		if(nodes[right_up_corner_node].x == nodes[last_node].x)
		{
			calculate_ES_out_right(element_number, A);
			calculate_P_1_out_right(element_number, A);
			calculate_P_2_out_right(element_number, A);
			calculate_SP_out_right(element_number, A);
		}//вертикальная правая граница
		if(nodes[left_low_corner_node].y == nodes[0].y) 
		{
			calculate_ES_out_low(element_number, A);
			calculate_P_1_out_low(element_number, A);
			calculate_P_2_out_low(element_number, A);
			calculate_SP_out_low(element_number, A);
		}//горизонтальная нижняя граница
		if(nodes[right_up_corner_node].y == nodes[last_node].y)	
		{
			calculate_ES_out_up(element_number, A);
			calculate_P_1_out_up(element_number, A);
			calculate_P_2_out_up(element_number, A);
			calculate_SP_out_up(element_number, A);
		}//горизонтальная верхняя граница
	}

	void OuterBoundaries::calculate_ES_out_left(int element_number, Matrix& A)
	{
		vector <vector<double>> S_out, E_out;

		S_out.resize(n_func_u);
		E_out.resize(n_func_u);

		for (int i = 0; i < n_func_u; i++)
		{
			initialize_vector(S_out[i], n_func_u);
			initialize_vector(E_out[i], n_func_u);
		}

		Element element = elements[element_number];
		double lambda = calculate_lambda(element.number_of_area);
		double hx = get_hx(element_number);
		double hy = get_hy(element_number);

		double a =  0.5 * hy * lambda; //якобиан*lambda

		double jacobian = 0.5 * hy;
		double st = jacobian * sigma;

		double n_vec[2] = {-1, 0};
		double n_vec2[2] = {1, 0};

		for(int i = 0; i < n_func_u; i++)
		{
			for(int j = 0; j < n_func_u; j++)
			{
				for(int k = 0; k < 3; k++)
				{
					double p_etta = gauss_points_1[k];
					E_out[i][j] += gauss_weights_1[k] * 
						(phix[j](-1, p_etta) * n_vec[0] * dphixksi[i](-1, p_etta) / hx +
						phix[j](-1, p_etta) * n_vec[1] * dphixetta[i](-1, p_etta) / hy +
						phiy[j](-1, p_etta) * n_vec[0] * dphiyksi[i](-1, p_etta) / hx +
						phiy[j](-1, p_etta) * n_vec[1] * dphiyetta[i](-1, p_etta) / hy -
						phix[i](-1, p_etta) * n_vec[0] * dphixksi[j](-1, p_etta) / hx -
						phix[i](-1, p_etta) * n_vec[1] * dphixetta[j](-1, p_etta) / hy -
						phiy[i](-1, p_etta) * n_vec[0] * dphiyksi[j](-1, p_etta) / hx -
						phiy[i](-1, p_etta) * n_vec[1] * dphiyetta[j](-1, p_etta) / hy);

					S_out[i][j] += gauss_weights_1[k] * 
						(phix[j](-1, p_etta) * phix[i](-1, p_etta) * n_vec2[0] +
						phix[j](-1, p_etta) * phix[i](-1, p_etta) * n_vec2[1] +
						phiy[j](-1, p_etta) * phiy[i](-1, p_etta) * n_vec2[0] +
						phiy[j](-1, p_etta) * phiy[i](-1, p_etta) * n_vec2[1]);
				} 
				E_out[i][j] *= a;
				S_out[i][j] *= st;
			}
		}

		for(int i = 0; i < n_func_u; i++)
		{
			int id_i = element.edges[i];
			for(int j = 0; j < n_func_u; j++)
			{				
				int id_j = element.edges[j];
				A.add_element(id_i, id_j, E_out[i][j] + S_out[i][j]); 
			}
		}
	}
	void OuterBoundaries::calculate_ES_out_right(int element_number, Matrix& A)
	{
		vector <vector<double>> S_out, E_out;

		S_out.resize(n_func_u);
		E_out.resize(n_func_u);

		for (int i = 0; i < n_func_u; i++)
		{
			initialize_vector(S_out[i], n_func_u);
			initialize_vector(E_out[i], n_func_u);
		}

		Element element = elements[element_number];
		double lambda = calculate_lambda(element.number_of_area);
		double hx = get_hx(element_number);
		double hy = get_hy(element_number);

		double a =  0.5 * hy * lambda; //якобиан*lambda

		double jacobian = 0.5 * hy;
		double st = jacobian * sigma;

		double n_vec[2] = {1, 0};
		double n_vec2[2] = {1, 0};

		for(int i = 0; i < n_func_u; i++)
		{
			for(int j = 0; j < n_func_u; j++)
			{
				for(int k = 0; k < 3; k++)
				{
					double p_etta = gauss_points_1[k];
					E_out[i][j] += gauss_weights_1[k] * 
						(phix[j](1, p_etta) * n_vec[0] * dphixksi[i](1, p_etta) / hx +
						phix[j](1, p_etta) * n_vec[1] * dphixetta[i](1, p_etta) / hy +
						phiy[j](1, p_etta) * n_vec[0] * dphiyksi[i](1, p_etta) / hx +
						phiy[j](1, p_etta) * n_vec[1] * dphiyetta[i](1, p_etta) / hy -
						phix[i](1, p_etta) * n_vec[0] * dphixksi[j](1, p_etta) / hx -
						phix[i](1, p_etta) * n_vec[1] * dphixetta[j](1, p_etta) / hy -
						phiy[i](1, p_etta) * n_vec[0] * dphiyksi[j](1, p_etta) / hx -
						phiy[i](1, p_etta) * n_vec[1] * dphiyetta[j](1, p_etta) / hy);

					S_out[i][j] += gauss_weights_1[k] * 
						(phix[j](1, p_etta) * phix[i](1, p_etta) * n_vec2[0] +
						phix[j](1, p_etta) * phix[i](1, p_etta) * n_vec2[1] +
						phiy[j](1, p_etta) * phiy[i](1, p_etta) * n_vec2[0] +
						phiy[j](1, p_etta) * phiy[i](1, p_etta) * n_vec2[1]);
				} 
				E_out[i][j] *= a;
				S_out[i][j] *= st;
			}
		}

		for(int i = 0; i < n_func_u; i++)
		{
			int id_i = element.edges[i];
			for(int j = 0; j < n_func_u; j++)
			{				
				int id_j = element.edges[j];
				A.add_element(id_i, id_j, E_out[i][j] + S_out[i][j]); 
			}
		}
	}
	void OuterBoundaries::calculate_ES_out_low(int element_number, Matrix& A)
	{
		vector <vector<double>> S_out, E_out;

		S_out.resize(n_func_u);
		E_out.resize(n_func_u);

		for (int i = 0; i < n_func_u; i++)
		{
			initialize_vector(S_out[i], n_func_u);
			initialize_vector(E_out[i], n_func_u);
		}

		Element element = elements[element_number];
		double lambda = calculate_lambda(element.number_of_area);
		double hx = get_hx(element_number);
		double hy = get_hy(element_number);

		double a =  0.5 * hx * lambda; //якобиан*lambda

		double jacobian = 0.5 * hx;
		double st = jacobian * sigma;

		double n_vec[2] = {0, -1};
		double n_vec2[2] = {0, 1};

		for(int i = 0; i < n_func_u; i++)
		{
			for(int j = 0; j < n_func_u; j++)
			{
				for(int k = 0; k < 3; k++)
				{
					double p_ksi = gauss_points_1[k];
					E_out[i][j] += gauss_weights_1[k] * 
						(phix[j](p_ksi, -1) * n_vec[0] * dphixksi[i](p_ksi, -1) / hx +
						phix[j](p_ksi, -1) * n_vec[1] * dphixetta[i](p_ksi, -1) / hy +
						phiy[j](p_ksi, -1) * n_vec[0] * dphiyksi[i](p_ksi, -1) / hx +
						phiy[j](p_ksi, -1) * n_vec[1] * dphiyetta[i](p_ksi, -1) / hy -
						phix[i](p_ksi, -1) * n_vec[0] * dphixksi[j](p_ksi, -1) / hx -
						phix[i](p_ksi, -1) * n_vec[1] * dphixetta[j](p_ksi, -1) / hy -
						phiy[i](p_ksi, -1) * n_vec[0] * dphiyksi[j](p_ksi, -1) / hx -
						phiy[i](p_ksi, -1) * n_vec[1] * dphiyetta[j](p_ksi, -1) / hy);

					S_out[i][j] += gauss_weights_1[k] * 
						(phix[j](p_ksi, -1) * phix[i](p_ksi, -1) * n_vec2[0] +
						phix[j](p_ksi, -1) * phix[i](p_ksi, -1) * n_vec2[1] +
						phiy[j](p_ksi, -1) * phiy[i](p_ksi, -1) * n_vec2[0] +
						phiy[j](p_ksi, -1) * phiy[i](p_ksi, -1) * n_vec2[1]);
				} 
				E_out[i][j] *= a;
				S_out[i][j] *= st;
			}
		}

		for(int i = 0; i < n_func_u; i++)
		{
			int id_i = element.edges[i];
			for(int j = 0; j < n_func_u; j++)
			{				
				int id_j = element.edges[j];
				A.add_element(id_i, id_j, E_out[i][j] + S_out[i][j]); 
			}
		}
	}
	void OuterBoundaries::calculate_ES_out_up(int element_number, Matrix& A)
	{
		vector <vector<double>> S_out, E_out;

		S_out.resize(n_func_u);
		E_out.resize(n_func_u);

		for (int i = 0; i < n_func_u; i++)
		{
			initialize_vector(S_out[i], n_func_u);
			initialize_vector(E_out[i], n_func_u);
		}

		Element element = elements[element_number];
		double lambda = calculate_lambda(element.number_of_area);
		double hx = get_hx(element_number);
		double hy = get_hy(element_number);

		double a =  0.5 * hx * lambda; //якобиан*lambda

		double jacobian = 0.5 * hx;
		double st = jacobian * sigma;

		double n_vec[2] = {0, 1};
		double n_vec2[2] = {0, 1};

		for(int i = 0; i < n_func_u; i++)
		{
			for(int j = 0; j < n_func_u; j++)
			{
				for(int k = 0; k < 3; k++)
				{
					double p_ksi = gauss_points_1[k];
					E_out[i][j] += gauss_weights_1[k] * 
						(phix[j](p_ksi, 1) * n_vec[0] * dphixksi[i](p_ksi, 1) / hx +
						phix[j](p_ksi, 1) * n_vec[1] * dphixetta[i](p_ksi, 1) / hy +
						phiy[j](p_ksi, 1) * n_vec[0] * dphiyksi[i](p_ksi, 1) / hx +
						phiy[j](p_ksi, 1) * n_vec[1] * dphiyetta[i](p_ksi, 1) / hy -
						phix[i](p_ksi, 1) * n_vec[0] * dphixksi[j](p_ksi, 1) / hx -
						phix[i](p_ksi, 1) * n_vec[1] * dphixetta[j](p_ksi, 1) / hy -
						phiy[i](p_ksi, 1) * n_vec[0] * dphiyksi[j](p_ksi, 1) / hx -
						phiy[i](p_ksi, 1) * n_vec[1] * dphiyetta[j](p_ksi, 1) / hy);

					S_out[i][j] += gauss_weights_1[k] * 
						(phix[j](p_ksi, 1) * phix[i](p_ksi, 1) * n_vec2[0] +
						phix[j](p_ksi, 1) * phix[i](p_ksi, 1) * n_vec2[1] +
						phiy[j](p_ksi, 1) * phiy[i](p_ksi, 1) * n_vec2[0] +
						phiy[j](p_ksi, 1) * phiy[i](p_ksi, 1) * n_vec2[1]);
				} 
				E_out[i][j] *= a;
				S_out[i][j] *= st;
			}
		}

		for(int i = 0; i < n_func_u; i++)
		{
			int id_i = element.edges[i];
			for(int j = 0; j < n_func_u; j++)
			{				
				int id_j = element.edges[j];
				A.add_element(id_i, id_j, E_out[i][j] + S_out[i][j]); 
			}
		}
	}

	void OuterBoundaries::calculate_P_1_out_left(int element_number, Matrix& A)
	{
		vector <vector<double>> P_1_out;

		P_1_out.resize(n_func_u);

		for (int i = 0; i < n_func_u; i++)
			initialize_vector(P_1_out[i], n_func_p);

		Element element = elements[element_number];
		double rho = calculate_rho(element.number_of_area);
		double hy = get_hy(element_number);

		double a =  0.5 * hy * (1 / rho); //якобиан*1/rho

		double n_vec[2] = {-1, 0};

		for(int i = 0; i < n_func_u; i++)
		{
			for(int j = 0; j < n_func_p; j++)
			{
				for(int k = 0; k < 3; k++)
				{
					double p_etta = gauss_points_1[k];
					P_1_out[i][j] += gauss_weights_1[k]  * psi[j](-1, p_etta) *
						(phix[i](-1, p_etta) * n_vec[0] + phiy[i](-1, p_etta) * n_vec[1]);

				} 
				P_1_out[i][j] *= a;
			}
		}

		int n_edges = elements.size() * 4;

		for(int i = 0; i < n_func_u; i++)
		{
			int id_i = element.edges[i];
			for(int j = 0; j < n_func_p; j++)
			{				
				int id_j = element.nodes[j] + n_edges;
				A.add_element(id_i, id_j, P_1_out[i][j]); 
			}
		}
	}
	void OuterBoundaries::calculate_P_1_out_right(int element_number, Matrix& A)
	{
		vector <vector<double>> P_1_out;

		P_1_out.resize(n_func_u);

		for (int i = 0; i < n_func_u; i++)
			initialize_vector(P_1_out[i], n_func_p);

		Element element = elements[element_number];
		double rho = calculate_rho(element.number_of_area);
		double hy = get_hy(element_number);

		double a =  0.5 * hy * (1 / rho); //якобиан*1/rho

		double n_vec[2] = {1, 0};

		for(int i = 0; i < n_func_u; i++)
		{
			for(int j = 0; j < n_func_p; j++)
			{
				for(int k = 0; k < 3; k++)
				{
					double p_etta = gauss_points_1[k];
					P_1_out[i][j] += gauss_weights_1[k]  * psi[j](1, p_etta) *
						(phix[i](1, p_etta) * n_vec[0] + phiy[i](1, p_etta) * n_vec[1]);

				} 
				P_1_out[i][j] *= a;
			}
		}

		int n_edges = elements.size() * 4;

		for(int i = 0; i < n_func_u; i++)
		{
			int id_i = element.edges[i];
			for(int j = 0; j < n_func_p; j++)
			{				
				int id_j = element.nodes[j] + n_edges;
				A.add_element(id_i, id_j, P_1_out[i][j]); 
			}
		}
	}
	void OuterBoundaries::calculate_P_1_out_low(int element_number, Matrix& A)
	{
		vector <vector<double>> P_1_out;

		P_1_out.resize(n_func_u);

		for (int i = 0; i < n_func_u; i++)
			initialize_vector(P_1_out[i], n_func_p);

		Element element = elements[element_number];
		double rho = calculate_rho(element.number_of_area);
		double hx = get_hx(element_number);

		double a =  0.5 * hx * (1 / rho); //якобиан*1/rho

		double n_vec[2] = {0, -1};

		for(int i = 0; i < n_func_u; i++)
		{
			for(int j = 0; j < n_func_p; j++)
			{
				for(int k = 0; k < 3; k++)
				{
					double p_ksi = gauss_points_1[k];
					P_1_out[i][j] += gauss_weights_1[k]  * psi[j](p_ksi, -1) *
						(phix[i](p_ksi, -1) * n_vec[0] + phiy[i](p_ksi, -1) * n_vec[1]);

				} 
				P_1_out[i][j] *= a;
			}
		}

		int n_edges = elements.size() * 4;

		for(int i = 0; i < n_func_u; i++)
		{
			int id_i = element.edges[i];
			for(int j = 0; j < n_func_p; j++)
			{				
				int id_j = element.nodes[j] + n_edges;
				A.add_element(id_i, id_j, P_1_out[i][j]); 
			}
		}
	}
	void OuterBoundaries::calculate_P_1_out_up(int element_number, Matrix& A)
	{
		vector <vector<double>> P_1_out;

		P_1_out.resize(n_func_u);

		for (int i = 0; i < n_func_u; i++)
			initialize_vector(P_1_out[i], n_func_p);

		Element element = elements[element_number];
		double rho = calculate_rho(element.number_of_area);
		double hx = get_hx(element_number);

		double a =  0.5 * hx * (1 / rho); //якобиан*1/rho

		double n_vec[2] = {0, 1};

		for(int i = 0; i < n_func_u; i++)
		{
			for(int j = 0; j < n_func_p; j++)
			{
				for(int k = 0; k < 3; k++)
				{
					double p_ksi = gauss_points_1[k];
					P_1_out[i][j] += gauss_weights_1[k]  * psi[j](p_ksi, 1) *
						(phix[i](p_ksi, 1) * n_vec[0] + phiy[i](p_ksi, 1) * n_vec[1]);

				} 
				P_1_out[i][j] *= a;
			}
		}

		int n_edges = elements.size() * 4;

		for(int i = 0; i < n_func_u; i++)
		{
			int id_i = element.edges[i];
			for(int j = 0; j < n_func_p; j++)
			{				
				int id_j = element.nodes[j] + n_edges;
				A.add_element(id_i, id_j, P_1_out[i][j]); 
			}
		}
	}

	void OuterBoundaries::calculate_P_2_out_left(int element_number, Matrix& A)
	{
		vector <vector<double>> P_2_out;

		P_2_out.resize(n_func_p);

		for (int i = 0; i < n_func_p; i++)
			initialize_vector(P_2_out[i], n_func_u);

		Element element = elements[element_number];
		double hy = get_hy(element_number);

		double a =  -0.5 * hy; //-якобиан

		double n_vec[2] = {-1, 0};

		for(int i = 0; i < n_func_p; i++)
		{
			for(int j = 0; j < n_func_u; j++)
			{
				for(int k = 0; k < 3; k++)
				{
					double p_etta = gauss_points_1[k];
					P_2_out[i][j] += gauss_weights_1[k]  * psi[i](-1, p_etta) *
						(phix[j](-1, p_etta) * n_vec[0] + phiy[j](-1, p_etta) * n_vec[1]);

				} 
				P_2_out[i][j] *= a;
			}
		}

		int n_edges = elements.size() * 4;

		for(int i = 0; i < n_func_p; i++)
		{
			int id_i = element.nodes[i] + n_edges;
			for(int j = 0; j < n_func_u; j++)
			{				
				int id_j = element.edges[j];
				A.add_element(id_i, id_j, P_2_out[i][j]); 
			}
		}
	}
	void OuterBoundaries::calculate_P_2_out_right(int element_number, Matrix& A)
	{
		vector <vector<double>> P_2_out;

		P_2_out.resize(n_func_p);

		for (int i = 0; i < n_func_p; i++)
			initialize_vector(P_2_out[i], n_func_u);

		Element element = elements[element_number];
		double hy = get_hy(element_number);

		double a =  -0.5 * hy; //-якобиан

		double n_vec[2] = {1, 0};

		for(int i = 0; i < n_func_p; i++)
		{
			for(int j = 0; j < n_func_u; j++)
			{
				for(int k = 0; k < 3; k++)
				{
					double p_etta = gauss_points_1[k];
					P_2_out[i][j] += gauss_weights_1[k]  * psi[i](1, p_etta) *
						(phix[j](1, p_etta) * n_vec[0] + phiy[j](1, p_etta) * n_vec[1]);

				} 
				P_2_out[i][j] *= a;
			}
		}

		int n_edges = elements.size() * 4;

		for(int i = 0; i < n_func_p; i++)
		{
			int id_i = element.nodes[i] + n_edges;
			for(int j = 0; j < n_func_u; j++)
			{				
				int id_j = element.edges[j];
				A.add_element(id_i, id_j, P_2_out[i][j]); 
			}
		}
	}
	void OuterBoundaries::calculate_P_2_out_low(int element_number, Matrix& A)
	{
		vector <vector<double>> P_2_out;

		P_2_out.resize(n_func_p);

		for (int i = 0; i < n_func_p; i++)
			initialize_vector(P_2_out[i], n_func_u);

		Element element = elements[element_number];
		double hx = get_hx(element_number);

		double a =  -0.5 * hx; //-якобиан

		double n_vec[2] = {0, -1};

		for(int i = 0; i < n_func_p; i++)
		{
			for(int j = 0; j < n_func_u; j++)
			{
				for(int k = 0; k < 3; k++)
				{
					double p_ksi = gauss_points_1[k];
					P_2_out[i][j] += gauss_weights_1[k]  * psi[i](p_ksi, -1) *
						(phix[j](p_ksi, -1) * n_vec[0] + phiy[j](p_ksi, -1) * n_vec[1]);

				} 
				P_2_out[i][j] *= a;
			}
		}

		int n_edges = elements.size() * 4;

		for(int i = 0; i < n_func_p; i++)
		{
			int id_i = element.nodes[i] + n_edges;
			for(int j = 0; j < n_func_u; j++)
			{				
				int id_j = element.edges[j];
				A.add_element(id_i, id_j, P_2_out[i][j]); 
			}
		}
	}
	void OuterBoundaries::calculate_P_2_out_up(int element_number, Matrix& A)
	{
		vector <vector<double>> P_2_out;

		P_2_out.resize(n_func_p);

		for (int i = 0; i < n_func_p; i++)
			initialize_vector(P_2_out[i], n_func_u);

		Element element = elements[element_number];
		double hx = get_hx(element_number);

		double a =  -0.5 * hx; //-якобиан

		double n_vec[2] = {0, 1};

		for(int i = 0; i < n_func_p; i++)
		{
			for(int j = 0; j < n_func_u; j++)
			{
				for(int k = 0; k < 3; k++)
				{
					double p_ksi = gauss_points_1[k];
					P_2_out[i][j] += gauss_weights_1[k]  * psi[i](p_ksi, 1) *
						(phix[j](p_ksi, 1) * n_vec[0] + phiy[j](p_ksi, 1) * n_vec[1]);

				} 
				P_2_out[i][j] *= a;
			}
		}

		int n_edges = elements.size() * 4;

		for(int i = 0; i < n_func_p; i++)
		{
			int id_i = element.nodes[i] + n_edges;
			for(int j = 0; j < n_func_u; j++)
			{				
				int id_j = element.edges[j];
				A.add_element(id_i, id_j, P_2_out[i][j]); 
			}
		}
	}

	void OuterBoundaries::calculate_SP_out_left(int element_number, Matrix& A)
	{
		vector <vector<double>> SP_out;

		SP_out.resize(n_func_p);

		for (int i = 0; i < n_func_p; i++)
			initialize_vector(SP_out[i], n_func_p);

		Element element = elements[element_number];
		double hy = get_hy(element_number);

		double n_vec[2] = {0, -1};

		double jacobian = 0.5 * hy;
		double st = jacobian * mu2;

		for(int i = 0; i < n_func_p; i++)
		{
			for(int j = 0; j < n_func_p; j++)
			{
				for(int k = 0; k < 3; k++)
				{
					double p_etta = gauss_points_1[k];
					SP_out[i][j] += gauss_weights_1[k] * psi[i](-1, p_etta) * psi[j](-1, p_etta);

				} 
				SP_out[i][j] *= st;
			}
		}

		int n_edges = elements.size() * 4;

		for(int i = 0; i < n_func_p; i++)
		{
			int id_i = element.nodes[i] + n_edges;
			for(int j = 0; j < n_func_p; j++)
			{				
				int id_j = element.nodes[j] + n_edges;
				A.add_element(id_i, id_j, SP_out[i][j]); 
			}
		}
	}
	void OuterBoundaries::calculate_SP_out_right(int element_number, Matrix& A)
	{
		vector <vector<double>> SP_out;

		SP_out.resize(n_func_p);

		for (int i = 0; i < n_func_p; i++)
			initialize_vector(SP_out[i], n_func_p);

		Element element = elements[element_number];
		double hy = get_hy(element_number);

		double n_vec[2] = {0, 1};

		double jacobian = 0.5 * hy;
		double st = jacobian * mu2;

		for(int i = 0; i < n_func_p; i++)
		{
			for(int j = 0; j < n_func_p; j++)
			{
				for(int k = 0; k < 3; k++)
				{
					double p_etta = gauss_points_1[k];
					SP_out[i][j] += gauss_weights_1[k] * psi[i](1, p_etta) * psi[j](1, p_etta);

				} 
				SP_out[i][j] *= st;
			}
		}

		int n_edges = elements.size() * 4;

		for(int i = 0; i < n_func_p; i++)
		{
			int id_i = element.nodes[i] + n_edges;
			for(int j = 0; j < n_func_p; j++)
			{				
				int id_j = element.nodes[j] + n_edges;
				A.add_element(id_i, id_j, SP_out[i][j]); 
			}
		}
	}
	void OuterBoundaries::calculate_SP_out_low(int element_number, Matrix& A)
	{
		vector <vector<double>> SP_out;

		SP_out.resize(n_func_p);

		for (int i = 0; i < n_func_p; i++)
			initialize_vector(SP_out[i], n_func_p);

		Element element = elements[element_number];
		double hx = get_hx(element_number);

		double n_vec[2] = {-1, 0};

		double jacobian = 0.5 * hx;
		double st = jacobian * mu2;

		for(int i = 0; i < n_func_p; i++)
		{
			for(int j = 0; j < n_func_p; j++)
			{
				for(int k = 0; k < 3; k++)
				{
					double p_ksi = gauss_points_1[k];
					SP_out[i][j] += gauss_weights_1[k] * psi[i](p_ksi, -1) * psi[j](p_ksi, -1);

				} 
				SP_out[i][j] *= st;
			}
		}

		int n_edges = elements.size() * 4;

		for(int i = 0; i < n_func_p; i++)
		{
			int id_i = element.nodes[i] + n_edges;
			for(int j = 0; j < n_func_p; j++)
			{				
				int id_j = element.nodes[j] + n_edges;
				A.add_element(id_i, id_j, SP_out[i][j]); 
			}
		}
	}
	void OuterBoundaries::calculate_SP_out_up(int element_number, Matrix& A)
	{
		vector <vector<double>> SP_out;

		SP_out.resize(n_func_p);

		for (int i = 0; i < n_func_p; i++)
			initialize_vector(SP_out[i], n_func_p);

		Element element = elements[element_number];
		double hx = get_hx(element_number);

		double n_vec[2] = {1, 0};

		double jacobian = 0.5 * hx;
		double st = jacobian * mu2;

		for(int i = 0; i < n_func_p; i++)
		{
			for(int j = 0; j < n_func_p; j++)
			{
				for(int k = 0; k < 3; k++)
				{
					double p_ksi = gauss_points_1[k];
					SP_out[i][j] += gauss_weights_1[k] * psi[i](p_ksi, 1) * psi[j](p_ksi, 1);

				} 
				SP_out[i][j] *= st;
			}
		}

		int n_edges = elements.size() * 4;

		for(int i = 0; i < n_func_p; i++)
		{
			int id_i = element.nodes[i] + n_edges;
			for(int j = 0; j < n_func_p; j++)
			{				
				int id_j = element.nodes[j] + n_edges;
				A.add_element(id_i, id_j, SP_out[i][j]); 
			}
		}
	}
}
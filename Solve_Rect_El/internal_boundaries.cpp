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
	void InternalBoundaries::initialize_penalty_parameters()
	{
		sigma = 1;
		mu2 = 1;
	}

	void InternalBoundaries::calculate_internal_boundaries(int element_number, Matrix& A)
	{
		Element element = elements[element_number];
		int neighbor_element;

		for(int k = 0; k < 4; k++)
		{
			neighbor_element = element.neighbors[k];
			//существующий соседний элемент с бОльшим номером
			if(neighbor_element > element_number)
			{
				//левый/правый сосед->вертикальная граница
				if(k == 0 || k == 1)
				{
					calculate_ES_vertical(element_number, neighbor_element, A);
					calculate_P_1_vertical(element_number, neighbor_element, A);
					calculate_P_2_vertical(element_number, neighbor_element, A);
					calculate_SP_vertical(element_number, neighbor_element, A);
				}
				else
				{
					calculate_ES_horizontal(element_number, neighbor_element, A);
					calculate_P_1_horizontal(element_number, neighbor_element, A);
					calculate_P_2_horizontal(element_number, neighbor_element, A);
					calculate_SP_horizontal(element_number, neighbor_element, A);
				}
			}
		}
	}

	void InternalBoundaries::calculate_ES_horizontal(int element_number1, int element_number2, matrix::Matrix& A)
	{
		vector <vector<double>> AK, AN, BK, BN, SNN, SNK, SKN, SKK, ES;

		AK.resize(n_func_u); 
		AN.resize(n_func_u);
		BK.resize(n_func_u);
		BN.resize(n_func_u);
		SNN.resize(n_func_u);
		SNK.resize(n_func_u);
		SKN.resize(n_func_u);
		SKK.resize(n_func_u);
		ES.resize(n_func_u * 2);

		for (int i = 0; i < n_func_u; i++)
		{
			initialize_vector(AK[i], n_func_u);
			initialize_vector(AN[i], n_func_u);
			initialize_vector(BK[i], n_func_u);
			initialize_vector(BN[i], n_func_u);
			initialize_vector(SNN[i], n_func_u);
			initialize_vector(SNK[i], n_func_u);
			initialize_vector(SKN[i], n_func_u);
			initialize_vector(SKK[i], n_func_u);
			initialize_vector(ES[i], n_func_u * 2);
			initialize_vector(ES[i + n_func_u], n_func_u * 2);
		}
			
		Element element = elements[element_number1];
		Element element_2 = elements[element_number2];
		double lambda = calculate_lambda(element.number_of_area);
		double lambda_2 = calculate_lambda(element_2.number_of_area);
		double hx = get_hx(element_number1);
		double hy = get_hy(element_number1);
		double hy_2 = get_hy(element_number2);

		double a1 =  0.25 * hx * lambda; //якобиан*0.5*lambda
		double a2 = 0.25 * hx * lambda_2;

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
					AN[i][j] += gauss_weights_1[k] * 
						(phix[j](p_ksi, 1) * n_vec[0] * dphixksi[i](p_ksi, 1) / hx +
						phix[j](p_ksi, 1) * n_vec[1] * dphixetta[i](p_ksi, 1) / hy +
						phiy[j](p_ksi, 1) * n_vec[0] * dphiyksi[i](p_ksi, 1) / hx +
						phiy[j](p_ksi, 1) * n_vec[1] * dphiyetta[i](p_ksi, 1) / hy -
						phix[i](p_ksi, 1) * n_vec[0] * dphixksi[j](p_ksi, 1) / hx -
						phix[i](p_ksi, 1) * n_vec[1] * dphixetta[j](p_ksi, 1) / hy -
						phiy[i](p_ksi, 1) * n_vec[0] * dphiyksi[j](p_ksi, 1) / hx -
						phiy[i](p_ksi, 1) * n_vec[1] * dphiyetta[j](p_ksi, 1) / hy);
					AK[i][j] += gauss_weights_1[k] * 
						(phix[j](p_ksi, -1) * n_vec[0] * dphixksi[i](p_ksi, -1) / hx +
						phix[j](p_ksi, -1) * n_vec[1] * dphixetta[i](p_ksi, -1) / hy_2 +
						phiy[j](p_ksi, -1) * n_vec[0] * dphiyksi[i](p_ksi, -1) / hx +
						phiy[j](p_ksi, -1) * n_vec[1] * dphiyetta[i](p_ksi, -1) / hy_2 -
						phix[i](p_ksi, -1) * n_vec[0] * dphixksi[j](p_ksi, -1) / hx -
						phix[i](p_ksi, -1) * n_vec[1] * dphixetta[j](p_ksi, -1) / hy_2 -
						phiy[i](p_ksi, -1) * n_vec[0] * dphiyksi[j](p_ksi, -1) / hx -
						phiy[i](p_ksi, -1) * n_vec[1] * dphiyetta[j](p_ksi, -1) / hy_2);
					BN[i][j] += gauss_weights_1[k] *
						(phix[j](p_ksi, -1) * n_vec[0] * dphixksi[i](p_ksi, 1) / hx +
						phix[j](p_ksi, -1) * n_vec[1] * dphixetta[i](p_ksi, 1) / hy +
						phiy[j](p_ksi, -1) * n_vec[0] * dphiyksi[i](p_ksi, 1) / hx +
						phiy[j](p_ksi, -1) * n_vec[1] * dphiyetta[i](p_ksi, 1) / hy -
						phix[i](p_ksi, 1) * n_vec[0] * dphixksi[j](p_ksi, -1) / hx -
						phix[i](p_ksi, 1) * n_vec[1] * dphixetta[j](p_ksi, -1) / hy_2 -
						phiy[i](p_ksi, 1) * n_vec[0] * dphiyksi[j](p_ksi, -1) / hx -
						phiy[i](p_ksi, 1) * n_vec[1] * dphiyetta[j](p_ksi, -1) / hy_2);
					BK[i][j] += gauss_weights_1[k] *
						(phix[j](p_ksi, 1) * n_vec[0] * dphixksi[i](p_ksi, -1) / hx +
						phix[j](p_ksi, 1) * n_vec[1] * dphixetta[i](p_ksi, -1) / hy_2 +
						phiy[j](p_ksi, 1) * n_vec[0] * dphiyksi[i](p_ksi, -1) / hx +
						phiy[j](p_ksi, 1) * n_vec[1] * dphiyetta[i](p_ksi, -1) / hy_2 -
						phix[i](p_ksi, -1) * n_vec[0] * dphixksi[j](p_ksi, 1) / hx -
						phix[i](p_ksi, -1) * n_vec[1] * dphixetta[j](p_ksi, 1) / hy -
						phiy[i](p_ksi, -1) * n_vec[0] * dphiyksi[j](p_ksi, 1) / hx -
						phiy[i](p_ksi, -1) * n_vec[1] * dphiyetta[j](p_ksi, 1) / hy);
					SNN[i][j] += gauss_weights_1[k] * 
						(phix[j](p_ksi, 1) * phix[i](p_ksi, 1) * n_vec2[0] +
						phix[j](p_ksi, 1) * phix[i](p_ksi, 1) * n_vec2[1] +
						phiy[j](p_ksi, 1) * phiy[i](p_ksi, 1) * n_vec2[0] +
						phiy[j](p_ksi, 1) * phiy[i](p_ksi, 1) * n_vec2[1]);
					SNK[i][j] += gauss_weights_1[k] * 
						(phix[j](p_ksi, -1) * phix[i](p_ksi, 1) * n_vec2[0] +
						phix[j](p_ksi, -1) * phix[i](p_ksi, 1) * n_vec2[1] +
						phiy[j](p_ksi, -1) * phiy[i](p_ksi, 1) * n_vec2[0] +
						phiy[j](p_ksi, -1) * phiy[i](p_ksi, 1) * n_vec2[1]);
					SKK[i][j] += gauss_weights_1[k] * 
						(phix[j](p_ksi, -1) * phix[i](p_ksi, -1) * n_vec2[0] +
						phix[j](p_ksi, -1) * phix[i](p_ksi, -1) * n_vec2[1] +
						phiy[j](p_ksi, -1) * phiy[i](p_ksi, -1) * n_vec2[0] +
						phiy[j](p_ksi, -1) * phiy[i](p_ksi, -1) * n_vec2[1]);
					SKN[i][j] += gauss_weights_1[k] * 
						(phix[j](p_ksi, 1) * phix[i](p_ksi, -1) * n_vec2[0] +
						phix[j](p_ksi, 1) * phix[i](p_ksi, -1) * n_vec2[1] +
						phiy[j](p_ksi, 1) * phiy[i](p_ksi, -1) * n_vec2[0] +
						phiy[j](p_ksi, 1) * phiy[i](p_ksi, -1) * n_vec2[1]);

				} 
				AN[i][j] *= a1;
				AK[i][j] *= -a2;
				BN[i][j] *= -a1;
				BK[i][j] *= a2;
				SNN[i][j] *= st;
				SNK[i][j] *= -st;
				SKK[i][j] *= st;
				SKN[i][j] *= -st;
			}
		}

		for(int i = 0; i < n_func_u; i++)
			for(int j = 0; j < n_func_u; j++)
			{
				ES[i][j] = AN[i][j] + SNN[i][j];
				ES[i + n_func_u][j + n_func_u] = AK[i][j] + SKK[i][j];
				ES[i][j + n_func_u] = BN[i][j] + SNK[i][j];
				ES[i + n_func_u][j] = BK[i][j] + SKN[i][j];
			}

		add_ES_to_global(element_number1, element_number2, A, ES);
	}
	void InternalBoundaries::calculate_ES_vertical(int element_number1, int element_number2, matrix::Matrix& A)
	{
		vector <vector<double>> AK, AN, BK, BN, SNN, SNK, SKN, SKK, ES;

		AK.resize(n_func_u);
		AN.resize(n_func_u);
		BK.resize(n_func_u);
		BN.resize(n_func_u);
		SNN.resize(n_func_u);
		SNK.resize(n_func_u);
		SKN.resize(n_func_u);
		SKK.resize(n_func_u);
		ES.resize(n_func_u * 2);

		for (int i = 0; i < n_func_u; i++)
		{
			initialize_vector(AK[i], n_func_u);
			initialize_vector(AN[i], n_func_u);
			initialize_vector(BK[i], n_func_u);
			initialize_vector(BN[i], n_func_u);
			initialize_vector(SNN[i], n_func_u);
			initialize_vector(SNK[i], n_func_u);
			initialize_vector(SKN[i], n_func_u);
			initialize_vector(SKK[i], n_func_u);
			initialize_vector(ES[i], n_func_u * 2);
			initialize_vector(ES[i + n_func_u], n_func_u * 2);
		}

		Element element = elements[element_number1];
		Element element_2 = elements[element_number2];
		double lambda = calculate_lambda(element.number_of_area);
		double lambda_2 = calculate_lambda(element_2.number_of_area);
		double hx = get_hx(element_number1);
		double hy = get_hy(element_number1);
		double hx_2 = get_hx(element_number2);

		double a1 =  0.25 * hy * lambda; //якобиан*0.5*lambda
		double a2 = 0.25 * hy * lambda_2;

		double jacobian = 0.5 * hy;
		double st = jacobian * sigma;

		double n_vec[2] = {1, 0};
		double n_vec2[2] = {1, 0};

		for(int i = 0; i < 4; i++)
		{
			for(int j = 0; j < 4; j++)
			{
				for(int k = 0; k < 3; k++)
				{
					double p_etta = gauss_points_1[k];
					AN[i][j] += gauss_weights_1[k] * 
						(phix[j](1, p_etta) * n_vec[0] * dphixksi[i](1, p_etta) / hx +
						phix[j](1, p_etta) * n_vec[1] * dphixetta[i](1, p_etta) / hy +
						phiy[j](1, p_etta) * n_vec[0] * dphiyksi[i](1, p_etta) / hx +
						phiy[j](1, p_etta) * n_vec[1] * dphiyetta[i](1, p_etta) / hy -
						phix[i](1, p_etta) * n_vec[0] * dphixksi[j](1, p_etta) / hx -
						phix[i](1, p_etta) * n_vec[1] * dphixetta[j](1, p_etta) / hy -
						phiy[i](1, p_etta) * n_vec[0] * dphiyksi[j](1, p_etta) / hx -
						phiy[i](1, p_etta) * n_vec[1] * dphiyetta[j](1, p_etta) / hy);
					AK[i][j] += gauss_weights_1[k] * 
						(phix[j](-1, p_etta) * n_vec[0] * dphixksi[i](-1, p_etta) / hx_2 +
						phix[j](-1, p_etta) * n_vec[1] * dphixetta[i](-1, p_etta) / hy +
						phiy[j](-1, p_etta) * n_vec[0] * dphiyksi[i](-1, p_etta) / hx_2 +
						phiy[j](-1, p_etta) * n_vec[1] * dphiyetta[i](-1, p_etta) / hy -
						phix[i](-1, p_etta) * n_vec[0] * dphixksi[j](-1, p_etta) / hx_2 -
						phix[i](-1, p_etta) * n_vec[1] * dphixetta[j](-1, p_etta) / hy -
						phiy[i](-1, p_etta) * n_vec[0] * dphiyksi[j](-1, p_etta) / hx_2 -
						phiy[i](-1, p_etta) * n_vec[1] * dphiyetta[j](-1, p_etta) / hy);
					BN[i][j] += gauss_weights_1[k] *
						(phix[j](-1, p_etta) * n_vec[0] * dphixksi[i](1, p_etta) / hx +
						phix[j](-1, p_etta) * n_vec[1] * dphixetta[i](1, p_etta) / hy +
						phiy[j](-1, p_etta) * n_vec[0] * dphiyksi[i](1, p_etta) / hx +
						phiy[j](-1, p_etta) * n_vec[1] * dphiyetta[i](1, p_etta) / hy -
						phix[i](1, p_etta) * n_vec[0] * dphixksi[j](-1, p_etta) / hx_2 -
						phix[i](1, p_etta) * n_vec[1] * dphixetta[j](-1, p_etta) / hy -
						phiy[i](1, p_etta) * n_vec[0] * dphiyksi[j](-1, p_etta) / hx_2 -
						phiy[i](1, p_etta) * n_vec[1] * dphiyetta[j](-1, p_etta) / hy);
					BK[i][j] += gauss_weights_1[k] *
						(phix[j](1, p_etta) * n_vec[0] * dphixksi[i](-1, p_etta) / hx_2 +
						phix[j](1, p_etta) * n_vec[1] * dphixetta[i](-1, p_etta) / hy +
						phiy[j](1, p_etta) * n_vec[0] * dphiyksi[i](-1, p_etta) / hx_2 +
						phiy[j](1, p_etta) * n_vec[1] * dphiyetta[i](-1, p_etta) / hy -
						phix[i](-1, p_etta) * n_vec[0] * dphixksi[j](1, p_etta) / hx -
						phix[i](-1, p_etta) * n_vec[1] * dphixetta[j](1, p_etta) / hy -
						phiy[i](-1, p_etta) * n_vec[0] * dphiyksi[j](1, p_etta) / hx -
						phiy[i](-1, p_etta) * n_vec[1] * dphiyetta[j](1, p_etta) / hy);
					SNN[i][j] += gauss_weights_1[k] * 
						(phix[j](1, p_etta) * phix[i](1, p_etta) * n_vec2[0] +
						phix[j](1, p_etta) * phix[i](1, p_etta) * n_vec2[1] +
						phiy[j](1, p_etta) * phiy[i](1, p_etta) * n_vec2[0] +
						phiy[j](1, p_etta) * phiy[i](1, p_etta) * n_vec2[1]);
					SNK[i][j] += gauss_weights_1[k] * 
						(phix[j](-1, p_etta) * phix[i](1, p_etta) * n_vec2[0] +
						phix[j](-1, p_etta) * phix[i](1, p_etta) * n_vec2[1] +
						phiy[j](-1, p_etta) * phiy[i](1, p_etta) * n_vec2[0] +
						phiy[j](-1, p_etta) * phiy[i](1, p_etta) * n_vec2[1]);
					SKK[i][j] += gauss_weights_1[k] * 
						(phix[j](-1, p_etta) * phix[i](-1, p_etta) * n_vec2[0] +
						phix[j](-1, p_etta) * phix[i](-1, p_etta) * n_vec2[1] +
						phiy[j](-1, p_etta) * phiy[i](-1, p_etta) * n_vec2[0] +
						phiy[j](-1, p_etta) * phiy[i](-1, p_etta) * n_vec2[1]);
					SKN[i][j] += gauss_weights_1[k] * 
						(phix[j](1, p_etta) * phix[i](-1, p_etta) * n_vec2[0] +
						phix[j](1, p_etta) * phix[i](-1, p_etta) * n_vec2[1] +
						phiy[j](1, p_etta) * phiy[i](-1, p_etta) * n_vec2[0] +
						phiy[j](1, p_etta) * phiy[i](-1, p_etta) * n_vec2[1]);

				} 
				AN[i][j] *= a1;
				AK[i][j] *= -a2;
				BN[i][j] *= -a1;
				BK[i][j] *= a2;
				SNN[i][j] *= st;
				SNK[i][j] *= -st;
				SKK[i][j] *= st;
				SKN[i][j] *= -st;
			}
		}

		for(int i = 0; i < n_func_u; i++)
			for(int j = 0; j < n_func_u; j++)
			{
				ES[i][j] = AN[i][j] + SNN[i][j];
				ES[i + n_func_u][j + n_func_u] = AK[i][j] + SKK[i][j];
				ES[i][j + n_func_u] = BN[i][j] + SNK[i][j];
				ES[i + n_func_u][j] = BK[i][j] + SKN[i][j];
			}

		add_ES_to_global(element_number1, element_number2, A, ES);
	}
	void InternalBoundaries::add_ES_to_global(int element_number, int neighbor_element_number, Matrix& A, vector <vector<double>>& ES)
	{
		int id_i, id_j;
		Element element = elements[element_number];
		Element neighbor_element = elements[neighbor_element_number];

		for(int i = 0; i < n_func_u; i++)
		{
			id_i = element.edges[i];
			for(int j = 0; j < n_func_u; j++)
			{				
				id_j = element.edges[j];
				A.add_element(id_i, id_j, ES[i][j]); 
			}

			for(int j = n_func_u; j < n_func_u * 2; j++)
			{				
				id_j = neighbor_element.edges[j - n_func_u];
				A.add_element(id_i, id_j, ES[i][j]); 
			}
		}

		for(int i = n_func_u; i < n_func_u * 2; i++)
		{
			id_i = neighbor_element.edges[i - n_func_u];
			for(int j = 0; j < n_func_u; j++)
			{				
				id_j = element.edges[j];
				A.add_element(id_i, id_j, ES[i][j]); 
			}

			for(int j = n_func_u; j < n_func_u * 2; j++)
			{				
				id_j = neighbor_element.edges[j - n_func_u];
				A.add_element(id_i, id_j, ES[i][j]); 
			}
		}
	}

	void InternalBoundaries::calculate_P_1_horizontal(int element_number1, int element_number2, matrix::Matrix& A)
	{
		vector <vector<double>> AK, AN, BK, BN, P_1;

		AK.resize(n_func_u);
		AN.resize(n_func_u);
		BK.resize(n_func_u);
		BN.resize(n_func_u);
		P_1.resize(n_func_u * 2);

		for (int i = 0; i < n_func_u; i++)
		{
			initialize_vector(AK[i], n_func_p);
			initialize_vector(AN[i], n_func_p);
			initialize_vector(BK[i], n_func_p);
			initialize_vector(BN[i], n_func_p);
			initialize_vector(P_1[i], n_func_p * 2);
			initialize_vector(P_1[i + n_func_u], n_func_p * 2);
		}

		Element element = elements[element_number1];
		Element element_2 = elements[element_number2];
		double rho = calculate_rho(element.number_of_area);
		double rho_2 = calculate_rho(element_2.number_of_area);
		double hx = get_hx(element_number1);

		double a1 =  0.25 * hx * (1 / rho); //якобиан*0.5*(1/rho)
		double a2 = 0.25 * hx * (1 / rho_2);

		double n_vec[2] = {0, 1};

		for(int i = 0; i < n_func_u; i++)
		{
			for(int j = 0; j < n_func_p; j++)
			{
				for(int k = 0; k < 3; k++)
				{
					double p_ksi = gauss_points_1[k];
					AN[i][j] += gauss_weights_1[k] * psi[j](p_ksi, 1) *
						(phix[i](p_ksi, 1) * n_vec[0] + phiy[i](p_ksi, 1) * n_vec[1]);
					AK[i][j] += gauss_weights_1[k] * psi[j](p_ksi, -1) *
						(phix[i](p_ksi, -1) * n_vec[0] + phiy[i](p_ksi, -1) * n_vec[1]);
					BN[i][j] += gauss_weights_1[k] * psi[j](p_ksi, -1) *
						(phix[i](p_ksi, 1) * n_vec[0] + phiy[i](p_ksi, 1) * n_vec[1]);
					BK[i][j] += gauss_weights_1[k] * psi[j](p_ksi, 1) *
						(phix[i](p_ksi, -1) * n_vec[0] + phiy[i](p_ksi, -1) * n_vec[1]);

				} 
				AN[i][j] *= a1;
				AK[i][j] *= -a2;
				BN[i][j] *= a1;
				BK[i][j] *= -a2;
			}
		}

		for(int i = 0; i < n_func_u; i++)
			for(int j = 0; j < n_func_p; j++)
			{
				P_1[i][j] = AN[i][j];
				P_1[i + n_func_u][j + n_func_p] = AK[i][j];
				P_1[i][j + n_func_p] = BN[i][j];
				P_1[i + n_func_u][j] = BK[i][j];
			}

		add_P_1_to_global(element_number1, element_number2, A, P_1);
	}
	void InternalBoundaries::calculate_P_1_vertical(int element_number1, int element_number2, matrix::Matrix& A)
	{
		vector <vector<double>> AK, AN, BK, BN, P_1;

		AK.resize(n_func_u);
		AN.resize(n_func_u);
		BK.resize(n_func_u);
		BN.resize(n_func_u);
		P_1.resize(n_func_u * 2);

		for (int i = 0; i < n_func_u; i++)
		{
			initialize_vector(AK[i], n_func_p);
			initialize_vector(AN[i], n_func_p);
			initialize_vector(BK[i], n_func_p);
			initialize_vector(BN[i], n_func_p);
			initialize_vector(P_1[i], n_func_p * 2);
			initialize_vector(P_1[i + n_func_u], n_func_p * 2);
		}

		Element element = elements[element_number1];
		Element element_2 = elements[element_number2];
		double hy = get_hy(element_number1);

		double a1 =  0.25 * hy; //якобиан*0.5*lambda
		double a2 = 0.25 * hy;

		double n_vec[2] = {1, 0};

		for(int i = 0; i < n_func_u; i++)
		{
			for(int j = 0; j < n_func_p; j++)
			{
				for(int k = 0; k < 3; k++)
				{
					double p_etta = gauss_points_1[k];
					AN[i][j] += gauss_weights_1[k] * psi[j](1, p_etta) *
						(phix[i](1, p_etta) * n_vec[0] + phiy[i](1, p_etta) * n_vec[1]);
					AK[i][j] += gauss_weights_1[k] * psi[j](-1, p_etta) *
						(phix[i](-1, p_etta) * n_vec[0] + phiy[i](-1, p_etta) * n_vec[1]);
					BN[i][j] += gauss_weights_1[k] * psi[j](-1, p_etta) *
						(phix[i](1, p_etta) * n_vec[0] + phiy[i](1, p_etta) * n_vec[1]);
					BK[i][j] += gauss_weights_1[k] * psi[j](1, p_etta) *
						(phix[i](-1, p_etta) * n_vec[0] + phiy[i](-1, p_etta) * n_vec[1]);

				} 
				AN[i][j] *= a1;
				AK[i][j] *= -a2;
				BN[i][j] *= -a1;
				BK[i][j] *= a2;
			}
		}

		for(int i = 0; i < n_func_u; i++)
			for(int j = 0; j < n_func_p; j++)
			{
				P_1[i][j] = AN[i][j];
				P_1[i + n_func_u][j + n_func_p] = AK[i][j];
				P_1[i][j + n_func_p] = BN[i][j];
				P_1[i + n_func_u][j] = BK[i][j];
			}

		add_P_1_to_global(element_number1, element_number2, A, P_1);
	}
	void InternalBoundaries::add_P_1_to_global(int element_number, int neighbor_element_number, Matrix& A, vector <vector<double>>& P_1)
	{
		int id_i, id_j;
		Element element = elements[element_number];
		Element neighbor_element = elements[neighbor_element_number];
		int n_edges = elements.size() * 4;

		for(int i = 0; i < n_func_u; i++)
		{
			id_i = element.edges[i];
			for(int j = 0; j < n_func_p; j++)
			{				
				id_j = element.nodes[j] + n_edges;
				A.add_element(id_i, id_j, P_1[i][j]); 
			}

			for(int j = n_func_p; j < n_func_p * 2; j++)
			{				
				id_j = neighbor_element.nodes[j - n_func_p] + n_edges;
				A.add_element(id_i, id_j, P_1[i][j]); 
			}
		}

		for(int i = n_func_u; i <  n_func_u * 2; i++)
		{
			id_i = neighbor_element.edges[i - n_func_u];
			for(int j = 0; j < n_func_p; j++)
			{				
				id_j = element.nodes[j] + n_edges;
				A.add_element(id_i, id_j, P_1[i][j]); 
			}

			for(int j = n_func_p; j < n_func_p * 2; j++)
			{				
				id_j = neighbor_element.nodes[j - n_func_p] + n_edges;
				A.add_element(id_i, id_j, P_1[i][j]); 
			}
		}
	}

	void InternalBoundaries::calculate_P_2_horizontal(int element_number1, int element_number2, matrix::Matrix& A)
	{
		vector <vector<double>> AK, AN, BK, BN, P_2;

		AK.resize(n_func_p);
		AN.resize(n_func_p);
		BK.resize(n_func_p);
		BN.resize(n_func_p);
		P_2.resize(n_func_p * 2);

		for (int i = 0; i < n_func_p; i++)
		{
			initialize_vector(AK[i], n_func_u);
			initialize_vector(AN[i], n_func_u);
			initialize_vector(BK[i], n_func_u);
			initialize_vector(BN[i], n_func_u);
			initialize_vector(P_2[i], n_func_u * 2);
			initialize_vector(P_2[i + n_func_p], n_func_u * 2);
		}

		Element element = elements[element_number1];
		Element element_2 = elements[element_number2];
		double hx = get_hx(element_number1);

		double a1 =  0.25 * hx; //якобиан*0.5
		double a2 = 0.25 * hx;

		double n_vec[2] = {0, 1};

		for(int i = 0; i < n_func_p; i++)
		{
			for(int j = 0; j < n_func_u; j++)
			{
				for(int k = 0; k < 3; k++)
				{
					double p_ksi = gauss_points_1[k];
					AN[i][j] += gauss_weights_1[k] * psi[i](p_ksi, 1) *
						(phix[j](p_ksi, 1) * n_vec[0] + phiy[j](p_ksi, 1) * n_vec[1]);
					AK[i][j] += gauss_weights_1[k] * psi[i](p_ksi, -1) *
						(phix[j](p_ksi, -1) * n_vec[0] + phiy[j](p_ksi, -1) * n_vec[1]);
					BN[i][j] += gauss_weights_1[k] * psi[i](p_ksi, 1) *
						(phix[j](p_ksi, -1) * n_vec[0] + phiy[j](p_ksi, -1) * n_vec[1]);
					BK[i][j] += gauss_weights_1[k] * psi[i](p_ksi, -1) *
						(phix[j](p_ksi, 1) * n_vec[0] + phiy[j](p_ksi, 1) * n_vec[1]);

				} 
				AN[i][j] *= -a1;
				AK[i][j] *= a2;
				BN[i][j] *= -a1;
				BK[i][j] *= a2;
			}
		}

		for(int i = 0; i < n_func_p; i++)
			for(int j = 0; j < n_func_u; j++)
			{
				P_2[i][j] = AN[i][j];
				P_2[i + n_func_p][j + n_func_u] = AK[i][j];
				P_2[i][j + n_func_u] = BN[i][j];
				P_2[i + n_func_p][j] = BK[i][j];
			}

		add_P_2_to_global(element_number1, element_number2, A, P_2);
	}
	void InternalBoundaries::calculate_P_2_vertical(int element_number1, int element_number2, matrix::Matrix& A)
	{
		vector <vector<double>> AK, AN, BK, BN, P_2;

		AK.resize(n_func_p);
		AN.resize(n_func_p);
		BK.resize(n_func_p);
		BN.resize(n_func_p);
		P_2.resize(n_func_p * 2);

		for (int i = 0; i < n_func_p; i++)
		{
			initialize_vector(AK[i], n_func_u);
			initialize_vector(AN[i], n_func_u);
			initialize_vector(BK[i], n_func_u);
			initialize_vector(BN[i], n_func_u);
			initialize_vector(P_2[i], n_func_u * 2);
			initialize_vector(P_2[i + n_func_p], n_func_u * 2);
		}

		Element element = elements[element_number1];
		Element element_2 = elements[element_number2];
		double hy = get_hy(element_number1);

		double a1 =  0.25 * hy; //якобиан*0.5
		double a2 = 0.25 * hy;

		double n_vec[2] = {1, 0};

		for(int i = 0; i < n_func_p; i++)
		{
			for(int j = 0; j < n_func_u; j++)
			{
				for(int k = 0; k < 3; k++)
				{
					double p_etta = gauss_points_1[k];
					AN[i][j] += gauss_weights_1[k] * psi[i](1, p_etta) *
						(phix[j](1, p_etta) * n_vec[0] + phiy[j](1, p_etta) * n_vec[1]);
					AK[i][j] += gauss_weights_1[k] * psi[i](-1, p_etta) *
						(phix[j](-1, p_etta) * n_vec[0] + phiy[j](-1, p_etta) * n_vec[1]);
					BN[i][j] += gauss_weights_1[k] * psi[i](1, p_etta) *
						(phix[j](-1, p_etta) * n_vec[0] + phiy[j](-1, p_etta) * n_vec[1]);
					BK[i][j] += gauss_weights_1[k] * psi[i](-1, p_etta) *
						(phix[j](1, p_etta) * n_vec[0] + phiy[j](1, p_etta) * n_vec[1]);

				} 
				AN[i][j] *= -a1;
				AK[i][j] *= a2;
				BN[i][j] *= -a1;
				BK[i][j] *= a2;
			}
		}

		for(int i = 0; i < n_func_p; i++)
			for(int j = 0; j < n_func_u; j++)
			{
				P_2[i][j] = AN[i][j];
				P_2[i + n_func_p][j + n_func_u] = AK[i][j];
				P_2[i][j + n_func_u] = BN[i][j];
				P_2[i + n_func_p][j] = BK[i][j];
			}

		add_P_2_to_global(element_number1, element_number2, A, P_2);
	}
	void InternalBoundaries::add_P_2_to_global(int element_number, int neighbor_element_number, Matrix& A, vector <vector<double>>& P_2)
	{
		int id_i, id_j;
		Element element = elements[element_number];
		Element neighbor_element = elements[neighbor_element_number];
		int n_edges = elements.size() * 4;

		for(int i = 0; i < n_func_p; i++)
		{
			id_i = element.nodes[i] + n_edges;
			for(int j = 0; j < n_func_u; j++)
			{				
				id_j = element.edges[j];
				A.add_element(id_i, id_j, P_2[i][j]); 
			}

			for(int j = n_func_u; j < n_func_u * 2; j++)
			{				
				id_j = neighbor_element.edges[j - n_func_u];
				A.add_element(id_i, id_j, P_2[i][j]); 
			}
		}

		for(int i = n_func_p; i < n_func_p * 2; i++)
		{
			id_i = neighbor_element.nodes[i - n_func_p] + n_edges;
			for(int j = 0; j < n_func_u; j++)
			{				
				id_j = element.edges[j];
				A.add_element(id_i, id_j, P_2[i][j]); 
			}

			for(int j = n_func_u; j < n_func_u * 2; j++)
			{				
				id_j = neighbor_element.edges[j - n_func_u];
				A.add_element(id_i, id_j, P_2[i][j]); 
			}
		}
	}

	void InternalBoundaries::calculate_SP_horizontal(int element_number1, int element_number2, matrix::Matrix& A)
	{
		vector <vector<double>> SNN, SNK, SKN, SKK, SP;

		SNN.resize(n_func_p);
		SNK.resize(n_func_p);
		SKN.resize(n_func_p);
		SKK.resize(n_func_p);
		SP.resize(n_func_p * 2);

		for (int i = 0; i < n_func_p; i++)
		{
			initialize_vector(SNN[i], n_func_p);
			initialize_vector(SNK[i], n_func_p);
			initialize_vector(SKN[i], n_func_p);
			initialize_vector(SKK[i], n_func_p);
			initialize_vector(SP[i], n_func_p * 2);
			initialize_vector(SP[i + n_func_p], n_func_p * 2);
		}

		double hx = get_hx(element_number1);

		double n_vec[2] = {0, 1};

		double jacobian = 0.5 * hx;
		double st = jacobian * mu2;

		for(int i = 0; i < 4; i++)
		{
			for(int j = 0; j < 4; j++)
			{
				for(int k = 0; k < 3; k++)
				{
					double p_ksi = gauss_points_1[k];
					SNN[i][j] += gauss_weights_1[k] * psi[i](p_ksi, 1) * psi[j](p_ksi, 1);
					SKK[i][j] += gauss_weights_1[k] * psi[i](p_ksi, -1) * psi[j](p_ksi, -1);
					SNK[i][j] += gauss_weights_1[k] * psi[i](p_ksi, 1) * psi[j](p_ksi, -1);
					SKN[i][j] += gauss_weights_1[k] * psi[i](p_ksi, -1) * psi[j](p_ksi, 1);

				} 
				SNN[i][j] *= st;
				SNK[i][j] *= -st;
				SKK[i][j] *= st;
				SKN[i][j] *= -st;
			}
		}

		for(int i = 0; i < n_func_p; i++)
			for(int j = 0; j < n_func_p; j++)
			{
				SP[i][j] = SNN[i][j];
				SP[i + n_func_p][j + n_func_p] = SKK[i][j];
				SP[i][j + n_func_p] = SNK[i][j];
				SP[i + n_func_p][j] = SKN[i][j];
			}

		add_SP_to_global(element_number1, element_number2, A, SP);
	}
	void InternalBoundaries::calculate_SP_vertical(int element_number1, int element_number2, matrix::Matrix& A)
	{
		vector <vector<double>> SNN, SNK, SKN, SKK, SP;

		SNN.resize(n_func_p);
		SNK.resize(n_func_p);
		SKN.resize(n_func_p);
		SKK.resize(n_func_p);
		SP.resize(n_func_p * 2);

		for (int i = 0; i < n_func_p; i++)
		{
			initialize_vector(SNN[i], n_func_p);
			initialize_vector(SNK[i], n_func_p);
			initialize_vector(SKN[i], n_func_p);
			initialize_vector(SKK[i], n_func_p);
			initialize_vector(SP[i], n_func_p * 2);
			initialize_vector(SP[i + n_func_p], n_func_p * 2);
		}

		double hy = get_hy(element_number1);

		double n_vec[2] = {1, 0};

		double jacobian = 0.5 * hy;
		double st = jacobian * mu2;

		for(int i = 0; i < n_func_p; i++)
		{
			for(int j = 0; j < n_func_p; j++)
			{
				for(int k = 0; k < 3; k++)
				{
					double p_etta = gauss_points_1[k];
					SNN[i][j] += gauss_weights_1[k] * psi[i](1, p_etta) * psi[j](1, p_etta);
					SKK[i][j] += gauss_weights_1[k] * psi[i](-1, p_etta) * psi[j](-1, p_etta);
					SNK[i][j] += gauss_weights_1[k] * psi[i](1, p_etta) * psi[j](-1, p_etta);
					SKN[i][j] += gauss_weights_1[k] * psi[i](-1, p_etta) * psi[j](1, p_etta);

				} 
				SNN[i][j] *= st;
				SNK[i][j] *= -st;
				SKK[i][j] *= st;
				SKN[i][j] *= -st;
			}
		}

		for(int i = 0; i < n_func_p; i++)
			for(int j = 0; j < n_func_p; j++)
			{
				SP[i][j] = SNN[i][j];
				SP[i + n_func_p][j + n_func_p] = SKK[i][j];
				SP[i][j + n_func_p] = SNK[i][j];
				SP[i + n_func_p][j] = SKN[i][j];
			}

		add_SP_to_global(element_number1, element_number2, A, SP);
	}
	void InternalBoundaries::add_SP_to_global(int element_number, int neighbor_element_number, Matrix& A, vector <vector<double>>& SP)
	{
		int id_i, id_j;
		Element element = elements[element_number];
		Element neighbor_element = elements[neighbor_element_number];
		int n_edges = elements.size() * 4;

		for(int i = 0; i < n_func_p; i++)
		{
			id_i = element.nodes[i] + n_edges;
			for(int j = 0; j < n_func_p; j++)
			{				
				id_j = element.nodes[j] + n_edges;
				A.add_element(id_i, id_j, SP[i][j]); 
			}

			for(int j = n_func_p; j < n_func_p * 2; j++)
			{				
				id_j = neighbor_element.nodes[j - n_func_p] + n_edges;
				A.add_element(id_i, id_j, SP[i][j]); 
			}
		}

		for(int i = n_func_p; i < n_func_p * 2; i++)
		{
			id_i = neighbor_element.nodes[i - n_func_p] + n_edges;
			for(int j = 0; j < n_func_p; j++)
			{				
				id_j = element.nodes[j] + n_edges;
				A.add_element(id_i, id_j, SP[i][j]); 
			}

			for(int j = n_func_p; j < n_func_p * 2; j++)
			{				
				id_j = neighbor_element.nodes[j - n_func_p] + n_edges;
				A.add_element(id_i, id_j, SP[i][j]); 
			}
		}
	}
}



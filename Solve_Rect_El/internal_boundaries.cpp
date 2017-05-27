#include "boundaries.h"
#include "element.h"
#include "myfunctions.h"
#include "point.h"

#include "rapidjson/document.h"
#include <fstream>

#include <array>

using namespace std;
using namespace element;
using namespace partition;
using namespace parameters;
using namespace basis;
using namespace matrix;
using namespace point;
using namespace integration;

using namespace rapidjson;

namespace boundaries
{
	void InternalBoundaries::initialize_penalty_parameters(std::string fileName)
	{
		ifstream fileIn(fileName);
		string json;
		string parameterName;

		while (!fileIn.eof())
		{
			string line;
			fileIn >> line;
			json += line;
		}

		fileIn.close();

		Document docIn;
		docIn.Parse(json.c_str());

		auto& jPenaltyParameters = docIn["penaltyParameters"];
		gamma = jPenaltyParameters["gamma"].GetDouble();
		sigma = jPenaltyParameters["sigma"].GetDouble();
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
				EdgeOrient orient;

				//левый/правый сосед->вертикальная граница
				if(k == 0 || k == 1)
				{
					orient = VERTICAL;
				}
				else
				{
					orient = HORIZONTAL;
				}

				calculate_ES(element_number, neighbor_element, A, orient);
				calculate_P_1(element_number, neighbor_element, A, orient);
				calculate_P_2(element_number, neighbor_element, A, orient);
				calculate_SP(element_number, neighbor_element, A, orient);
			}
		}
	}

	void InternalBoundaries::calculate_ES(int element_number1, int element_number2, matrix::Matrix& A, EdgeOrient orient)
	{		
		double EKK[n_func_u][n_func_u];
		double ENN[n_func_u][n_func_u];
		double EKN[n_func_u][n_func_u];
		double ENK[n_func_u][n_func_u];
		double SNN[n_func_u][n_func_u];
		double SNK[n_func_u][n_func_u];
		double SKN[n_func_u][n_func_u];
		double SKK[n_func_u][n_func_u];
		vector <vector<double>> ES(n_func_u * 2);

		for (int i = 0; i < n_func_u; i++)
		{
			memset(&EKK[i][0], 0, sizeof(double) * n_func_u);
			memset(&ENN[i][0], 0, sizeof(double) * n_func_u);
			memset(&EKN[i][0], 0, sizeof(double) * n_func_u);
			memset(&ENK[i][0], 0, sizeof(double) * n_func_u);
			memset(&SNN[i][0], 0, sizeof(double) * n_func_u);
			memset(&SNK[i][0], 0, sizeof(double) * n_func_u);
			memset(&SKN[i][0], 0, sizeof(double) * n_func_u);
			memset(&SKK[i][0], 0, sizeof(double) * n_func_u);
			initialize_vector(ES[i], n_func_u * 2);
			initialize_vector(ES[i + n_func_u], n_func_u * 2);
		}
			
		Element element = elements[element_number1];
		Element element_2 = elements[element_number2];
		double lambdaN = calculate_lambda(element.number_of_area);
		double lambdaK = calculate_lambda(element_2.number_of_area);
		double lambda_s = min(lambdaN, lambdaK);
		double hxN = get_hx(element_number1);
		double hyN = get_hy(element_number1);
		double hxK = get_hx(element_number2);
		double hyK = get_hy(element_number2);

		double a, jacobian, h, p_ksi, p_etta, signKsi, signEtta;
		Point nN;

		if(orient == HORIZONTAL)
		{ 
			a = 0.25 * hxN; //якобиан*0.5
			jacobian = 0.5 * hxN;
			h = min(hyN, hyK);
			nN.x = 0; nN.y = 1;
			signKsi = 1;
			signEtta = -1;
			p_etta = 1;
		}
		else
		{
			a = 0.25 * hyN; //якобиан*0.5
			jacobian = 0.5 * hyN;
			h = min(hxN, hxK);
			nN.x = 1; nN.y = 0;
			signKsi = -1;
			signEtta = 1;
			p_ksi = 1;
		}

		int kK = element.order, kN = element_2.order;
		int k = max(kK, kN);
		double c = gamma * k * k / h;
		double st = jacobian * c * lambda_s;

		for(int i = 0; i < n_func_u; i++)
		{
			for(int j = 0; j < n_func_u; j++)
			{
				for(int k = 0; k < n_ip1D; k++)
				{
					if (orient == HORIZONTAL)
						p_ksi = gauss_points_1[k];
					else
						p_etta = gauss_points_1[k];

					ENN[i][j] += gauss_weights_1[k] * 
							   (phix[j](p_ksi, p_etta) * nN.x * dphixksi[i](p_ksi, p_etta) / hxN +
								phix[j](p_ksi, p_etta) * nN.y * dphixetta[i](p_ksi, p_etta) / hyN +
								phiy[j](p_ksi, p_etta) * nN.x * dphiyksi[i](p_ksi, p_etta) / hxN +
								phiy[j](p_ksi, p_etta) * nN.y * dphiyetta[i](p_ksi, p_etta) / hyN +
								phix[i](p_ksi, p_etta) * nN.x * dphixksi[j](p_ksi, p_etta) / hxN +
								phix[i](p_ksi, p_etta) * nN.y * dphixetta[j](p_ksi, p_etta) / hyN +
								phiy[i](p_ksi, p_etta) * nN.x * dphiyksi[j](p_ksi, p_etta) / hxN +
								phiy[i](p_ksi, p_etta) * nN.y * dphiyetta[j](p_ksi, p_etta) / hyN);
					EKK[i][j] += gauss_weights_1[k] * 
							   (phix[j](p_ksi * signKsi, p_etta * signEtta) * nN.x * dphixksi[i](p_ksi * signKsi, p_etta * signEtta) / hxK +
								phix[j](p_ksi * signKsi, p_etta * signEtta) * nN.y * dphixetta[i](p_ksi * signKsi, p_etta * signEtta) / hyK +
								phiy[j](p_ksi * signKsi, p_etta * signEtta) * nN.x * dphiyksi[i](p_ksi * signKsi, p_etta * signEtta) / hxK +
								phiy[j](p_ksi * signKsi, p_etta * signEtta) * nN.y * dphiyetta[i](p_ksi * signKsi, p_etta * signEtta) / hyK +
								phix[i](p_ksi * signKsi, p_etta * signEtta) * nN.x * dphixksi[j](p_ksi * signKsi, p_etta * signEtta) / hxK +
								phix[i](p_ksi * signKsi, p_etta * signEtta) * nN.y * dphixetta[j](p_ksi * signKsi, p_etta * signEtta) / hyK +
								phiy[i](p_ksi * signKsi, p_etta * signEtta) * nN.x * dphiyksi[j](p_ksi * signKsi, p_etta * signEtta) / hxK +
								phiy[i](p_ksi * signKsi, p_etta * signEtta) * nN.y * dphiyetta[j](p_ksi * signKsi, p_etta * signEtta) / hyK);
					ENK[i][j] += gauss_weights_1[k] *
							 (-(phix[j](p_ksi * signKsi, p_etta * signEtta) * nN.x * dphixksi[i](p_ksi, p_etta) / hxN -
								phix[j](p_ksi * signKsi, p_etta * signEtta) * nN.y * dphixetta[i](p_ksi, p_etta) / hyN -
								phiy[j](p_ksi * signKsi, p_etta * signEtta) * nN.x * dphiyksi[i](p_ksi, p_etta) / hxN -
								phiy[j](p_ksi * signKsi, p_etta * signEtta) * nN.y * dphiyetta[i](p_ksi, p_etta) / hyN) * lambdaN +
							   (phix[i](p_ksi, p_etta) * nN.x * dphixksi[j](p_ksi * signKsi, p_etta * signEtta) / hxK +
								phix[i](p_ksi, p_etta) * nN.y * dphixetta[j](p_ksi * signKsi, p_etta * signEtta) / hyK +
								phiy[i](p_ksi, p_etta) * nN.x * dphiyksi[j](p_ksi * signKsi, p_etta * signEtta) / hxK +
								phiy[i](p_ksi, p_etta) * nN.y * dphiyetta[j](p_ksi * signKsi, p_etta * signEtta) / hyK) * lambdaK);
					EKN[i][j] += gauss_weights_1[k] *
							  ((phix[j](p_ksi, p_etta) * nN.x * dphixksi[i](p_ksi * signKsi, p_etta * signEtta) / hxK +
								phix[j](p_ksi, p_etta) * nN.y * dphixetta[i](p_ksi * signKsi, p_etta * signEtta) / hyK +
								phiy[j](p_ksi, p_etta) * nN.x * dphiyksi[i](p_ksi * signKsi, p_etta * signEtta) / hxK +
								phiy[j](p_ksi, p_etta) * nN.y * dphiyetta[i](p_ksi * signKsi, p_etta * signEtta) / hyK) * lambdaK -
							   (phix[i](p_ksi * signKsi, p_etta * signEtta) * nN.x * dphixksi[j](p_ksi, p_etta) / hxN -
								phix[i](p_ksi * signKsi, p_etta * signEtta) * nN.y * dphixetta[j](p_ksi, p_etta) / hyN -
								phiy[i](p_ksi * signKsi, p_etta * signEtta) * nN.x * dphiyksi[j](p_ksi, p_etta) / hxN -
								phiy[i](p_ksi * signKsi, p_etta * signEtta) * nN.y * dphiyetta[j](p_ksi, p_etta) / hyN) * lambdaN);
					SNN[i][j] += gauss_weights_1[k] * 
								(phix[j](p_ksi, p_etta) * phix[i](p_ksi, p_etta) * nN.x +
								 phix[j](p_ksi, p_etta) * phix[i](p_ksi, p_etta) * nN.y +
								 phiy[j](p_ksi, p_etta) * phiy[i](p_ksi, p_etta) * nN.x +
								 phiy[j](p_ksi, p_etta) * phiy[i](p_ksi, p_etta) * nN.y);
					SNK[i][j] += gauss_weights_1[k] * 
								(phix[j](p_ksi * signKsi, p_etta * signEtta) * phix[i](p_ksi, p_etta) * nN.x +
								 phix[j](p_ksi * signKsi, p_etta * signEtta) * phix[i](p_ksi, p_etta) * nN.y +
								 phiy[j](p_ksi * signKsi, p_etta * signEtta) * phiy[i](p_ksi, p_etta) * nN.x +
								 phiy[j](p_ksi * signKsi, p_etta * signEtta) * phiy[i](p_ksi, p_etta) * nN.y);
					SKK[i][j] += gauss_weights_1[k] * 
								(phix[j](p_ksi * signKsi, p_etta * signEtta) * phix[i](p_ksi * signKsi, p_etta * signEtta) * nN.x +
								 phix[j](p_ksi * signKsi, p_etta * signEtta) * phix[i](p_ksi * signKsi, p_etta * signEtta) * nN.y +
								 phiy[j](p_ksi * signKsi, p_etta * signEtta) * phiy[i](p_ksi * signKsi, p_etta * signEtta) * nN.x +
								 phiy[j](p_ksi * signKsi, p_etta * signEtta) * phiy[i](p_ksi * signKsi, p_etta * signEtta) * nN.y);
					SKN[i][j] += gauss_weights_1[k] * 
								(phix[j](p_ksi, p_etta) * phix[i](p_ksi * signKsi, p_etta * signEtta) * nN.x +
								 phix[j](p_ksi, p_etta) * phix[i](p_ksi * signKsi, p_etta * signEtta) * nN.y +
								 phiy[j](p_ksi, p_etta) * phiy[i](p_ksi * signKsi, p_etta * signEtta) * nN.x +
								 phiy[j](p_ksi, p_etta) * phiy[i](p_ksi * signKsi, p_etta * signEtta) * nN.y);

				} 
				ENN[i][j] *= a * lambdaN;
				EKK[i][j] *= -a * lambdaK;
				ENK[i][j] *= a;
				EKN[i][j] *= a;
				SNN[i][j] *= st;
				SNK[i][j] *= -st;
				SKK[i][j] *= st;
				SKN[i][j] *= -st;
			}
		}

		for(int i = 0; i < n_func_u; i++)
			for(int j = 0; j < n_func_u; j++)
			{
				ES[i][j] = -ENN[i][j] + SNN[i][j];
				ES[i + n_func_u][j + n_func_u] = -EKK[i][j] + SKK[i][j];
				ES[i][j + n_func_u] = -ENK[i][j] + SNK[i][j];
				ES[i + n_func_u][j] = -EKN[i][j] + SKN[i][j];
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
			id_i = element.dof_u[i];
			for(int j = 0; j < n_func_u; j++)
			{				
				id_j = element.dof_u[j];
				A.add_element(id_i, id_j, ES[i][j]); 
			}

			for(int j = n_func_u; j < n_func_u * 2; j++)
			{				
				id_j = neighbor_element.dof_u[j - n_func_u];
				A.add_element(id_i, id_j, ES[i][j]); 
			}
		}

		for(int i = n_func_u; i < n_func_u * 2; i++)
		{
			id_i = neighbor_element.dof_u[i - n_func_u];
			for(int j = 0; j < n_func_u; j++)
			{				
				id_j = element.dof_u[j];
				A.add_element(id_i, id_j, ES[i][j]); 
			}

			for(int j = n_func_u; j < n_func_u * 2; j++)
			{				
				id_j = neighbor_element.dof_u[j - n_func_u];
				A.add_element(id_i, id_j, ES[i][j]); 
			}
		}
	}

	void InternalBoundaries::calculate_P_1(int element_number1, int element_number2, matrix::Matrix& A, EdgeOrient orient)
	{
		double PKK[n_func_u][n_func_p];
		double PNN[n_func_u][n_func_p];
		double PKN[n_func_u][n_func_p];
		double PNK[n_func_u][n_func_p];
		vector <vector<double>> P_1(n_func_u * 2);

		for (int i = 0; i < n_func_u; i++)
		{
			memset(&PKK[i][0], 0, sizeof(double) * n_func_p);
			memset(&PNN[i][0], 0, sizeof(double) * n_func_p);
			memset(&PKN[i][0], 0, sizeof(double) * n_func_p);
			memset(&PNK[i][0], 0, sizeof(double) * n_func_p);
			initialize_vector(P_1[i], n_func_p * 2);
			initialize_vector(P_1[i + n_func_u], n_func_p * 2);
		}


		Element element = elements[element_number1];
		Element element_2 = elements[element_number2];
		double rho = calculate_rho(element.number_of_area);
		double rho_2 = calculate_rho(element_2.number_of_area);
		
		double h, p_ksi, p_etta, signKsi, signEtta;
		Point nN;

		if (orient == HORIZONTAL)
		{
			h = get_hx(element_number1);
			nN.x = 0; nN.y = 1;
			signKsi = 1;
			signEtta = -1;
			p_etta = 1;
		}
		else
		{
			h = get_hy(element_number1);
			nN.x = 1; nN.y = 0;
			signKsi = -1;
			signEtta = 1;
			p_ksi = 1;
		}

		double a1 = 0.25 * h * (1 / rho); //якобиан*0.5*(1/rho)
		double a2 = 0.25 * h * (1 / rho_2);

		for(int i = 0; i < n_func_u; i++)
		{
			for(int j = 0; j < n_func_p; j++)
			{
				for(int k = 0; k < n_ip1D; k++)
				{
					if (orient == HORIZONTAL)
						p_ksi = gauss_points_1[k];
					else
						p_etta = gauss_points_1[k];

					PNN[i][j] += gauss_weights_1[k] * psi[j](p_ksi, p_etta) *
							   (phix[i](p_ksi, p_etta) * nN.x + phiy[i](p_ksi, p_etta) * nN.y);
					PKK[i][j] += gauss_weights_1[k] * psi[j](p_ksi * signKsi, p_etta * signEtta) *
							   (phix[i](p_ksi * signKsi, p_etta * signEtta) * nN.x + phiy[i](p_ksi * signKsi, p_etta * signEtta) * nN.y);
					PNK[i][j] += gauss_weights_1[k] * psi[j](p_ksi * signKsi, p_etta * signEtta) *
							   (phix[i](p_ksi, p_etta) * nN.x + phiy[i](p_ksi, p_etta) * nN.y);
					PKN[i][j] += gauss_weights_1[k] * psi[j](p_ksi, p_etta) *
							   (phix[i](p_ksi * signKsi, p_etta * signEtta) * nN.x + phiy[i](p_ksi * signKsi, p_etta * signEtta) * nN.y);

				} 
				PNN[i][j] *= a1;
				PKK[i][j] *= -a2;
				PNK[i][j] *= a2;
				PKN[i][j] *= -a1;
			}
		}

		for(int i = 0; i < n_func_u; i++)
			for(int j = 0; j < n_func_p; j++)
			{
				P_1[i][j] = PNN[i][j];
				P_1[i + n_func_u][j + n_func_p] = PKK[i][j];
				P_1[i][j + n_func_p] = PNK[i][j];
				P_1[i + n_func_u][j] = PKN[i][j];
			}

		add_P_1_to_global(element_number1, element_number2, A, P_1);
	}
	void InternalBoundaries::add_P_1_to_global(int element_number, int neighbor_element_number, Matrix& A, vector <vector<double>>& P_1)
	{
		int id_i, id_j;
		Element element = elements[element_number];
		Element neighbor_element = elements[neighbor_element_number];

		for(int i = 0; i < n_func_u; i++)
		{
			id_i = element.dof_u[i];
			for(int j = 0; j < n_func_p; j++)
			{				
				id_j = element.dof_p[j];
				A.add_element(id_i, id_j, P_1[i][j]); 
			}

			for(int j = n_func_p; j < n_func_p * 2; j++)
			{				
				id_j = neighbor_element.dof_p[j - n_func_p];
				A.add_element(id_i, id_j, P_1[i][j]); 
			}
		}

		for(int i = n_func_u; i <  n_func_u * 2; i++)
		{
			id_i = neighbor_element.dof_u[i - n_func_u];
			for(int j = 0; j < n_func_p; j++)
			{				
				id_j = element.dof_p[j];
				A.add_element(id_i, id_j, P_1[i][j]); 
			}

			for(int j = n_func_p; j < n_func_p * 2; j++)
			{				
				id_j = neighbor_element.dof_p[j - n_func_p];
				A.add_element(id_i, id_j, P_1[i][j]); 
			}
		}
	}

	void InternalBoundaries::calculate_P_2(int element_number1, int element_number2, matrix::Matrix& A, EdgeOrient orient)
	{
		double PKK[n_func_p][n_func_u];
		double PNN[n_func_p][n_func_u];
		double PKN[n_func_p][n_func_u];
		double PNK[n_func_p][n_func_u];
		vector <vector<double>> P_2(n_func_p * 2);

		for (int i = 0; i < n_func_p; i++)
		{
			memset(&PKK[i][0], 0, sizeof(double) * n_func_u);
			memset(&PNN[i][0], 0, sizeof(double) * n_func_u);
			memset(&PKN[i][0], 0, sizeof(double) * n_func_u);
			memset(&PNK[i][0], 0, sizeof(double) * n_func_u);
			initialize_vector(P_2[i], n_func_u * 2);
			initialize_vector(P_2[i + n_func_p], n_func_u * 2);
		}


		Element element = elements[element_number1];
		Element element_2 = elements[element_number2];

		double h, p_ksi, p_etta, signKsi, signEtta;
		Point nN;

		if (orient == HORIZONTAL)
		{
			h = get_hx(element_number1);
			nN.x = 0; nN.y = 1;
			signKsi = 1;
			signEtta = -1;
			p_etta = 1;
		}
		else
		{
			h = get_hy(element_number1);
			nN.x = 1; nN.y = 0;
			signKsi = -1;
			signEtta = 1;
			p_ksi = 1;
		}

		double a = 0.25 * h; //якобиан*0.5

		for(int i = 0; i < n_func_p; i++)
		{
			for(int j = 0; j < n_func_u; j++)
			{
				for(int k = 0; k < n_ip1D; k++)
				{
					if (orient == HORIZONTAL)
						p_ksi = gauss_points_1[k];
					else
						p_etta = gauss_points_1[k];
					PNN[i][j] += gauss_weights_1[k] * psi[i](p_ksi, p_etta) *
							   (phix[j](p_ksi, p_etta) * nN.x + phiy[j](p_ksi, p_etta) * nN.y);
					PKK[i][j] += gauss_weights_1[k] * psi[i](p_ksi * signKsi, p_etta * signEtta) *
							   (phix[j](p_ksi * signKsi, p_etta * signEtta) * nN.x + phiy[j](p_ksi * signKsi, p_etta * signEtta) * nN.y);
					PNK[i][j] += gauss_weights_1[k] * psi[i](p_ksi, p_etta) *
							   (phix[j](p_ksi * signKsi, p_etta * signEtta) * nN.x + phiy[j](p_ksi * signKsi, p_etta * signEtta) * nN.y);
					PKN[i][j] += gauss_weights_1[k] * psi[i](p_ksi * signKsi, p_etta * signEtta) *
							   (phix[j](p_ksi, p_etta) * nN.x + phiy[j](p_ksi, p_etta) * nN.y);

				} 
				PNN[i][j] *= a;
				PKK[i][j] *= -a;
				PNK[i][j] *= -a;
				PKN[i][j] *= a;
			}
		}

		for(int i = 0; i < n_func_p; i++)
			for(int j = 0; j < n_func_u; j++)
			{
				P_2[i][j] = -PNN[i][j];
				P_2[i + n_func_p][j + n_func_u] = -PKK[i][j];
				P_2[i][j + n_func_u] = -PNK[i][j];
				P_2[i + n_func_p][j] = -PKN[i][j];
			}

		add_P_2_to_global(element_number1, element_number2, A, P_2);
	}
	void InternalBoundaries::add_P_2_to_global(int element_number, int neighbor_element_number, Matrix& A, vector <vector<double>>& P_2)
	{
		int id_i, id_j;
		Element element = elements[element_number];
		Element neighbor_element = elements[neighbor_element_number];

		for(int i = 0; i < n_func_p; i++)
		{
			id_i = element.dof_p[i];
			for(int j = 0; j < n_func_u; j++)
			{				
				id_j = element.dof_u[j];
				A.add_element(id_i, id_j, P_2[i][j]); 
			}

			for(int j = n_func_u; j < n_func_u * 2; j++)
			{				
				id_j = neighbor_element.dof_u[j - n_func_u];
				A.add_element(id_i, id_j, P_2[i][j]); 
			}
		}

		for(int i = n_func_p; i < n_func_p * 2; i++)
		{
			id_i = neighbor_element.dof_p[i - n_func_p];
			for(int j = 0; j < n_func_u; j++)
			{				
				id_j = element.dof_u[j];
				A.add_element(id_i, id_j, P_2[i][j]); 
			}

			for(int j = n_func_u; j < n_func_u * 2; j++)
			{				
				id_j = neighbor_element.dof_u[j - n_func_u];
				A.add_element(id_i, id_j, P_2[i][j]); 
			}
		}
	}

	void InternalBoundaries::calculate_SP(int element_number1, int element_number2, matrix::Matrix& A, EdgeOrient orient)
	{
		double SNN[n_func_p][n_func_p], SNK[n_func_p][n_func_p], SKN[n_func_p][n_func_p], SKK[n_func_p][n_func_p];
		vector <vector<double>> SP(n_func_p * 2);

		for (int i = 0; i < n_func_p; i++)
		{
			memset(&SNN[i][0], 0, sizeof(double) * n_func_p);
			memset(&SNK[i][0], 0, sizeof(double) * n_func_p);
			memset(&SKN[i][0], 0, sizeof(double) * n_func_p);
			memset(&SKK[i][0], 0, sizeof(double) * n_func_p);
			initialize_vector(SP[i], n_func_p * 2);
			initialize_vector(SP[i + n_func_p], n_func_p * 2);
		}

		Element element1 = elements[element_number1];
		Element element2 = elements[element_number2];
		double lambdaN = calculate_lambda(element1.number_of_area);
		double lambdaK = calculate_lambda(element2.number_of_area);
		double lambda = min(lambdaN, lambdaK);

		double h, p_ksi, p_etta, signKsi, signEtta;
		Point nN;

		if (orient == HORIZONTAL)
		{
			h = get_hx(element_number1);
			nN.x = 0; nN.y = 1;
			signKsi = 1;
			signEtta = -1;
			p_etta = 1;
		}
		else
		{
			h = get_hy(element_number1);
			nN.x = 1; nN.y = 0;
			signKsi = -1;
			signEtta = 1;
			p_ksi = 1;
		}

		double jacobian = 0.5 * h; //якобиан
		double st = jacobian * 1.0 / lambda * sigma;

		for(int i = 0; i < n_func_p; i++)
		{
			for(int j = 0; j < n_func_p; j++)
			{
				for(int k = 0; k < n_ip1D; k++)
				{
					if (orient == HORIZONTAL)
						p_ksi = gauss_points_1[k];
					else
						p_etta = gauss_points_1[k];
					SNN[i][j] += gauss_weights_1[k] * psi[i](p_ksi, p_etta) * psi[j](p_ksi, p_etta);
					SKK[i][j] += gauss_weights_1[k] * psi[i](p_ksi * signKsi, p_etta * signEtta) * psi[j](p_ksi * signKsi, p_etta * signEtta);
					SNK[i][j] += gauss_weights_1[k] * psi[i](p_ksi, p_etta) * psi[j](p_ksi * signKsi, p_etta * signEtta);
					SKN[i][j] += gauss_weights_1[k] * psi[i](p_ksi * signKsi, p_etta * signEtta) * psi[j](p_ksi, p_etta);

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

		for(int i = 0; i < n_func_p; i++)
		{
			id_i = element.dof_p[i];
			for(int j = 0; j < n_func_p; j++)
			{				
				id_j = element.dof_p[j];
				A.add_element(id_i, id_j, SP[i][j]); 
			}

			for(int j = n_func_p; j < n_func_p * 2; j++)
			{				
				id_j = neighbor_element.dof_p[j - n_func_p];
				A.add_element(id_i, id_j, SP[i][j]); 
			}
		}

		for(int i = n_func_p; i < n_func_p * 2; i++)
		{
			id_i = neighbor_element.dof_p[i - n_func_p];
			for(int j = 0; j < n_func_p; j++)
			{				
				id_j = element.dof_p[j];
				A.add_element(id_i, id_j, SP[i][j]); 
			}

			for(int j = n_func_p; j < n_func_p * 2; j++)
			{				
				id_j = neighbor_element.dof_p[j - n_func_p];
				A.add_element(id_i, id_j, SP[i][j]); 
			}
		}
	}
}



#include "boundaries.h"
#include "element.h"
#include "myfunctions.h"
#include "point.h"

#include "rapidjson/document.h"
#include <fstream>

using namespace std;
using namespace element;
using namespace partition;
using namespace parameters;
using namespace basis;
using namespace matrix;
using namespace point;

using namespace rapidjson;

namespace boundaries
{

	void OuterBoundaries::initialize_penalty_parameters(std::string fileName)
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
	}

	void OuterBoundaries::calculate_outer_boundaries(int element_number, matrix::Matrix& A)
	{
		Element element = elements[element_number];

		int last_node = nodes.size() - 1;
		int left_low_corner_node = element.nodes[0];
		int right_up_corner_node = element.nodes[3];

		EdgeSide edgeSide;

		if(nodes[left_low_corner_node].x == nodes[0].x)
		{
			edgeSide = LEFT;
			calculate_ES_out(element_number, A, edgeSide);
			calculate_P_1_out(element_number, A, edgeSide);
			calculate_P_2_out(element_number, A, edgeSide);
			calculate_SP_out(element_number, A, edgeSide);
		}//вертикальная левая граница
		if(nodes[right_up_corner_node].x == nodes[last_node].x)
		{
			edgeSide = RIGHT;
			calculate_ES_out(element_number, A, edgeSide);
			calculate_P_1_out(element_number, A, edgeSide);
			calculate_P_2_out(element_number, A, edgeSide);
			calculate_SP_out(element_number, A, edgeSide);
		}//вертикальная правая граница
		if(nodes[left_low_corner_node].y == nodes[0].y) 
		{
			edgeSide = LOW;
			calculate_ES_out(element_number, A, edgeSide);
			calculate_P_1_out(element_number, A, edgeSide);
			calculate_P_2_out(element_number, A, edgeSide);
			calculate_SP_out(element_number, A, edgeSide);
		}//горизонтальная нижняя граница
		if(nodes[right_up_corner_node].y == nodes[last_node].y)	
		{
			edgeSide = UP;
			calculate_ES_out(element_number, A, edgeSide);
			calculate_P_1_out(element_number, A, edgeSide);
			calculate_P_2_out(element_number, A, edgeSide);
			calculate_SP_out(element_number, A, edgeSide);
		}//горизонтальная верхняя граница
	}

	void OuterBoundaries::calculate_ES_out(int element_number, Matrix& A, EdgeSide edgeSide)
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

		double h, p_ksi, p_etta;
		Point nEl, nEl2;

		if (edgeSide == LEFT || edgeSide == RIGHT)
		{			
			h = hy;
			if (edgeSide == LEFT)
			{
				nEl.x = -1; 
				p_ksi = -1;
			}
			else
			{
				nEl.x = 1;
				p_ksi = 1;
			}
			nEl.y = 0;
			nEl2.x = 1;
			nEl2.y = 0;
		}
		else
		{			
			h = hx;
			if (edgeSide == LOW)
			{
				nEl.y = -1;
				p_etta = -1;
			}
			else
			{
				nEl.y = 1;
				p_etta = 1;
			}
			nEl.x = 0;
			nEl2.x = 0;
			nEl2.y = 1;
		}

		double jacobian = 0.5 * h;
		double a = jacobian * lambda;
		int k = 1;
		double c = gamma * k * k / h;
		double st = jacobian * lambda * c;


		for(int i = 0; i < n_func_u; i++)
		{
			for(int j = 0; j < n_func_u; j++)
			{
				for(int k = 0; k < 3; k++)
				{
					if (edgeSide == LOW || edgeSide == UP)
						p_ksi = gauss_points_1[k];
					else
						p_etta = gauss_points_1[k];

					E_out[i][j] += gauss_weights_1[k] * 
								  (phix[j](p_ksi, p_etta) * nEl.x * dphixksi[i](p_ksi, p_etta) / hx +
								   phix[j](p_ksi, p_etta) * nEl.y * dphixetta[i](p_ksi, p_etta) / hy +
								   phiy[j](p_ksi, p_etta) * nEl.x * dphiyksi[i](p_ksi, p_etta) / hx +
								   phiy[j](p_ksi, p_etta) * nEl.y * dphiyetta[i](p_ksi, p_etta) / hy -
								   phix[i](p_ksi, p_etta) * nEl.x * dphixksi[j](p_ksi, p_etta) / hx -
								   phix[i](p_ksi, p_etta) * nEl.y * dphixetta[j](p_ksi, p_etta) / hy -
								   phiy[i](p_ksi, p_etta) * nEl.x * dphiyksi[j](p_ksi, p_etta) / hx -
								   phiy[i](p_ksi, p_etta) * nEl.y * dphiyetta[j](p_ksi, p_etta) / hy);

					S_out[i][j] += gauss_weights_1[k] * 
								  (phix[j](p_ksi, p_etta) * phix[i](p_ksi, p_etta) * nEl2.x +
								   phix[j](p_ksi, p_etta) * phix[i](p_ksi, p_etta) * nEl2.y +
								   phiy[j](p_ksi, p_etta) * phiy[i](p_ksi, p_etta) * nEl2.x +
								   phiy[j](p_ksi, p_etta) * phiy[i](p_ksi, p_etta) * nEl2.y);
				} 
				E_out[i][j] *= a;
				S_out[i][j] *= st;
			}
		}

		for(int i = 0; i < n_func_u; i++)
		{
			int id_i = element.dof_u[i];
			for(int j = 0; j < n_func_u; j++)
			{				
				int id_j = element.dof_u[j];
				A.add_element(id_i, id_j, E_out[i][j] + S_out[i][j]); 
			}
		}
	}
	void OuterBoundaries::calculate_P_1_out(int element_number, Matrix& A, EdgeSide edgeSide)
	{
		vector <vector<double>> P_1_out;

		P_1_out.resize(n_func_u);

		for (int i = 0; i < n_func_u; i++)
			initialize_vector(P_1_out[i], n_func_p);

		Element element = elements[element_number];
		double rho = calculate_rho(element.number_of_area);

		double a, jacobian, h, p_ksi, p_etta;
		Point nEl;

		if (edgeSide == LEFT || edgeSide == RIGHT)
		{
			h = get_hy(element_number);

			if (edgeSide == LEFT)
			{
				nEl.x = -1;
				p_ksi = -1;
			}
			else
			{
				nEl.x = 1;
				p_ksi = 1;
			}
			nEl.y = 0;
		}
		else
		{
			h = get_hx(element_number);

			if (edgeSide == LOW)
			{
				nEl.y = -1;
				p_etta = -1;
			}
			else
			{
				nEl.y = 1;
				p_etta = 1;
			}
			nEl.x = 0;
		}

		jacobian = 0.5 * h;
		a = jacobian * (1 / rho);

		for(int i = 0; i < n_func_u; i++)
		{
			for(int j = 0; j < n_func_p; j++)
			{
				for(int k = 0; k < 3; k++)
				{
					if (edgeSide == LOW || edgeSide == UP)
						p_ksi = gauss_points_1[k];
					else
						p_etta = gauss_points_1[k];

					P_1_out[i][j] += gauss_weights_1[k]  * psi[j](p_ksi, p_etta) *
									(phix[i](p_ksi, p_etta) * nEl.x + phiy[i](p_ksi, p_etta) * nEl.y);

				} 
				P_1_out[i][j] *= a;
			}
		}

		for(int i = 0; i < n_func_u; i++)
		{
			int id_i = element.dof_u[i];
			for(int j = 0; j < n_func_p; j++)
			{				
				int id_j = element.dof_p[j];
				A.add_element(id_i, id_j, P_1_out[i][j]); 
			}
		}
	}
	void OuterBoundaries::calculate_P_2_out(int element_number, Matrix& A, EdgeSide edgeSide)
	{
		vector <vector<double>> P_2_out;

		P_2_out.resize(n_func_p);

		for (int i = 0; i < n_func_p; i++)
			initialize_vector(P_2_out[i], n_func_u);

		Element element = elements[element_number];

		double a, jacobian, h, p_ksi, p_etta;
		Point nEl;

		if (edgeSide == LEFT || edgeSide == RIGHT)
		{
			h = get_hy(element_number);

			if (edgeSide == LEFT)
			{
				nEl.x = -1;
				p_ksi = -1;
			}
			else
			{
				nEl.x = 1;
				p_ksi = 1;
			}
			nEl.y = 0;
		}
		else
		{
			h = get_hx(element_number);

			if (edgeSide == LOW)
			{
				nEl.y = -1;
				p_etta = -1;
			}
			else
			{
				nEl.y = 1;
				p_etta = 1;
			}
			nEl.x = 0;
		}

		jacobian = 0.5 * h;
		a = -jacobian;


		for(int i = 0; i < n_func_p; i++)
		{
			for(int j = 0; j < n_func_u; j++)
			{
				for(int k = 0; k < 3; k++)
				{
					if (edgeSide == LOW || edgeSide == UP)
						p_ksi = gauss_points_1[k];
					else
						p_etta = gauss_points_1[k];

					P_2_out[i][j] += gauss_weights_1[k]  * psi[i](p_ksi, p_etta) *
									(phix[j](p_ksi, p_etta) * nEl.x + phiy[j](p_ksi, p_etta) * nEl.y);

				} 
				P_2_out[i][j] *= a;
			}
		}

		for(int i = 0; i < n_func_p; i++)
		{
			int id_i = element.dof_p[i];
			for(int j = 0; j < n_func_u; j++)
			{				
				int id_j = element.dof_u[j];
				A.add_element(id_i, id_j, P_2_out[i][j]); 
			}
		}
	}
	void OuterBoundaries::calculate_SP_out(int element_number, Matrix& A, EdgeSide edgeSide)
	{
		vector <vector<double>> SP_out;

		SP_out.resize(n_func_p);

		for (int i = 0; i < n_func_p; i++)
			initialize_vector(SP_out[i], n_func_p);

		Element element = elements[element_number];
		double lambda = calculate_lambda(element.number_of_area);

		double h, p_ksi, p_etta;
		Point nEl;

		if (edgeSide == LEFT || edgeSide == RIGHT)
		{
			h = get_hy(element_number);

			if (edgeSide == LEFT)
			{
				nEl.x = -1;
				p_ksi = -1;
			}
			else
			{
				nEl.x = 1;
				p_ksi = 1;
			}
			nEl.y = 0;
		}
		else
		{
			h = get_hx(element_number);

			if (edgeSide == LOW)
			{
				nEl.y = -1;
				p_etta = -1;
			}
			else
			{
				nEl.y = 1;
				p_etta = 1;
			}
			nEl.x = 0;
		}

		double jacobian = 0.5 * h;
		double st = jacobian * 1.0 / lambda * gamma;


		for(int i = 0; i < n_func_p; i++)
		{
			for(int j = 0; j < n_func_p; j++)
			{
				for(int k = 0; k < 3; k++)
				{
					if (edgeSide == LOW || edgeSide == UP)
						p_ksi = gauss_points_1[k];
					else
						p_etta = gauss_points_1[k];
					SP_out[i][j] += gauss_weights_1[k] * psi[i](p_ksi, p_etta) * psi[j](p_ksi, p_etta);

				} 
				SP_out[i][j] *= st;
			}
		}

		for(int i = 0; i < n_func_p; i++)
		{
			int id_i = element.dof_p[i];
			for(int j = 0; j < n_func_p; j++)
			{				
				int id_j = element.dof_p[j];
				A.add_element(id_i, id_j, SP_out[i][j]); 
			}
		}
	}

}
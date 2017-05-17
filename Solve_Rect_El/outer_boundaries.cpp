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
using namespace integration;

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
		sigma = jPenaltyParameters["sigma"].GetDouble();
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
			int id_i = element.dof_u[i];
			for(int j = 0; j < n_func_u; j++)
			{
				int id_j = element.dof_u[j];
				double e_outij = 0, s_outij = 0;
				for(int k = 0; k < n_ip1D; k++)
				{
					if (edgeSide == LOW || edgeSide == UP)
						p_ksi = gauss_points_1[k];
					else
						p_etta = gauss_points_1[k];

					e_outij += gauss_weights_1[k] *
								  (phix[j](p_ksi, p_etta) * nEl.x * dphixksi[i](p_ksi, p_etta) / hx +
								   phix[j](p_ksi, p_etta) * nEl.y * dphixetta[i](p_ksi, p_etta) / hy +
								   phiy[j](p_ksi, p_etta) * nEl.x * dphiyksi[i](p_ksi, p_etta) / hx +
								   phiy[j](p_ksi, p_etta) * nEl.y * dphiyetta[i](p_ksi, p_etta) / hy +
								   phix[i](p_ksi, p_etta) * nEl.x * dphixksi[j](p_ksi, p_etta) / hx +
								   phix[i](p_ksi, p_etta) * nEl.y * dphixetta[j](p_ksi, p_etta) / hy +
								   phiy[i](p_ksi, p_etta) * nEl.x * dphiyksi[j](p_ksi, p_etta) / hx +
								   phiy[i](p_ksi, p_etta) * nEl.y * dphiyetta[j](p_ksi, p_etta) / hy);

					s_outij += gauss_weights_1[k] *
								  (phix[j](p_ksi, p_etta) * phix[i](p_ksi, p_etta) * nEl2.x +
								   phix[j](p_ksi, p_etta) * phix[i](p_ksi, p_etta) * nEl2.y +
								   phiy[j](p_ksi, p_etta) * phiy[i](p_ksi, p_etta) * nEl2.x +
								   phiy[j](p_ksi, p_etta) * phiy[i](p_ksi, p_etta) * nEl2.y);
				} 
				e_outij *= a;
				s_outij *= st;
				A.add_element(id_i, id_j, -e_outij + s_outij);
			}
		}
	}
	void OuterBoundaries::calculate_P_1_out(int element_number, Matrix& A, EdgeSide edgeSide)
	{
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
			int id_i = element.dof_u[i];
			for(int j = 0; j < n_func_p; j++)
			{
				int id_j = element.dof_p[j];
				double p_1_outij = 0;
				for(int k = 0; k < n_ip1D; k++)
				{
					if (edgeSide == LOW || edgeSide == UP)
						p_ksi = gauss_points_1[k];
					else
						p_etta = gauss_points_1[k];

					p_1_outij += gauss_weights_1[k]  * psi[j](p_ksi, p_etta) *
									(phix[i](p_ksi, p_etta) * nEl.x + phiy[i](p_ksi, p_etta) * nEl.y);

				} 
				p_1_outij *= a;
				A.add_element(id_i, id_j, p_1_outij);
			}
		}
	}
	void OuterBoundaries::calculate_P_2_out(int element_number, Matrix& A, EdgeSide edgeSide)
	{
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
			int id_i = element.dof_p[i];
			for(int j = 0; j < n_func_u; j++)
			{
				int id_j = element.dof_u[j];
				double p_2_outij = 0;
				for(int k = 0; k < n_ip1D; k++)
				{
					if (edgeSide == LOW || edgeSide == UP)
						p_ksi = gauss_points_1[k];
					else
						p_etta = gauss_points_1[k];

					p_2_outij += gauss_weights_1[k]  * psi[i](p_ksi, p_etta) *
									(phix[j](p_ksi, p_etta) * nEl.x + phiy[j](p_ksi, p_etta) * nEl.y);

				} 
				p_2_outij *= a;
				A.add_element(id_i, id_j, p_2_outij);
			}
		}
	}
	void OuterBoundaries::calculate_SP_out(int element_number, Matrix& A, EdgeSide edgeSide)
	{
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
		double st = jacobian * 1.0 / lambda * sigma;


		for(int i = 0; i < n_func_p; i++)
		{
			int id_i = element.dof_p[i];
			for(int j = 0; j < n_func_p; j++)
			{
				int id_j = element.dof_p[j];
				double sp_outij = 0;
				for(int k = 0; k < n_ip1D; k++)
				{
					if (edgeSide == LOW || edgeSide == UP)
						p_ksi = gauss_points_1[k];
					else
						p_etta = gauss_points_1[k];
					sp_outij += gauss_weights_1[k] * psi[i](p_ksi, p_etta) * psi[j](p_ksi, p_etta);

				} 
				sp_outij *= st;
				A.add_element(id_i, id_j, sp_outij);
			}
		}
	}

}
#include "boundary_conditions.h"
#include "element.h"
#include "partition.h"
#include "myfunctions.h"
#include "point.h"

#include "rapidjson/document.h"
#include <fstream>

using namespace std;
using namespace element;
using namespace partition;
using namespace parameters;
using namespace myvector;
using namespace basis;
using namespace point;
using namespace integration;

using namespace rapidjson;

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

	void BoundaryConditionsSupport::initialize_penalty_parameters(std::string fileName)
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

	void BoundaryConditionsSupport::calculate_all_boundaries1(MyVector& b)
	{
		int size_b = boundaries1.size();
		for(int i = 0; i < size_b; i++)
			calculate_boundaries1(i, b);
	}

	void BoundaryConditionsSupport::calculate_boundaries1_for_side(int number, myvector::MyVector & b, EdgeSide edgeSide)
	{
		Element element = elements[boundaries1[number].elem];
		double hy = get_hy(boundaries1[number].elem);
		double hx = get_hx(boundaries1[number].elem);

		double x0 = nodes[element.nodes[0]].x;
		double y0 = nodes[element.nodes[0]].y;

		double lambda = calculate_lambda(element.number_of_area);

		double p_etta, p_y, p_ksi, p_x, h;
		int formula;
		Point nEl;
	

		if (edgeSide == LOW || edgeSide == UP)
		{
			if (edgeSide == LOW)
			{
				nEl.y = -1;
				p_etta = -1;
				p_y = y0;
				formula = boundaries1[number].formula_number[2];
			}
			else
			{
				nEl.y = 1;
				p_etta = 1;
				p_y = y0 + hy;
				formula = boundaries1[number].formula_number[3];
			}
			h = hx;
			nEl.x = 0;
		}
		else
		{
			if (edgeSide == LEFT)
			{
				nEl.x = -1;
				p_ksi = -1;
				p_x = x0;
				formula = boundaries1[number].formula_number[0];
			}
			else
			{
				nEl.x = 1;
				p_ksi = 1;
				p_x = x0 + hx;
				formula = boundaries1[number].formula_number[1];
			}
			h = hy;
			nEl.y = 0;
		}


		double jacobian = h * 0.5;
		double Ugi, Pgi;
		double g_x, g_y;
		int k = element.order;
		double c = gamma * k * k / h;
		double mu = lambda * c;

		for (int i = 0; i < n_func_u; i++)
		{
			Ugi = 0;
			for (int j = 0; j < n_ip1D; j++)
			{
				if (edgeSide == LOW || edgeSide == UP)
				{
					p_ksi = gauss_points_1[j];
					p_x = p_ksi * hx + x0;
				}
				else
				{
					p_etta = gauss_points_1[j];
					p_y = p_etta * hy + y0;
				}

				g_x = gx(formula, p_x, p_y);
				g_y = gy(formula, p_x, p_y);
				Ugi += -gauss_weights_1[j] *
					  (g_x * nEl.x * dphixksi[i](p_ksi, p_etta) / hx +
					   g_x * nEl.y * dphixetta[i](p_ksi, p_etta) / hy +
					   g_y * nEl.x * dphiyksi[i](p_ksi, p_etta) / hx +
					   g_y * nEl.y * dphiyetta[i](p_ksi, p_etta) / hy) +
					   gauss_weights_1[j] * mu * (g_x * phix[i](p_ksi, p_etta) +
					   g_y * phiy[i](p_ksi, p_etta));
			}
			Ugi *= jacobian;
			b[element.dof_u[i]] += Ugi;
		}

		for (int i = 0; i < n_func_p; i++)
		{
			Pgi = 0;
			for (int j = 0; j < n_ip1D; j++)
			{
				if (edgeSide == LOW || edgeSide == UP)
				{
					p_ksi = gauss_points_1[j];
					p_x = p_ksi * hx + x0;
				}
				else
				{
					p_etta = gauss_points_1[j];
					p_y = p_etta * hy + y0;
				}
				g_x = gx(formula, p_x, p_y);
				g_y = gy(formula, p_x, p_y);
				Pgi -= gauss_weights_1[j] * psi[i](p_ksi, p_etta) *
					  (g_x * nEl.x + g_y * nEl.y);
			}
			Pgi *= jacobian;
			b[element.dof_p[i]] += Pgi;
		}
	}

	void BoundaryConditionsSupport::calculate_boundaries1(int number, MyVector& b)
	{
		if(boundaries1[number].edges[0] == 1) calculate_boundaries1_for_side(number, b, LEFT);
		if(boundaries1[number].edges[1] == 1) calculate_boundaries1_for_side(number, b, RIGHT);
		if(boundaries1[number].edges[2] == 1) calculate_boundaries1_for_side(number, b, LOW);
		if(boundaries1[number].edges[3] == 1) calculate_boundaries1_for_side(number, b, UP);
	}

}

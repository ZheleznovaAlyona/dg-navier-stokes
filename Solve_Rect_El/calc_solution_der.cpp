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
	double MainSolver::get_solution_in_point_ux(double x, double y, int element_number, MyVector qi)
	{
		vector <int> indexes;
		indexes.resize(elements[element_number].ndof_u);

		vector <double> qi_local;
		qi_local.resize(elements[element_number].ndof_u);

		double u_in_point;

		//собираем глобальные номера с элемента
		for (int j = 0; j < elements[element_number].ndof_u; j++)
			indexes[j] = elements[element_number].dof_u[j];

		//собираем локальный набор весов
		for (int j = 0; j < elements[element_number].ndof_u; j++)
			qi_local[j] = qi[indexes[j]];

		//вычисляем в решение в точке
		static auto& part = static_cast<Partition>(*this);
		u_in_point = 0;
		for (int j = 0; j < elements[element_number].ndof_u; j++)
			u_in_point += qi_local[j] * phix_i(j, x, y, element_number, part);

		return u_in_point;
	}

	double MainSolver::get_solution_in_point_uy(double x, double y, int element_number, MyVector qi)
	{
		vector <int> indexes;
		indexes.resize(elements[element_number].ndof_u);

		vector <double> qi_local;
		qi_local.resize(elements[element_number].ndof_u);

		double u_in_point;

		//собираем глобальные номера с элемента
		for (int j = 0; j < elements[element_number].ndof_u; j++)
			indexes[j] = elements[element_number].dof_u[j];

		//собираем локальный набор весов
		for (int j = 0; j < elements[element_number].ndof_u; j++)
			qi_local[j] = qi[indexes[j]];

		//вычисляем в решение в точке
		static auto& part = static_cast<Partition>(*this);
		u_in_point = 0;
		for (int j = 0; j < elements[element_number].ndof_u; j++)
			u_in_point += qi_local[j] * phiy_i(j, x, y, element_number, part);

		return u_in_point;
	}

	double MainSolver::get_solution_in_point_uxdx(double x, double y, int element_number, MyVector qi)
	{
		vector <int> indexes;
		indexes.resize(elements[element_number].ndof_u);

		vector <double> qi_local;
		qi_local.resize(elements[element_number].ndof_u);

		double du_in_point;

		//собираем глобальные номера с элемента
		for (int j = 0; j < elements[element_number].ndof_u; j++)
			indexes[j] = elements[element_number].dof_u[j];

		//собираем локальный набор весов
		for (int j = 0; j < elements[element_number].ndof_u; j++)
			qi_local[j] = qi[indexes[j]];

		//вычисляем в решение в точке
		static auto& part = static_cast<Partition>(*this);
		du_in_point = 0;
		for (int j = 0; j < elements[element_number].ndof_u; j++)
			du_in_point += qi_local[j] * phixdx_i(j, x, y, element_number, part);

		return du_in_point;
	}

	double MainSolver::get_solution_in_point_uydy(double x, double y, int element_number, MyVector qi)
	{
		vector <int> indexes;
		indexes.resize(elements[element_number].ndof_u);

		vector <double> qi_local;
		qi_local.resize(elements[element_number].ndof_u);

		double du_in_point;

		//собираем глобальные номера с элемента
		for (int j = 0; j < elements[element_number].ndof_u; j++)
			indexes[j] = elements[element_number].dof_u[j];

		//собираем локальный набор весов
		for (int j = 0; j < elements[element_number].ndof_u; j++)
			qi_local[j] = qi[indexes[j]];

		//вычисляем в решение в точке
		static auto& part = static_cast<Partition>(*this);
		du_in_point = 0;
		for (int j = 0; j < elements[element_number].ndof_u; j++)
			du_in_point += qi_local[j] * phiydy_i(j, x, y, element_number, part);

		return du_in_point;
	}

	double MainSolver::get_solution_in_point_p(double x, double y, int element_number, MyVector qi)
	{
		vector <int> indexes;
		indexes.resize(elements[element_number].ndof_p);

		vector <double> qi_local;
		qi_local.resize(elements[element_number].ndof_p);

		double p_in_point;

		//собираем глобальные номера с элемента
		for (int j = 0; j < elements[element_number].ndof_p; j++)
			indexes[j] = elements[element_number].dof_p[j];

		//собираем локальный набор весов
		for (int j = 0; j < elements[element_number].ndof_p; j++)
			qi_local[j] = qi[indexes[j]];

		//вычисляем в решение в точке
		static auto& part = static_cast<Partition>(*this);
		p_in_point = 0;
		for (int j = 0; j < elements[element_number].ndof_p; j++)
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
		indexes.resize(elements[0].ndof_u);

		vector <double> qi_local;
		qi_local.resize(elements[0].ndof_u);

		static auto& part = static_cast<Partition>(*this);

		for (int i = 0; i < size; i++)
		{
			//собираем глобальные номера с элемента
			for (int j = 0; j < elements[i].ndof_u; j++)
				indexes[j] = elements[i].dof_u[j];

			for (int j = 0; j < 4; j++)
				indexes_nodes[j] = elements[i].nodes[j];

			//собираем локальный набор весов
			for (int j = 0; j < elements[i].ndof_u; j++)
				qi_local[j] = qi[indexes[j]];

			//вычисляем в узлах элемента решение
			for (int j = 0; j < 4; j++)
			{
				u_local[j] = 0;
				x = nodes[indexes_nodes[j]].x;
				y = nodes[indexes_nodes[j]].y;
				for (int k = 0; k < elements[i].ndof_u; k++)
					u_local[j] += qi_local[k] * phix_i(k, x, y, i, part);
			}

			//кладём в результирующий вектор
			for (int j = 0; j < 4; j++)
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
		indexes.resize(elements[0].ndof_u);

		vector <double> qi_local;
		qi_local.resize(elements[0].ndof_u);

		static auto& part = static_cast<Partition>(*this);

		for (int i = 0; i < size; i++)
		{
			//собираем глобальные номера с элемента
			for (int j = 0; j < elements[i].ndof_u; j++)
				indexes[j] = elements[i].dof_u[j];

			//собираем локальный набор весов
			for (int j = 0; j < elements[i].ndof_u; j++)
				qi_local[j] = qi[indexes[j]];

			for (int j = 0; j < 4; j++)
				indexes_nodes[j] = elements[i].nodes[j];

			//вычисляем в узлах элемента решение
			for (int j = 0; j < 4; j++)
			{
				u_local[j] = 0;
				x = nodes[indexes_nodes[j]].x;
				y = nodes[indexes_nodes[j]].y;
				for (int k = 0; k < elements[i].ndof_u; k++)
					u_local[j] += qi_local[k] * phiy_i(k, x, y, i, part);
			}

			//кладём в результирующий вектор
			for (int j = 0; j < 4; j++)
				solution[indexes_nodes[j]] = u_local[j];
		}
	}

	void MainSolver::get_vector_solution_in_nodes_p(MyVector qi, MyVector &solution)
	{
		logger.send_message_P();
		int indexes_nodes[4];
		int size = elements.size();
		double p_local[4];
		double x, y;

		vector <int> indexes;
		indexes.resize(elements[0].ndof_p);

		vector <double> qi_local;
		qi_local.resize(elements[0].ndof_p);

		static auto& part = static_cast<Partition>(*this);

		for (int i = 0; i < size; i++)
		{
			//собираем глобальные номера с элемента
			for (int j = 0; j < elements[i].ndof_p; j++)
				indexes[j] = elements[i].dof_p[j];

			//собираем локальный набор весов
			for (int j = 0; j < elements[i].ndof_p; j++)
				qi_local[j] = qi[indexes[j]];

			for (int j = 0; j < 4; j++)
				indexes_nodes[j] = elements[i].nodes[j];

			//вычисляем в узлах элемента решение
			for (int j = 0; j < 4; j++)
			{
				p_local[j] = 0;
				x = nodes[indexes_nodes[j]].x;
				y = nodes[indexes_nodes[j]].y;
				for (int k = 0; k < elements[i].ndof_p; k++)
					p_local[j] += qi_local[k] * psi_i(k, x, y, i, part);
			}

			//кладём в результирующий вектор
			for (int j = 0; j < 4; j++)
				solution[indexes_nodes[j]] = p_local[j];
		}
	}
}
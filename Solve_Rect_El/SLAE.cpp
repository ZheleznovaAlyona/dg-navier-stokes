#include "SLAE.h"

using namespace boundary_conditions;
using namespace element;
using namespace partition;
using namespace point;
using namespace myvector;
using namespace matrix;
using namespace densematrix;
using namespace testingparameters;

namespace slae
{
void SLAE::initialize(int max_number_of_iterations,
					  int max_number_of_iterations_non_lin,
					  double epsilon,
					  int gmres_m,
					  ifstream& grid_f_in,
					  ifstream& elements_f_in,
					  string log_f,
					  ifstream& boundary1)
{
	boundary1 >> boundaries1;

	int tmp;

	input(grid_f_in, elements_f_in);

	n = nodes.size() + elements.size() * 4; //узлы и рёбра
	max_iter = max_number_of_iterations;
	max_iter_nonlinear = max_number_of_iterations_non_lin;
	eps = epsilon;
	m = gmres_m;

	//sigma = 1; 
	//mu1 
	//mu2 = 1;
	logger.open(log_f);

	tmp = count_unzero_matrix_elements();

	A.initialize(n, tmp);

	Ux_numerical.initialize(nodes.size());
	Uy_numerical.initialize(nodes.size());
	P_numerical.initialize(nodes.size());
	q_prev.initialize(n);

	A.LU_ggl.reserve(tmp); A.LU_ggu.reserve(tmp);
	for(int i = 0; i < tmp; i++)
	{
		A.LU_ggl.push_back(0.0);
		A.LU_ggu.push_back(0.0);
	}

	A.LU_di.reserve(n);
	for(int i = 0; i < n; i++)
		A.LU_di.push_back(0.0);

	A.yl.initialize(n); A.yu.initialize(n);

	//инициализация базиса

	//инициализация параметров интегрирования
}

void SLAE::reinitialize()
{
	Ux_numerical.make_zero();
	Uy_numerical.make_zero();
	P_numerical.make_zero();
	A.reinitialize();

	int gg_size = A.LU_ggl.size();

	for(int i = 0; i < gg_size; i++)
	{
		A.LU_ggl[i] = 0.0;
		A.LU_ggu[i] = 0.0;
	}

	for(int i = 0; i < n; i++)
		A.LU_di[i] = 0.0;

	if(Testing_parameters::instance().solver == 0) 
	{
		int gg_size = LU_ggl2.size();

		for(int i = 0; i < gg_size; i++)
		{
			LU_ggl2[i] = 0.0;
			LU_ggu2[i] = 0.0;
		}

		for(int i = 0; i < n; i++)
			LU_di2[i] = 0.0;
	}
}



int SLAE::count_unzero_matrix_elements()
{
	int count_uu = 0;
	for(unsigned int i = 0; i < elements.size(); i++)
	{
		//с собой
		count_uu += 4;
		//с соседями
		for(int j = 0; j < 4; j++)
			if(elements[i].neighbors[j] != -1)
				count_uu += 4;
	}
	count_uu *= 4;
	count_uu -= elements.size() * 4;
	//так как нужно для одного треугольника ввиду симметричности портрета,
	//то необходимо полученное количество поделить на 2
	count_uu /= 2;

	int count_pp = 0;
	for(unsigned int i = 0; i < elements.size(); i++)
	{
		//с собой
		count_pp += 4;
		//с соседями
		for(int j = 0; j < 4; j++)
			if(elements[i].neighbors[j] != -1)
				count_pp += 4;
	}
	count_pp *= 4;
	count_pp -= nodes.size();
	//так как нужно для одного треугольника ввиду симметричности портрета,
	//то необходимо полученное количество поделить на 2
	count_pp /= 2;

	int count_up = 0;
	for(unsigned int i = 0; i < elements.size(); i++)
	{
		//с собой
		count_up += 4;
		//с соседями
		for(int j = 0; j < 4; j++)
			if(elements[i].neighbors[j] != -1)
				count_up += 4;
	}
	count_up *= 4;

	return count_uu + count_pp + count_up;
}

int SLAE::create_unzero_elements_list(int element_number,
									  vector <int> &list,
									  int dof_num_i,
									  int dof_num_j,
									  int *dof_i,
									  int *dof_j,
									  bool dof_j_edge)
{
	int neighbor;

	//свои по строкам
	for(int i = 0; i < dof_num_i; i++)
		list.push_back(dof_i[i]);

	//свои
	for(int i = 0; i < dof_num_j; i++)
		list.push_back(dof_j[i]);

	//соседей
	for(int j = 0; j < 4; j++)
	{
		neighbor = elements[element_number].neighbors[j];
		if(neighbor != -1)
		for(int i = 0; i < dof_num_j; i++)
		{
			if(dof_j_edge) list.push_back(elements[neighbor].edges[i]);
			else list.push_back(elements[neighbor].nodes[i]);
		}			
	}

	return list.size();
}

#pragma region локальные матрицы и векторы

void SLAE::calculate_locals(int element_number, MyVector q_calc)
{
	Element element = elements[element_number];
	int id_i, id_j;
	int n_edges = elements.size() * 4;
	
	calculate_G(element_number);
	calculate_C(element_number, q_calc);
	calculate_P1(element_number);
	calculate_P2(element_number);
	calculate_F(element_number);

	for(int i = 0; i < 4; i++)
	{
		id_i = element.edges[i];
		for(int j = 0; j < 4; j++)
		{				
			id_j = element.edges[j];
			A.add_element(id_i, id_j, G[i][j] + C[i][j]); 
		}
		A.b[id_i] += F[i];
	}

	for(int i = 0; i < 4; i++)
	{
		id_i = element.edges[i];
		for(int j = 0; j < 4; j++)
		{				
			id_j = element.nodes[j] + n_edges;
			add_element_to_global_matrix(id_i, id_j, P1[i][j]); 
		}
	}

	for(int i = 0; i < 4; i++)
	{
		id_i = element.nodes[i] + n_edges;
		for(int j = 0; j < 4; j++)
		{				
			id_j = element.edges[j];
			add_element_to_global_matrix(id_i, id_j, P2[i][j]); 
		}
	}
}

void SLAE::calculate_G(int element_number)
{
	double hx, hy, hx2, hy2, g1, g2;
	Element element = elements[element_number];
	double lambda = calculate_lambda(element.number_of_area);

	hx = get_hx(element_number);
	hy = get_hy(element_number);
	hx2 = hx * hx;
	hy2 = hy * hy;

	double jacobian = hx * hy / 4.0;

	for(int i = 0; i < 4; i++)
	{
		for(int j = i; j < 4; j++)
		{
			g1 = 0;
			g2 = 0;
			for(int k = 0; k < 9; k++)
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

	for(int i = 1; i < 4; i++)
		for(int j = 0; j < i; j++)
			G[i][j] = G[j][i];
}
void SLAE::calculate_C(int element_number, MyVector q_calc)
{
	double hx, hy, c1, c2;
	Element element = elements[element_number];

	hx = get_hx(element_number);
	hy = get_hy(element_number);

	double jacobian = hx * hy / 4.0;
	double u_x, u_y;
	double x0 = nodes[element.nodes[0]].x;
	double y0 = nodes[element.nodes[0]].y;

	for(int i = 0; i < 4; i++)
	{
		for(int j = 0; j < 4; j++)
		{
			C[i][j] = 0;
			for(int k = 0; k < 9; k++)
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
						    dphiyksi[j](p_ksi, p_etta) / hx * phiy[i](p_ksi, p_etta)-
							dphixksi[i](p_ksi, p_etta) / hx * phix[j](p_ksi, p_etta) -
						    dphiyksi[i](p_ksi, p_etta) / hx * phiy[j](p_ksi, p_etta));
				c2 = u_y * (dphixetta[j](p_ksi, p_etta) / hy * phix[i](p_ksi, p_etta) +
						    dphiyetta[j](p_ksi, p_etta) / hy * phiy[i](p_ksi, p_etta)-
							dphixetta[i](p_ksi, p_etta) / hy * phix[j](p_ksi, p_etta) -
						    dphiyetta[i](p_ksi, p_etta) / hy * phiy[j](p_ksi, p_etta));			
				C[i][j] += gauss_weights[k] * 0.5 * (c1 + c2);
			}
			C[i][j] *= jacobian;
		}
	}
}
void SLAE::calculate_P1(int element_number)
{
	Element element = elements[element_number];
	double hx = get_hx(element_number);
	double hy = get_hy(element_number);

	double jacobian = hx * hy / 4.0;

	double rho = calculate_rho(element.number_of_area);
	rho = - 1 / rho;

	for(int i = 0; i < 4; i++)
	{
		for(int j = 0; j < 4; j++)
		{
			P1[i][j] = 0;
			for(int k = 0; k < 9; k++)
			{
				double p_ksi = gauss_points[0][k], p_etta = gauss_points[1][k];
				P1[i][j] += gauss_weights[k] * psi[j](p_ksi, p_etta) *
						   	(dphixksi[i](p_ksi, p_etta) / hx + 
							dphiyetta[i](p_ksi, p_etta) / hy);
			}
			P1[i][j] *= jacobian * rho;
		}
	}
}
void SLAE::calculate_P2(int element_number)
{
	Element element = elements[element_number];
	double hx = get_hx(element_number);
	double hy = get_hy(element_number);

	double jacobian = hx * hy / 4.0;

	for(int i = 0; i < 4; i++)
	{
		for(int j = 0; j < 4; j++)
		{
			P2[i][j] = 0;
			for(int k = 0; k < 9; k++)
			{
				double p_ksi = gauss_points[0][k], p_etta = gauss_points[1][k];
				P2[i][j] += gauss_weights[k] * psi[i](p_ksi, p_etta) *
						   	(dphixksi[j](p_ksi, p_etta) / hx + 
							dphiyetta[j](p_ksi, p_etta) / hy);
			}
			P2[i][j] *= jacobian;
		}
	}
}

void SLAE::calculate_F(int element_number)
{
	double f_x, f_y;
	Element element = elements[element_number];
	int number_of_area = element.number_of_area;

	double x0 = nodes[element.nodes[0]].x;
	double y0 = nodes[element.nodes[0]].y;

	double hx = get_hx(element_number);
	double hy = get_hy(element_number);

	double jacobian = hx * hy / 4.0;

	for(int i = 0;  i < 4; i++)
	{
		F[i] = 0;
		for(int j = 0; j < 9; j++)
		{
			double p_ksi = gauss_points[0][j], p_etta = gauss_points[1][j];
			double p_x = p_ksi * hx + x0, p_y = p_etta * hy + y0;
			f_x = calculate_fx(number_of_area, p_x, p_y);	
			f_y = calculate_fy(number_of_area, p_x, p_y);
			F[i] += gauss_weights[j] * (f_x * phix[i](p_ksi, p_etta) 
					+ f_y * phiy[i](p_ksi, p_etta));
		}
		 F[i] *= jacobian;
	}
}

#pragma endregion



#pragma region для глобальной матрицы

void SLAE::create_portret()
{
	vector <int> unzero_elements_list;
	vector <int> *lists;	
	int unzero_elements_lists_size; 
	int current_number;
	int n_edges = elements.size() * 4;

	lists = new vector <int>[n];
	int count_elements = elements.size();

	for(int i = 0; i < count_elements; i++)
	{
		//общий принцип сборки портрета для uu, pp, pu
		//--------------------------------------------
		//1. собираем в список:
		//*для к.э. глобальные номера dof по i и затем
		//*ненулевые для к.э. глобальные номера dof по j;
		//2. идём по первым элементам списка (dof по i)
		//и выбираем для каждого номера (для p номер=индекс цикла + число рёбер),
		//меньшие его, т.е. те, которые будут располагаться левее соответствующей 
		//диагонали, потому что портрет строим по строкам
		//а затем кладём в соответствующий список
		//2.a. сортируем список по возрастанию

		//структура СЛАУ:
		/*
		 ---------
		| UU | UP |
		-----|----
		| PU | PP |
		 ---------
		*/

		//блок UU
		//1
		unzero_elements_lists_size = create_unzero_elements_list(i, 
																 unzero_elements_list, 
																 4, 
																 4, 
																 elements[i].edges, 
																 elements[i].edges, 
																 true);
		//2
		for(int j = 0; j < 4; j++)
		{
			current_number = unzero_elements_list[j];
			for(int k = 4; k < unzero_elements_lists_size; k++)
				if(unzero_elements_list[k] < current_number)
					lists[current_number].push_back(unzero_elements_list[k]);
			    // 2.a
				sort(lists[current_number].begin(), lists[current_number].end());
		}
		unzero_elements_list.clear();

		//блок PP
		//1
		unzero_elements_lists_size = create_unzero_elements_list(i, 
																 unzero_elements_list, 
																 4, 
																 4, 
																 elements[i].nodes, 
																 elements[i].nodes, 
																 false);
		//2
		for(int j = 0; j < 4; j++)
		{
			current_number = unzero_elements_list[j];
			for(int k = 4; k < unzero_elements_lists_size; k++)
				if(unzero_elements_list[k] < current_number)
					lists[current_number + n_edges].push_back(unzero_elements_list[k] + n_edges);	
			    //2.a 
				//можно не сортировать, потому что потом всё равно добавятся ещё элементы из PU
		}
		unzero_elements_list.clear();

		//блок PU
		//1
		unzero_elements_lists_size = create_unzero_elements_list(i,
																 unzero_elements_list, 
																 4, 
																 4, 
																 elements[i].nodes, 
																 elements[i].edges, 
																 true);
		//2
		for(int j = 0; j < 4; j++)
		{
			current_number = unzero_elements_list[j];
			for(int k = 4; k < unzero_elements_lists_size; k++)
				if(unzero_elements_list[k] < current_number + n_edges)
					lists[current_number + n_edges].push_back(unzero_elements_list[k]);	
			    // 2.a
				sort(lists[current_number + n_edges].begin(), 
					lists[current_number + n_edges].end());
		}
		unzero_elements_list.clear();
	}

	A.ig[0] = 0;

	for(int i = 0; i < n; i++)
	{
		if(!lists[i].empty())
			A.ig[i + 1] = A.ig[i] + lists[i].size();
		else A.ig[i + 1] = A.ig[i];
	}

	int k = 0;
	for(int i = 0; i < n; i++)
	{
		if(!lists[i].empty())
		{
			for(unsigned int j = 0; j < lists[i].size(); j++)
			{
				A.jg[k] = lists[i][j];
				k++;
			}
			lists[i].clear();
		}
	}

	delete[] lists;
}

void SLAE::calculate_global_matrix(MyVector q_calc)
{
	int size = elements.size();

	//локаьные матрицы и вектор правой части
	for(int el_i = 0; el_i < size; el_i++)
		calculate_locals(el_i, q_calc);

	//матрицы межэлементных границ
	for(int el_i = 0; el_i < size; el_i++)
		calculate_internal_boundaries(el_i, A);

	//расчёт матриц границ области
		for(int el_i = 0; el_i < size; el_i++)
			calculate_outer_boundaries(el_i, A);

	//учёт первых краевых условий
	calculate_all_boundaries1(A.b);
}

#pragma endregion

double SLAE::get_solution_in_point_ux(double x, double y, int element_number, MyVector qi)
{
	int indexes[4];
	double u_in_point, qi_local[4];

	//собираем глобальные номера с элемента
	for(int j = 0; j < 4; j++)
		indexes[j] = elements[element_number].edges[j];

	//собираем локальный набор весов
	for(int j = 0; j < 4; j++)
		qi_local[j] = qi[indexes[j]];

	//вычисляем в решение в точке
	u_in_point = 0;
	for(int j = 0; j < 4; j++)
		u_in_point += qi_local[j] * phix_i(j, x, y, element_number);

	return u_in_point;
}

double SLAE::get_solution_in_point_uy(double x, double y, int element_number, MyVector qi)
{
	int indexes[4];
	double u_in_point, qi_local[4];

	//собираем глобальные номера с элемента
	for(int j = 0; j < 4; j++)
		indexes[j] = elements[element_number].edges[j];

	//собираем локальный набор весов
	for(int j = 0; j < 4; j++)
		qi_local[j] = qi[indexes[j]];

	//вычисляем в решение в точке
	u_in_point = 0;
	for(int j = 0; j < 4; j++)
		u_in_point += qi_local[j] * phiy_i(j, x, y, element_number);

	return u_in_point;
}

double SLAE::get_solution_in_point_uxdx(double x, double y, int element_number, MyVector qi)
{
	int indexes[4];
	double du_in_point, qi_local[4];

	//собираем глобальные номера с элемента
	for(int j = 0; j < 4; j++)
		indexes[j] = elements[element_number].edges[j];

	//собираем локальный набор весов
	for(int j = 0; j < 4; j++)
		qi_local[j] = qi[indexes[j]];

	//вычисляем в решение в точке
	du_in_point = 0;
	for(int j = 0; j < 4; j++)
		du_in_point += qi_local[j] * phixdx_i(j, x, y, element_number);

	return du_in_point;
}

double SLAE::get_solution_in_point_uydy(double x, double y, int element_number, MyVector qi)
{
	int indexes[4];
	double du_in_point, qi_local[4];

	//собираем глобальные номера с элемента
	for(int j = 0; j < 4; j++)
		indexes[j] = elements[element_number].edges[j];

	//собираем локальный набор весов
	for(int j = 0; j < 4; j++)
		qi_local[j] = qi[indexes[j]];

	//вычисляем в решение в точке
	du_in_point = 0;
	for(int j = 0; j < 4; j++)
		du_in_point += qi_local[j] * phiydy_i(j, x, y, element_number);
	
	return du_in_point;
}

double SLAE::get_solution_in_point_p(double x, double y, int element_number, MyVector qi)
{
	int indexes[4], n_edges = elements.size() * 4;
	double p_in_point, qi_local[4];

	//собираем глобальные номера с элемента
	for(int j = 0; j < 4; j++)
		indexes[j] = elements[element_number].nodes[j] + n_edges;

	//собираем локальный набор весов
	for(int j = 0; j < 4; j++)
		qi_local[j] = qi[indexes[j]];

	//вычисляем в решение в точке
	p_in_point = 0;
	for(int j = 0; j < 4; j++)
		p_in_point += qi_local[j] * psi_i(j, x, y, element_number);

	return p_in_point;
}

double SLAE::get_solution_in_point2_ux(double x, double y, MyVector qi)
{
	int element_number = search_element(x, y);
	return get_solution_in_point_ux(x, y, element_number, qi);
}

double SLAE::get_solution_in_point2_uy(double x, double y, MyVector qi)
{
	int element_number = search_element(x, y);
	return get_solution_in_point_uy(x, y, element_number, qi);
}

double SLAE::get_solution_in_point2_p(double x, double y, MyVector qi)
{
	int element_number = search_element(x, y);
	return get_solution_in_point_p(x, y, element_number, qi);
}

int SLAE::search_element(double x, double y)
{
	double x_left, x_right, y_low, y_up;
	int size = elements.size();

	for(int i = 0; i < size; i++)
	{
		x_left = nodes[elements[i].nodes[0]].x;
		x_right = nodes[elements[i].nodes[1]].x;
		y_low = nodes[elements[i].nodes[0]].y;
		y_up = nodes[elements[i].nodes[3]].y;
		if(x_left <= x && x <= x_right && y_low <= y && y <= y_up)
			return i;
	}
}

void SLAE::get_vector_solution_in_nodes_ux(MyVector qi, MyVector &solution)
{
	int indexes[4], indexes_nodes[4];
	int size = elements.size();
	double u_local[4], qi_local[4];
	double x, y;

	for(int i = 0; i < size; i++)
	{
		//собираем глобальные номера с элемента
		for(int j = 0; j < 4; j++)
			indexes[j] = elements[i].edges[j];

		for(int j = 0; j < 4; j++)
			indexes_nodes[j] = elements[i].nodes[j];

		//собираем локальный набор весов
		for(int j = 0; j < 4; j++)
			qi_local[j] = qi[indexes[j]];

		//вычисляем в узлах элемента решение
		for(int j = 0; j < 4; j++)
		{
			u_local[j] = 0;
			x = nodes[indexes_nodes[j]].x;
			y = nodes[indexes_nodes[j]].y;
			for(int k = 0; k < 4; k++)
				u_local[j] += qi_local[k] * phix_i(k, x, y, i);
		}

		//кладём в результирующий вектор
		for(int j = 0; j < 4; j++)
			solution[indexes_nodes[j]] = u_local[j];
	}
}

void SLAE::get_vector_solution_in_nodes_uy(MyVector qi, MyVector &solution)
{
	int indexes[4], indexes_nodes[4];
	int size = elements.size();
	double u_local[4], qi_local[4];
	double x, y;

	for(int i = 0; i < size; i++)
	{
		//собираем глобальные номера с элемента
		for(int j = 0; j < 4; j++)
			indexes[j] = elements[i].edges[j];

		//собираем локальный набор весов
		for(int j = 0; j < 4; j++)
			qi_local[j] = qi[indexes[j]];

		for(int j = 0; j < 4; j++)
			indexes_nodes[j] = elements[i].nodes[j];

		//вычисляем в узлах элемента решение
		for(int j = 0; j < 4; j++)
		{
			u_local[j] = 0;
			x = nodes[indexes_nodes[j]].x;
			y = nodes[indexes_nodes[j]].y;
			for(int k = 0; k < 4; k++)
				u_local[j] += qi_local[k] * phiy_i(k, x, y, i);
		}

		//кладём в результирующий вектор
		for(int j = 0; j < 4; j++)
			solution[indexes_nodes[j]] = u_local[j];
	}
}

void SLAE::get_vector_solution_in_nodes_p(MyVector qi, MyVector &solution)
{
	int indexes[4], n_edges = elements.size() * 4;
	int size = elements.size();
	double p_local[4], qi_local[4];
	double x, y;

	for(int i = 0; i < size; i++)
	{
		//собираем глобальные номера с элемента
		for(int j = 0; j < 4; j++)
			indexes[j] = elements[i].nodes[j];

		//собираем локальный набор весов
		for(int j = 0; j < 4; j++)
			qi_local[j] = qi[indexes[j]+ n_edges];

		//вычисляем в узлах элемента решение
		for(int j = 0; j < 4; j++)
		{
			p_local[j] = 0;
			x = nodes[indexes[j]].x;
			y = nodes[indexes[j]].y;
			for(int k = 0; k < 4; k++)
				p_local[j] += qi_local[k] * psi_i(k, x, y, i);
		}

		//кладём в результирующий вектор
		for(int j = 0; j < 4; j++)
			solution[indexes[j]] = p_local[j];
	}
}


#pragma region решатель

void SLAE::solve_min_sqr_problem(MyVector d, DenseMatrix H, MyVector &result)
{
	int m2 = H.n_columns;
	DenseMatrix H_previous(m2 + 1, m2), H2(m, m);
	MyVector d1(m2 + 1), d2(m);
	double ci, si, tmp;

	double *tmp1, *tmp2;
	tmp1 = new double[m];
	tmp2 = new double[m];
	double tmp11, tmp22;

	H_previous = H;
	d1 = d;
	for(int i = 0; i < m; i++)
	{
		tmp = sqrt(H_previous[i][i] * H_previous[i][i] +
									H[i][i + 1] * H[i][i + 1]);
		ci = H_previous[i][i] / tmp;
		si = H[i][i + 1] / tmp;

		#pragma region H_prev=R*H_prev
		//расчитываем заранее элементы строк, где блок синусов-косинусов
		for(int l = 0; l < m; l++)
		{
			tmp1[l] = H_previous[l][i] * ci + H_previous[l][i + 1] * si;
			tmp2[l] = -H_previous[l][i] * si + H_previous[l][i + 1] * ci;
		}

		//заполняем строки,где блок синусов-косинусов
		for(int l = 0; l < m; l++)
		{

			H_previous[l][i] = tmp1[l];
			H_previous[l][i + 1] = tmp2[l];
		}

		#pragma endregion

		#pragma region d1=R*d1
		tmp11 = 0;
		//рассчитываем элементы вектора, где блок синусов-косинусов
		tmp11 = d1[i] * ci + d1[i + 1] * si;
		tmp22 = -d1[i] * si + d1[i + 1] * ci;

		//заполняем элементы вектора, где блок синусов-косинусов
		//остальные не изменяются
		d1[i] = tmp11;
		d1[i + 1] = tmp22;
		#pragma endregion
	}

	//исключаем m+1-ю строку и m+1-й элемент
	for(int j = 0; j < m; j++)
	{
		for(int i = 0; i < m; i++)
			H2[j][i] = H_previous[j][i];
		d2[j] = d1[j];
	}

	//находим неизвестный вектор из СЛАУ H2*result=d2
	for(int i = m-1; i >= 0; i--)
	{
		result[i] = d2[i];
		for(int j = i + 1; j < m; j++)
		{		
			result[i] -= result[j] * H2[j][i];
		}
		result[i] = result[i] / H2[i][i];
	}

}

void SLAE::GMRES(MyVector U_begin, MyVector &solution)
{
	MyVector r(n), x(n), d(m + 1), z(m), w(n), tmp(n), f(n), rr2(n);
	DenseMatrix V(n, m), H(m + 1, m);

	double norm_r, norm_f;
	bool continue_;

	x = U_begin;
	f = A.b;
	A.LU();
	r = f - A * x;
	A.LYF(r); r = A.yl;
	norm_r = r.norm();
	A.LYF(f); norm_f = A.yl.norm();
	x = A.Uv(U_begin);

	for(int k_iter = 1; k_iter <= max_iter && norm_r / norm_f > eps; k_iter++)
	{
		d.make_zero();
		V[0] = r / norm_r;

		continue_ = true;
		for(int j = 1; j <= m && continue_; j++)
		{
			A.UXY(V[j - 1]);
			tmp = A * A.yu;
			A.LYF(tmp); w = A.yl;

			for(int l = 1; l <= j; l++)
			{
				H[j - 1][l - 1] = scal(V[l - 1], w);
				w = w - V[l - 1] * H[j - 1][l - 1];
			}

			H[j - 1][j] = w.norm();
			if(abs(H[j - 1][j]) < 1E-14)
			{
				m = j;
				continue_ = false;
			}
			else
			{
				if(m != j)
					V[j] = w / H[j - 1][j];
			}
		}

		d[0] = norm_r;
		solve_min_sqr_problem(d, H, z);
		x = x + V * z;
		A.UXY(x);
		tmp = f - A * A.yu;
		A.LYF(tmp); r = A.yl;
		norm_r = r.norm();
		logger.send_current_information(norm_r / norm_f, k_iter);
		printf("%d\tr=%.10e\n", k_iter, norm_r / norm_f);
	}
	A.UXY(x); x = A.yu;
	solution = x;
}

void SLAE::BCGStab(MyVector U_begin, MyVector &solution)
{
	double  rkr0, ak, gk, bk;
	MyVector r(n), f(n), x(n), r0(n), z(n), p(n), v(n), v1(n), rr2(n);
	double r_norm, f_norm;

	A.LU();

	x = U_begin;
	f = A.b;
	f_norm = f.norm();
	r0 = f - A * x;
	A.LYF(r0); r0 = A.yl;
	r_norm = r0.norm() / f_norm;

	A.UXY(r0); z = A.yu;
	r = r0;

	logger.send_current_information(r_norm, 0);

	for(int k_it = 1; k_it <= max_iter && r_norm > eps; k_it++)
	{
		//найдем L^(-1)AU^(-1)zk
		A.UXY(z); v = A.yu;// v = U(-1)zk
		v1 = A * v; // v1 = AU^(-1)zk

		A.LYF(v1); v = A.yl;// v = L^(-1)AU^(-1)zk

		rkr0 = scal(r, r0);
		ak = rkr0 / scal(v, r0); // ak = (r,r0)/ ( L^(-1)AU^(-1)zk,r0)

		p = r - v * ak; // pk = r - ak*L^(-1)AU^(-1)zk

		//найдем L^(-1)AU^(-1)pk
		A.UXY(p); v1 = A.yu; // v1 = U^(-1)pk
		A.LYF(A * v1); v1 = A.yl; // v1 = L^(-1)AU^(-1)pk

		gk = scal(v1, p) / scal(v1, v1); // gk = (L^(-1)AU^(-1)pk,pk)/ ( L^(-1)AU^(-1)pk,L^(-1)AU^(-1)pk)

		//x = x + ak * zk + gk * pk;
		x = x + z * ak + p * gk;

		r = p - v1 * gk; //r = pk - gk * v1;

		bk = scal(r, r0) / rkr0  * (ak / gk);

		//zk = rk + bk * zk - gk * bk * v;
		z = r + z * bk - v * (gk * bk);

		r_norm = r.norm() / f_norm;
		printf("%d\tr=%.10e\n", k_it, r_norm);
		logger.send_current_information(r_norm, k_it);
	}
	A.UXY(x); x = A.yu;
	solution = x;
}

void SLAE::BCG(MyVector U_begin, MyVector &solution)
{
    MyVector r(n), r_(n), p(n), p_(n), f(n);
    MyVector v1(n), v2(n), v3(n), x(n);
    double alpha, betta, old_r_norm = 1.e+20, sc1, sc2;
	double r_norm, r_norm_;
	int k_it;

	A.LU();

	f = A.b;

    k_it = 0;
    r_norm = old_r_norm / 10;
    r_norm_ = 0;

	x = U_begin;
	v1 = A * x;
	v1 = f - v1;
    
	A.LYF(v1); r = A.yl;
	r_ = r;
	p = r;
	p_ = r_;
	int flag = 0;

    while(flag == 0 && k_it < max_iter)
    {
        sc1 = scal(r, r_);
		A.UXY(p); v1 = A.yu;
		v2 = A * v1;
		A.LYF(v2); v3 = A.yl;

        sc2 = scal(p_, v3);

        alpha = sc1 / sc2;
		x = x +  v1 * alpha;
		r = r - v3 * alpha;

		A.LYFt(p_); v1 = A.yl;

		v2 = A  /  v1;

		A.UXYt(v2); v3 = A.yu;

		r_ = r_ - v3 * alpha;
        sc2 = scal(r, r_);

        betta = sc2 / sc1;

        old_r_norm = r_norm;
		r_norm = r.norm();
		printf("%d\tr=%.10e\n", k_it, r_norm);

		logger.send_current_information(r_norm, k_it);
        if( r_norm < eps) flag = 1;

        if(flag == 0)
        {
			p = r + p * betta;
			p_ = r_ + p_ * betta;

            k_it++;
        }
    }

    printf ("\nk_iterations: %ld	r:%.6e", k_it, r_norm);
    if(flag == 0) printf ("\nexit r\n");
    if(flag == 1) printf ("\nsolution end\n");

    if(flag == 2) printf ("\nsolution not end!\n- change vector R_\n");
	solution = x;

}

void SLAE::Solve(MyVector U_begin, double &normL2u, double &normL2p)
{
	MyVector q(n);

	switch(Testing_parameters::instance().solver)
	{
	case 1:
		{
			//use_LU = false;
			BCGStab(U_begin, q);
		}
		break;
	case 2:
		{
			//use_LU = false;
			GMRES(U_begin, q);
		}
		break;
	case 3:
		{
			//use_LU = false;
			double eps2 = eps;
			eps = 1e-7;
			GMRES(U_begin, q);
			eps = eps2;
			BCG(q, q);
		}
		break;
	};

	q_prev = q;

}

#pragma endregion

#pragma region нелинейная

void SLAE::si_print(ofstream& log_f, int iteration_number, double &normL2u, double &normL2p)
{
	string f_name_s, f_name_i;

	log_f << "---" << iteration_number << "---" << endl;
	f_name_s = string("s_") + to_string(iteration_number) + ".txt";
	f_name_i = string("i_") + to_string(iteration_number) + ".txt";
	ofstream solution_f_out(f_name_s), info_f_out(f_name_i);

	output(solution_f_out, info_f_out, normL2u, normL2p);
	solution_f_out.close();
	info_f_out.close();
}

double SLAE::find_relaxation_parameter(MyVector q_current, MyVector q_previous, double &residual_previous)
{
	MyVector qNonL(n);
	double w = 1, residual;

	do
	{
		qNonL = q_current * w + q_previous * (1 - w);
		reinitialize();
		calculate_global_matrix(qNonL);
		residual = (A.b - A * qNonL).norm() / A.b.norm();
		printf("w = %.3lf -> residual = %.4le\n", w, residual);
		if(w > 0.1) w -= 0.1;
		else w -= 0.01;		
	}
	while(residual > residual_previous && w > 0.01);

	residual_previous = residual;
	return w;
}

void SLAE::simple_iterations()
{
	//FILE *in_f;
	//in_f = fopen("tmp.txt", "w");
	int k_it = 0;
	double w, residual_previous;
	MyVector q0(n), q1(n);
	double normL2u, normL2p;
	MyVector sol_0;
	sol_0.initialize(n);

	create_portret();

	printf("\n%d iteration is in process\n", k_it);

	calculate_global_matrix(sol_0);
	//for(int i = 0; i < A.ggl.size(); i++)
		//fprintf(in_f, "%.20lf\n", A.ggl[i]);
	
	Solve(sol_0, normL2u, normL2p);
	q0 = q_prev;

	get_vector_solution_in_nodes_ux(q0, Ux_numerical);
	get_vector_solution_in_nodes_uy(q0, Uy_numerical);
	get_vector_solution_in_nodes_p(q0, P_numerical);

	normL2u = diff_normL2_u(q0);
	normL2p = diff_normL2_p(q0);

	printf("%le\n", normL2u);
	printf("%le\n", normL2p);

	si_print(logger.log_f, k_it, normL2u, normL2p);

	do
	{
		k_it++;

		printf("\n%d iteration is in process\n", k_it);

		reinitialize();
		calculate_global_matrix(q0);
		Solve(q0, normL2u, normL2p);
		q1 = q_prev;

		if(k_it != 1 && k_it != 2)
		{
			w = find_relaxation_parameter(q1, q0, residual_previous);
			q1 = q1 * w + q0 * (1 - w);
		}
		else
		{
			reinitialize();
			calculate_global_matrix(q1);
			residual_previous = (A.b - A * q1).norm() / A.b.norm();
			///q1.output(in_f);
			printf("\n%.20lf residual\n", residual_previous);
		}

		q0 = q1;
		q_prev = q1;

		get_vector_solution_in_nodes_ux(q0, Ux_numerical);
		get_vector_solution_in_nodes_uy(q0, Uy_numerical);
		get_vector_solution_in_nodes_p(q0, P_numerical);

		normL2u = diff_normL2_u(q0);
		normL2p = diff_normL2_p(q0);

		printf("%le\n", normL2u);
		printf("%le\n", normL2p);

		si_print(logger.log_f, k_it, normL2u, normL2p);
	} while(k_it <= max_iter_nonlinear && residual_previous > 1e-6);

}
#pragma endregion

double SLAE::diff_normL2_p(MyVector q_solution)
{
	double diff_local, diff;
	int size = elements.size();
	double p, function;
	double x0, y0, hx, hy;
	double jacobian;

	diff = 0;
	for(int i = 0; i < size; i++)
	{
		x0 = nodes[elements[i].nodes[0]].x;
		y0 = nodes[elements[i].nodes[0]].y;
		hx = get_hx(i);
		hy = get_hy(i);
		jacobian = hx * hy / 4.0;

		diff_local = 0;
		for(int k = 0; k < 9; k++)
		{
			double p_x = hx * gauss_points[0][k] + x0;
			double p_y = hy * gauss_points[1][k] + y0;
			p = get_solution_in_point_p(p_x, p_y, i, q_solution);
			function = calculate_p_analytic(elements[i].number_of_area, p_x, p_y);
			function -= p;
			diff_local += gauss_weights[k] * function * function;
		}
		diff_local *= jacobian;

		diff += diff_local;
	}

	return sqrt(diff / nodes.size());
}

double SLAE::diff_normL2_u(MyVector q_solution)
{
	double diff_local, diff;
	int size = elements.size();
	double function;
	double x0, y0, hx, hy;
	double jacobian;
	double ux, uy, uxdx, uydy;
	double ux_an, uy_an, uxdx_an, uydy_an;

	diff = 0;
	for(int i = 0; i < size; i++)
	{
		x0 = nodes[elements[i].nodes[0]].x;
		y0 = nodes[elements[i].nodes[0]].y;
		hx = get_hx(i);
		hy = get_hy(i);
		jacobian = hx * hy / 4.0;

		diff_local = 0;
		for(int k = 0; k < 9; k++)
		{
			double p_x = hx * gauss_points[0][k] + x0;
			double p_y = hy * gauss_points[1][k] + y0;

			ux = get_solution_in_point_ux(p_x, p_y, i, q_solution);
			uy = get_solution_in_point_uy(p_x, p_y, i, q_solution);
			uxdx = get_solution_in_point_uxdx(p_x, p_y, i, q_solution);
			uydy = get_solution_in_point_uydy(p_x, p_y, i, q_solution);
			ux_an = calculate_ux_analytic(elements[i].number_of_area, p_x, p_y);
			uy_an = calculate_uy_analytic(elements[i].number_of_area, p_x, p_y);
			uxdx_an = calculate_uxdx_analytic(elements[i].number_of_area, p_x, p_y);
			uydy_an = calculate_uydy_analytic(elements[i].number_of_area, p_x, p_y);

			function = (ux - ux_an) * (ux - ux_an) + (uy - uy_an) * (uy - uy_an);
			function += (uxdx - uxdx_an + uydy - uydy_an) * 
						(uxdx - uxdx_an + uydy - uydy_an);

			diff_local += gauss_weights[k] * function;
		}
		diff_local *= jacobian;

		diff += diff_local;
	}

	return sqrt(diff / nodes.size());
}

void SLAE::run(ofstream& solution_f_out, ofstream& info_f_out)
{
	bool one_research = true;
	double normL2u, normL2p;
	MyVector sol_0;
	sol_0.initialize(n);

	create_portret();

	calculate_global_matrix(sol_0);
	Solve(sol_0, normL2u, normL2p);
	output(solution_f_out, info_f_out, normL2u, normL2p);
}
}
void main()
{
	ifstream l1_in("l1.txt"), grid_in("grid.txt"), elements_in("elements.txt");
	ofstream solution_out("solution.txt"), info_out("info.txt");

	SLAE my_SLAE = SLAE(1000, 100, 1e-12, 30, grid_in, elements_in, l1_in);
	//my_SLAE.run(solution_f_out, info_f_out);
	my_SLAE.simple_iterations();

	l1_in.close();
	grid_in.close();
	elements_in.close();
	solution_out.close();
	info_out.close();

	_getch();
}

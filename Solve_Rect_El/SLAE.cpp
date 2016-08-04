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

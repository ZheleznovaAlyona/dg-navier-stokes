#include "solver.h"

double pow_i(int i, double a)
{
	double res = 1;
	for(int j = 1; j <= i; j++)
		res *= a;
	return res;
}

void Partition::input(FILE *grid_f_in, FILE *elements_f_in)
{
	int tmp;
	Point point_tmp;
	Element element_tmp;

	//считываем узлы
	fscanf(grid_f_in, "%d", &tmp);
	nodes.reserve(tmp);
	for(int i = 0; i < tmp; i++)
	{
		fscanf(grid_f_in, "%lf", &point_tmp.x);
		fscanf(grid_f_in, "%lf", &point_tmp.y);
		nodes.push_back(point_tmp);
	}

	//считываем элементы
	tmp = tmp / 4; //количество элементов
	elements.reserve(tmp);
	for(int i = 0; i < tmp; i++)
	{
		fscanf(elements_f_in, "%d", &element_tmp.number_of_area);
		for(int j = 0; j < 4; j++)
			fscanf(elements_f_in, "%d", &element_tmp.nodes[j]);
		for(int j = 0; j < 4; j++)
			fscanf(elements_f_in, "%d", &element_tmp.neighbors[j]);
		for(int j = 0; j < 4; j++)
			fscanf(elements_f_in, "%d", &element_tmp.edges[j]);

		elements.push_back(element_tmp);
	}
}

void initialize_vector(vector <double> &v, int size)
{
	v.reserve(size);
	for(int i = 0; i < size; i++)
		v.push_back(0.0);
}

MyVector Matrix::Uv(MyVector v)
{
	int i, j, k, kol;
	int iend;
	MyVector new_vector = MyVector(v.ar.size());

	assert(v.ar.size() == n);
	for(i = 0; i < n; i++)
	{
		kol = ig[i+1] - ig[i];//количество ненулевых элементов столбца от первого
								//ненулевого элемента до диагонального элемента (не включая его)
		iend = ig[i+1];
		k = ig[i]; // адрес первого занятого элемента столбца

		new_vector[i] = v[i];//от главной диагонали (у U на диагонали 1)

		for(; k < iend; k++)//проходим по всем элементам i столбца
		{
			j = jg[k];
			new_vector[j] += ggu[k] * v[i];//от верхнего треугольника
		}
	}

	return new_vector;
}

void Matrix::initialize(int size1, int size2)
{
	n = size1; size = size2;

	ggl.reserve(size); ggu.reserve(size);
	for(int i = 0; i < size; i++)
	{
		ggl.push_back(0.0);
		ggu.push_back(0.0);
	}

	di.reserve(n); b.ar.reserve(n);
	for(int i = 0; i < n; i++)
	{
		di.push_back(0.0);
		b.ar.push_back(0.0);
	}

	ig.reserve(n + 1); jg.reserve(size);
	for(int i = 0; i < n + 1; i++)
		ig.push_back(0.0);
	for(int i = 0; i < size; i++)
		jg.push_back(0.0);
}

void Matrix::reinitialize()
{
	for(int i = 0; i < size; i++)
	{
		ggl[i] = 0.0;
		ggu[i] = 0.0;
	}

	for(int i = 0; i < n; i++)
		di[i] = 0.0;

	b.make_zero();
}

void SLAE::initialize(int max_number_of_iterations, int max_number_of_iterations_non_lin, double epsilon, int gmres_m, FILE *grid_f_in, FILE *elements_f_in, FILE *log_f, FILE *boundary1)
{
	int tmp;

	P.input(grid_f_in,elements_f_in);

	n = P.nodes.size() + P.elements.size() * 4; //узлы и рёбра
	max_iter = max_number_of_iterations;
	max_iter_nonlinear = max_number_of_iterations_non_lin;
	eps = epsilon;
	m = gmres_m;

	sigma = 1; mu1 = 1; mu2 = 1;

	tmp = count_unzero_matrix_elements();

	A.initialize(n, tmp);

	Ux_numerical.initialize(P.nodes.size());
	Uy_numerical.initialize(P.nodes.size());
	P_numerical.initialize(P.nodes.size());
	q_prev.initialize(n);

	LU_ggl.reserve(tmp); LU_ggu.reserve(tmp);
	for(int i = 0; i < tmp; i++)
	{
		LU_ggl.push_back(0.0);
		LU_ggu.push_back(0.0);
	}

	LU_di.reserve(n);
	for(int i = 0; i < n; i++)
		LU_di.push_back(0.0);

	yl.initialize(n); yu.initialize(n);
	logger.log_f = log_f;

	input_boundaries1(boundary1);

	phix[0] = [](double ksi, double etta) { return 0.5 * (1 - ksi); };
	phix[1] = [](double ksi, double etta) { return 0.5 * (1 + ksi); };
	phix[2] = [](double ksi, double etta) { return 0.0; };
	phix[3] = [](double ksi, double etta) { return 0.0; };

	phiy[0] = [](double ksi, double etta) { return 0.0; };
	phiy[1] = [](double ksi, double etta) { return 0.0; };
	phiy[2] = [](double ksi, double etta) { return 0.5 * (1 - etta); };
	phiy[3] = [](double ksi, double etta) { return 0.5 * (1 + etta); };

	dphixksi[0] = [](double ksi, double etta) { return -0.5; };
	dphixksi[1] = [](double ksi, double etta) { return 0.5; };
	dphixksi[2] = [](double ksi, double etta) { return 0.0; };
	dphixksi[3] = [](double ksi, double etta) { return  0.0; };

	dphiyksi[0] = [](double ksi, double etta) { return 0.0; };
	dphiyksi[1] = [](double ksi, double etta) { return 0.0; };
	dphiyksi[2] = [](double ksi, double etta) { return 0.0; };
	dphiyksi[3] = [](double ksi, double etta) { return 0.0; };

	dphixetta[0] = [](double ksi, double etta) { return 0.0; };
	dphixetta[1] = [](double ksi, double etta) { return 0.0; };
	dphixetta[2] = [](double ksi, double etta) { return 0.0; };
	dphixetta[3] = [](double ksi, double etta) { return 0.0; };

	dphiyetta[0] = [](double ksi, double etta) { return 0.0; };
	dphiyetta[1] = [](double ksi, double etta) { return 0.0; };
	dphiyetta[2] = [](double ksi, double etta) { return -0.5; };
	dphiyetta[3] = [](double ksi, double etta) { return 0.5; };

	psi[0] = [](double ksi, double etta) { return 0.25 * (1 - ksi) * (1 - etta); };
	psi[1] = [](double ksi, double etta) { return 0.25 * (1 + ksi) * (1 - etta); };
	psi[2] = [](double ksi, double etta) { return 0.25 * (1 - ksi) * (1 + etta); };
	psi[3] = [](double ksi, double etta) { return 0.25 * (1 + ksi) * (1 + etta); };

	dpsiksi[0] = [](double ksi, double etta) { return -0.25 * (1 - etta); };
	dpsiksi[1] = [](double ksi, double etta) { return 0.25 * (1 - etta); };
	dpsiksi[2] = [](double ksi, double etta) { return -0.25 * (1 + etta); };
	dpsiksi[3] = [](double ksi, double etta) { return 0.25 * (1 + etta); };

	dpsietta[0] = [](double ksi, double etta) { return -0.25 * (1 - ksi); };
	dpsietta[1] = [](double ksi, double etta) { return -0.25 * (1 + ksi); };
	dpsietta[2] = [](double ksi, double etta) { return 0.25 * (1 - ksi); };
	dpsietta[3] = [](double ksi, double etta) { return 0.25 * (1 + ksi); };


	double tmp_gauss_points[2][9] =
    {
        {0.0, 0.0, 0.0, sqrt(3.0 / 5.0), -sqrt(3.0 / 5.0), sqrt(3.0 / 5.0), sqrt(3.0 / 5.0), -sqrt(3.0 / 5.0), -sqrt(3.0 / 5.0)},
        {0.0, sqrt(3.0 / 5.0), -sqrt(3.0 / 5.0), 0.0, 0.0, sqrt(3.0 / 5.0), -sqrt(3.0 / 5.0), sqrt(3.0 / 5.0), -sqrt(3.0 / 5.0)}
    };
	for(int i = 0; i < 2; i++)
		for(int j = 0; j < 9; j++)
			gauss_points[i][j] = tmp_gauss_points[i][j];

    gauss_weights[0] = 64.0 / 81.0;
    gauss_weights[1] = gauss_weights[2] = gauss_weights[3] = gauss_weights[4] = 40.0 / 81.0;
    gauss_weights[5] = gauss_weights[6] = gauss_weights[7] = gauss_weights[8] = 25.0 / 81.0;

	gauss_points_1[0] = -sqrt(3.0 / 5.0); gauss_points_1[1] = 0.0; gauss_points_1[2] = sqrt(3.0 / 5.0);
	gauss_weights_1[0] = 5.0 / 9.0; gauss_weights_1[1] = 8.0 / 9.0; gauss_weights_1[2] =  5.0 / 9.0;
}

void SLAE::reinitialize()
{
	Ux_numerical.make_zero();
	Uy_numerical.make_zero();
	P_numerical.make_zero();
	A.reinitialize();

	int gg_size = LU_ggl.size();

	for(int i = 0; i < gg_size; i++)
	{
		LU_ggl[i] = 0.0;
		LU_ggu[i] = 0.0;
	}

	for(int i = 0; i < n; i++)
		LU_di[i] = 0.0;

	if(solver == 0) 
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

double SLAE::get_hx(int element_number)
{
	Element element = P.elements[element_number];
	return P.nodes[element.nodes[1]].x - P.nodes[element.nodes[0]].x;
}
double SLAE::get_hy(int element_number)
{
	Element element = P.elements[element_number];
	return P.nodes[element.nodes[2]].y - P.nodes[element.nodes[0]].y;
}

int SLAE::count_unzero_matrix_elements()
{
	int count_uu = 0;
	for(int i = 0; i < P.elements.size(); i++)
	{
		//с собой
		count_uu += 4;
		//с соседями
		for(int j = 0; j < 4; j++)
			if(P.elements[i].neighbors[j] != -1)
				count_uu += 4;
	}
	count_uu *= 4;
	count_uu -= P.elements.size() * 4;
	//так как нужно для одного треугольника ввиду симметричности портрета,
	//то необходимо полученное количество поделить на 2
	count_uu /= 2;

	int count_pp = 0;
	for(int i = 0; i < P.elements.size(); i++)
	{
		//с собой
		count_pp += 4;
		//с соседями
		for(int j = 0; j < 4; j++)
			if(P.elements[i].neighbors[j] != -1)
				count_pp += 4;
	}
	count_pp *= 4;
	count_pp -= P.nodes.size();
	//так как нужно для одного треугольника ввиду симметричности портрета,
	//то необходимо полученное количество поделить на 2
	count_pp /= 2;

	int count_up = 0;
	for(int i = 0; i < P.elements.size(); i++)
	{
		//с собой
		count_up += 4;
		//с соседями
		for(int j = 0; j < 4; j++)
			if(P.elements[i].neighbors[j] != -1)
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
		neighbor = P.elements[element_number].neighbors[j];
		if(neighbor != -1)
		for(int i = 0; i < dof_num_j; i++)
		{
			if(dof_j_edge) list.push_back(P.elements[neighbor].edges[i]);
			else list.push_back(P.elements[neighbor].nodes[i]);
		}			
	}

	return list.size();
}

#pragma region локальные матрицы и векторы

void SLAE::calculate_locals(int element_number, MyVector q_calc)
{
	Element element = P.elements[element_number];
	int id_i, id_j;
	int n_edges = P.elements.size() * 4;
	
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
			add_element_to_global_matrix(id_i, id_j, G[i][j] + C[i][j]); 
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
	Element element = P.elements[element_number];
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
	Element element = P.elements[element_number];

	hx = get_hx(element_number);
	hy = get_hy(element_number);

	double jacobian = hx * hy / 4.0;
	double u_x, u_y;
	double x0 = P.nodes[element.nodes[0]].x;
	double y0 = P.nodes[element.nodes[0]].y;

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
	Element element = P.elements[element_number];
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
	Element element = P.elements[element_number];
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
	Element element = P.elements[element_number];
	int number_of_area = element.number_of_area;

	double x0 = P.nodes[element.nodes[0]].x;
	double y0 = P.nodes[element.nodes[0]].y;

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

#pragma region границы

#pragma region внутренние границы

void SLAE::calculate_internal_boundaries(int element_number)
{
	Element element = P.elements[element_number];
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
				calculate_ES_vertical(element_number, neighbor_element);
				calculate_P_1_vertical(element_number, neighbor_element);
				calculate_P_2_vertical(element_number, neighbor_element);
				calculate_SP_vertical(element_number, neighbor_element);
			}
			else
			{
				calculate_ES_horizontal(element_number, neighbor_element);
				calculate_P_1_horizontal(element_number, neighbor_element);
				calculate_P_2_horizontal(element_number, neighbor_element);		
				calculate_SP_horizontal(element_number, neighbor_element);	
			}
			add_ES_to_global(element_number, neighbor_element);
			add_P_1_to_global(element_number, neighbor_element);
			add_P_2_to_global(element_number, neighbor_element);
			add_SP_to_global(element_number, neighbor_element);
		}
	}
}

void SLAE::calculate_ES_horizontal(int element_number1, int element_number2)
{
	double AK[4][4], AN[4][4], BK[4][4], BN[4][4];
	double SNN[4][4], SNK[4][4], SKN[4][4], SKK[4][4];
	Element element = P.elements[element_number1];
	Element element_2 = P.elements[element_number2];
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

	for(int i = 0; i < 4; i++)
	{
		for(int j = 0; j < 4; j++)
		{
			AN[i][j] = 0;
			AK[i][j] = 0;
			BN[i][j] = 0;
			BK[i][j] = 0;
			SNN[i][j] = 0;
			SNK[i][j] = 0;
			SKK[i][j] = 0;
			SKN[i][j] = 0;
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

	for(int i = 0; i < 4; i++)
		for(int j = 0; j < 4; j++)
		{
			E[i][j] = AN[i][j] + SNN[i][j];
			E[i + 4][j + 4] = AK[i][j] + SKK[i][j];
			E[i][j + 4] = BN[i][j] + SNK[i][j];
			E[i + 4][j] = BK[i][j] + SKN[i][j];
		}
}
void SLAE::calculate_ES_vertical(int element_number1, int element_number2)
{
	double AK[4][4], AN[4][4], BK[4][4], BN[4][4];
	double SNN[4][4], SNK[4][4], SKN[4][4], SKK[4][4];
	Element element = P.elements[element_number1];
	Element element_2 = P.elements[element_number2];
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
			AN[i][j] = 0;
			AK[i][j] = 0;
			BN[i][j] = 0;
			BK[i][j] = 0;
			SNN[i][j] = 0;
			SNK[i][j] = 0;
			SKK[i][j] = 0;
			SKN[i][j] = 0;
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

	for(int i = 0; i < 4; i++)
		for(int j = 0; j < 4; j++)
		{
			E[i][j] = AN[i][j] + SNN[i][j];
			E[i + 4][j + 4] = AK[i][j] + SKK[i][j];
			E[i][j + 4] = BN[i][j] + SNK[i][j];
			E[i + 4][j] = BK[i][j] + SKN[i][j];
		}
}
void SLAE::add_ES_to_global(int element_number, int neighbor_element_number)
{
	int id_i, id_j;
	Element element = P.elements[element_number];
	Element neighbor_element = P.elements[neighbor_element_number];

	for(int i = 0; i < 4; i++)
	{
		id_i = element.edges[i];
		for(int j = 0; j < 4; j++)
		{				
			id_j = element.edges[j];
			add_element_to_global_matrix(id_i, id_j, E[i][j]); 
		}

		for(int j = 4; j < 8; j++)
		{				
			id_j = neighbor_element.edges[j - 4];
			add_element_to_global_matrix(id_i, id_j, E[i][j]); 
		}
	}

	for(int i = 4; i < 8; i++)
	{
		id_i = neighbor_element.edges[i - 4];
		for(int j = 0; j < 4; j++)
		{				
			id_j = element.edges[j];
			add_element_to_global_matrix(id_i, id_j, E[i][j]); 
		}

		for(int j = 4; j < 8; j++)
		{				
			id_j = neighbor_element.edges[j - 4];
			add_element_to_global_matrix(id_i, id_j, E[i][j]); 
		}
	}
}

void SLAE::calculate_P_1_horizontal(int element_number1, int element_number2)
{
	double AK[4][4], AN[4][4], BK[4][4], BN[4][4];
	Element element = P.elements[element_number1];
	Element element_2 = P.elements[element_number2];
	double rho = calculate_rho(element.number_of_area);
	double rho_2 = calculate_rho(element_2.number_of_area);
	double hx = get_hx(element_number1);

	double a1 =  0.25 * hx * (1 / rho); //якобиан*0.5*(1/rho)
	double a2 = 0.25 * hx * (1 / rho_2);

	double n_vec[2] = {0, 1};

	for(int i = 0; i < 4; i++)
	{
		for(int j = 0; j < 4; j++)
		{
			AN[i][j] = 0;
			AK[i][j] = 0;
			BN[i][j] = 0;
			BK[i][j] = 0;
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

	for(int i = 0; i < 4; i++)
		for(int j = 0; j < 4; j++)
		{
			P_1[i][j] = AN[i][j];
			P_1[i + 4][j + 4] = AK[i][j];
			P_1[i][j + 4] = BN[i][j];
			P_1[i + 4][j] = BK[i][j];
		}
}
void SLAE::calculate_P_1_vertical(int element_number1, int element_number2)
{
	double AK[4][4], AN[4][4], BK[4][4], BN[4][4];
	Element element = P.elements[element_number1];
	Element element_2 = P.elements[element_number2];
	double rho = calculate_rho(element.number_of_area);
	double rho_2 = calculate_rho(element_2.number_of_area);
	double hy = get_hy(element_number1);

	double a1 =  0.25 * hy * (1 / rho); //якобиан*0.5*lambda
	double a2 = 0.25 * hy * (1 / rho_2);

	double n_vec[2] = {1, 0};

	for(int i = 0; i < 4; i++)
	{
		for(int j = 0; j < 4; j++)
		{
			AN[i][j] = 0;
			AK[i][j] = 0;
			BN[i][j] = 0;
			BK[i][j] = 0;
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

	for(int i = 0; i < 4; i++)
		for(int j = 0; j < 4; j++)
		{
			P_1[i][j] = AN[i][j];
			P_1[i + 4][j + 4] = AK[i][j];
			P_1[i][j + 4] = BN[i][j];
			P_1[i + 4][j] = BK[i][j];
		}
}
void SLAE::add_P_1_to_global(int element_number, int neighbor_element_number)
{
	int id_i, id_j;
	Element element = P.elements[element_number];
	Element neighbor_element = P.elements[neighbor_element_number];
	int n_edges = P.elements.size() * 4;

	for(int i = 0; i < 4; i++)
	{
		id_i = element.edges[i];
		for(int j = 0; j < 4; j++)
		{				
			id_j = element.nodes[j] + n_edges;
			add_element_to_global_matrix(id_i, id_j, P_1[i][j]); 
		}

		for(int j = 4; j < 8; j++)
		{				
			id_j = neighbor_element.nodes[j - 4] + n_edges;
			add_element_to_global_matrix(id_i, id_j, P_1[i][j]); 
		}
	}

	for(int i = 4; i < 8; i++)
	{
		id_i = neighbor_element.edges[i - 4];
		for(int j = 0; j < 4; j++)
		{				
			id_j = element.nodes[j] + n_edges;
			add_element_to_global_matrix(id_i, id_j, P_1[i][j]); 
		}

		for(int j = 4; j < 8; j++)
		{				
			id_j = neighbor_element.nodes[j - 4] + n_edges;
			add_element_to_global_matrix(id_i, id_j, P_1[i][j]); 
		}
	}
}

void SLAE::calculate_P_2_horizontal(int element_number1, int element_number2)
{
	double AK[4][4], AN[4][4], BK[4][4], BN[4][4];
	Element element = P.elements[element_number1];
	Element element_2 = P.elements[element_number2];
	double hx = get_hx(element_number1);

	double a1 =  0.25 * hx; //якобиан*0.5
	double a2 = 0.25 * hx;

	double n_vec[2] = {0, 1};

	for(int i = 0; i < 4; i++)
	{
		for(int j = 0; j < 4; j++)
		{
			AN[i][j] = 0;
			AK[i][j] = 0;
			BN[i][j] = 0;
			BK[i][j] = 0;
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

	for(int i = 0; i < 4; i++)
		for(int j = 0; j < 4; j++)
		{
			P_2[i][j] = AN[i][j];
			P_2[i + 4][j + 4] = AK[i][j];
			P_2[i][j + 4] = BN[i][j];
			P_2[i + 4][j] = BK[i][j];
		}
}
void SLAE::calculate_P_2_vertical(int element_number1, int element_number2)
{
	double AK[4][4], AN[4][4], BK[4][4], BN[4][4];
	Element element = P.elements[element_number1];
	Element element_2 = P.elements[element_number2];
	double rho = calculate_rho(element.number_of_area);
	double rho_2 = calculate_rho(element_2.number_of_area);
	double hy = get_hy(element_number1);

	double a1 =  0.25 * hy * (1 / rho); //якобиан*0.5*lambda
	double a2 = 0.25 * hy * (1 / rho_2);

	double n_vec[2] = {1, 0};

	for(int i = 0; i < 4; i++)
	{
		for(int j = 0; j < 4; j++)
		{
			AN[i][j] = 0;
			AK[i][j] = 0;
			BN[i][j] = 0;
			BK[i][j] = 0;
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

	for(int i = 0; i < 4; i++)
		for(int j = 0; j < 4; j++)
		{
			P_2[i][j] = AN[i][j];
			P_2[i + 4][j + 4] = AK[i][j];
			P_2[i][j + 4] = BN[i][j];
			P_2[i + 4][j] = BK[i][j];
		}
}
void SLAE::add_P_2_to_global(int element_number, int neighbor_element_number)
{
	int id_i, id_j;
	Element element = P.elements[element_number];
	Element neighbor_element = P.elements[neighbor_element_number];
	int n_edges = P.elements.size() * 4;

	for(int i = 0; i < 4; i++)
	{
		id_i = element.nodes[i] + n_edges;
		for(int j = 0; j < 4; j++)
		{				
			id_j = element.edges[j];
			add_element_to_global_matrix(id_i, id_j, P_2[i][j]); 
		}

		for(int j = 4; j < 8; j++)
		{				
			id_j = neighbor_element.edges[j - 4];
			add_element_to_global_matrix(id_i, id_j, P_2[i][j]); 
		}
	}

	for(int i = 4; i < 8; i++)
	{
		id_i = neighbor_element.nodes[i - 4] + n_edges;
		for(int j = 0; j < 4; j++)
		{				
			id_j = element.edges[j];
			add_element_to_global_matrix(id_i, id_j, P_2[i][j]); 
		}

		for(int j = 4; j < 8; j++)
		{				
			id_j = neighbor_element.edges[j - 4];
			add_element_to_global_matrix(id_i, id_j, P_2[i][j]); 
		}
	}
}

void SLAE::calculate_SP_horizontal(int element_number1, int element_number2)
{
	double SNN[4][4], SNK[4][4], SKN[4][4], SKK[4][4];
	double hx = get_hx(element_number1);

	double n_vec[2] = {0, 1};

	double jacobian = 0.5 * hx;
	double st = jacobian * mu2;

	for(int i = 0; i < 4; i++)
	{
		for(int j = 0; j < 4; j++)
		{
			SNN[i][j] = 0;
			SNK[i][j] = 0;
			SKK[i][j] = 0;
			SKN[i][j] = 0;
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

	for(int i = 0; i < 4; i++)
		for(int j = 0; j < 4; j++)
		{
			SP[i][j] = SNN[i][j];
			SP[i + 4][j + 4] = SKK[i][j];
			SP[i][j + 4] = SNK[i][j];
			SP[i + 4][j] = SKN[i][j];
		}
}
void SLAE::calculate_SP_vertical(int element_number1, int element_number2)
{
	double SNN[4][4], SNK[4][4], SKN[4][4], SKK[4][4];
	double hy = get_hy(element_number1);

	double n_vec[2] = {1, 0};

	double jacobian = 0.5 * hy;
	double st = jacobian * mu2;

	for(int i = 0; i < 4; i++)
	{
		for(int j = 0; j < 4; j++)
		{
			SNN[i][j] = 0;
			SNK[i][j] = 0;
			SKK[i][j] = 0;
			SKN[i][j] = 0;
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

	for(int i = 0; i < 4; i++)
		for(int j = 0; j < 4; j++)
		{
			SP[i][j] = SNN[i][j];
			SP[i + 4][j + 4] = SKK[i][j];
			SP[i][j + 4] = SNK[i][j];
			SP[i + 4][j] = SKN[i][j];
		}
}
void SLAE::add_SP_to_global(int element_number, int neighbor_element_number)
{
	int id_i, id_j;
	Element element = P.elements[element_number];
	Element neighbor_element = P.elements[neighbor_element_number];
	int n_edges = P.elements.size() * 4;

	for(int i = 0; i < 4; i++)
	{
		id_i = element.nodes[i] + n_edges;
		for(int j = 0; j < 4; j++)
		{				
			id_j = element.nodes[j] + n_edges;
			add_element_to_global_matrix(id_i, id_j, SP[i][j]); 
		}

		for(int j = 4; j < 8; j++)
		{				
			id_j = neighbor_element.nodes[j - 4] + n_edges;
			add_element_to_global_matrix(id_i, id_j, SP[i][j]); 
		}
	}

	for(int i = 4; i < 8; i++)
	{
		id_i = neighbor_element.nodes[i - 4] + n_edges;
		for(int j = 0; j < 4; j++)
		{				
			id_j = element.nodes[j] + n_edges;
			add_element_to_global_matrix(id_i, id_j, SP[i][j]); 
		}

		for(int j = 4; j < 8; j++)
		{				
			id_j = neighbor_element.nodes[j - 4] + n_edges;
			add_element_to_global_matrix(id_i, id_j, SP[i][j]); 
		}
	}
}


#pragma endregion

#pragma region внешние границы

void SLAE::calculate_outer_boundaries(int element_number)
{
	Element element = P.elements[element_number];

	int last_node = P.nodes.size() - 1;
	int left_low_corner_node = element.nodes[0];
	int right_up_corner_node = element.nodes[8];

	if(P.nodes[left_low_corner_node].x == P.nodes[0].x)
	{
		calculate_ES_out_left(element_number);
		calculate_P_1_out_left(element_number);
		calculate_P_2_out_left(element_number);
		calculate_SP_out_left(element_number);
	}//вертикальная левая граница
	if(P.nodes[right_up_corner_node].x == P.nodes[last_node].x)
	{
		calculate_ES_out_right(element_number);
		calculate_P_1_out_right(element_number);
		calculate_P_2_out_right(element_number);
		calculate_SP_out_right(element_number);
	}//вертикальная правая граница
	if(P.nodes[left_low_corner_node].y == P.nodes[0].y) 
	{
		calculate_ES_out_low(element_number);
		calculate_P_1_out_low(element_number);
		calculate_P_2_out_low(element_number);
		calculate_SP_out_low(element_number);
	}//горизонтальная нижняя граница
	if(P.nodes[right_up_corner_node].y == P.nodes[last_node].y)	
	{
		calculate_ES_out_up(element_number);
		calculate_P_1_out_up(element_number);
		calculate_P_2_out_up(element_number);
		calculate_SP_out_up(element_number);
	}//горизонтальная верхняя граница
}

void SLAE::calculate_ES_out_left(int element_number)
{
	double S_out[4][4];
	Element element = P.elements[element_number];
	double lambda = calculate_lambda(element.number_of_area);
	double hx = get_hx(element_number);
	double hy = get_hy(element_number);

	double a =  0.5 * hy * lambda; //якобиан*lambda

	double jacobian = 0.5 * hy;
	double st = jacobian * sigma;

	double n_vec[2] = {-1, 0};
	double n_vec2[2] = {1, 0};

	for(int i = 0; i < 4; i++)
	{
		for(int j = 0; j < 4; j++)
		{
			E_out[i][j] = 0;
			S_out[i][j] = 0;
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

	for(int i = 0; i < 4; i++)
	{
		int id_i = element.edges[i];
		for(int j = 0; j < 4; j++)
		{				
			int id_j = element.edges[j];
			add_element_to_global_matrix(id_i, id_j, E_out[i][j] + S_out[i][j]); 
		}
	}
}
void SLAE::calculate_ES_out_right(int element_number)
{
	double S_out[4][4];
	Element element = P.elements[element_number];
	double lambda = calculate_lambda(element.number_of_area);
	double hx = get_hx(element_number);
	double hy = get_hy(element_number);

	double a =  0.5 * hy * lambda; //якобиан*lambda

	double jacobian = 0.5 * hy;
	double st = jacobian * sigma;

	double n_vec[2] = {1, 0};
	double n_vec2[2] = {1, 0};

	for(int i = 0; i < 4; i++)
	{
		for(int j = 0; j < 4; j++)
		{
			E_out[i][j] = 0;
			S_out[i][j] = 0;
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

	for(int i = 0; i < 4; i++)
	{
		int id_i = element.edges[i];
		for(int j = 0; j < 4; j++)
		{				
			int id_j = element.edges[j];
			add_element_to_global_matrix(id_i, id_j, E_out[i][j] + S_out[i][j]); 
		}
	}
}
void SLAE::calculate_ES_out_low(int element_number)
{
	double S_out[4][4];
	Element element = P.elements[element_number];
	double lambda = calculate_lambda(element.number_of_area);
	double hx = get_hx(element_number);
	double hy = get_hy(element_number);

	double a =  0.5 * hx * lambda; //якобиан*lambda

	double jacobian = 0.5 * hx;
	double st = jacobian * sigma;

	double n_vec[2] = {0, -1};
	double n_vec2[2] = {0, 1};

	for(int i = 0; i < 4; i++)
	{
		for(int j = 0; j < 4; j++)
		{
			E_out[i][j] = 0;
			S_out[i][j] = 0;
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

	for(int i = 0; i < 4; i++)
	{
		int id_i = element.edges[i];
		for(int j = 0; j < 4; j++)
		{				
			int id_j = element.edges[j];
			add_element_to_global_matrix(id_i, id_j, E_out[i][j] + S_out[i][j]); 
		}
	}
}
void SLAE::calculate_ES_out_up(int element_number)
{
	double S_out[4][4];
	Element element = P.elements[element_number];
	double lambda = calculate_lambda(element.number_of_area);
	double hx = get_hx(element_number);
	double hy = get_hy(element_number);

	double a =  0.5 * hx * lambda; //якобиан*lambda

	double jacobian = 0.5 * hx;
	double st = jacobian * sigma;

	double n_vec[2] = {0, 1};
	double n_vec2[2] = {0, 1};

	for(int i = 0; i < 4; i++)
	{
		for(int j = 0; j < 4; j++)
		{
			E_out[i][j] = 0;
			S_out[i][j] = 0;
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

	for(int i = 0; i < 4; i++)
	{
		int id_i = element.edges[i];
		for(int j = 0; j < 4; j++)
		{				
			int id_j = element.edges[j];
			add_element_to_global_matrix(id_i, id_j, E_out[i][j] + S_out[i][j]); 
		}
	}
}

void SLAE::calculate_P_1_out_left(int element_number)
{
	Element element = P.elements[element_number];
	double rho = calculate_rho(element.number_of_area);
	double hy = get_hy(element_number);

	double a =  0.5 * hy * (1 / rho); //якобиан*1/rho

	double n_vec[2] = {-1, 0};

	for(int i = 0; i < 4; i++)
	{
		for(int j = 0; j < 4; j++)
		{
			P_1_out[i][j] = 0;
			for(int k = 0; k < 3; k++)
			{
				double p_etta = gauss_points_1[k];
				P_1_out[i][j] += gauss_weights_1[k]  * psi[j](-1, p_etta) *
								 (phix[i](-1, p_etta) * n_vec[0] + phiy[i](-1, p_etta) * n_vec[1]);

			} 
			P_1_out[i][j] *= a;
		}
	}

	int n_edges = P.elements.size() * 4;

	for(int i = 0; i < 4; i++)
	{
		int id_i = element.edges[i];
		for(int j = 0; j < 4; j++)
		{				
			int id_j = element.nodes[j] + n_edges;
			add_element_to_global_matrix(id_i, id_j, P_1_out[i][j]); 
		}
	}
}
void SLAE::calculate_P_1_out_right(int element_number)
{
	Element element = P.elements[element_number];
	double rho = calculate_rho(element.number_of_area);
	double hy = get_hy(element_number);

	double a =  0.5 * hy * (1 / rho); //якобиан*1/rho

	double n_vec[2] = {1, 0};

	for(int i = 0; i < 4; i++)
	{
		for(int j = 0; j < 4; j++)
		{
			P_1_out[i][j] = 0;
			for(int k = 0; k < 3; k++)
			{
				double p_etta = gauss_points_1[k];
				P_1_out[i][j] += gauss_weights_1[k]  * psi[j](1, p_etta) *
								 (phix[i](1, p_etta) * n_vec[0] + phiy[i](1, p_etta) * n_vec[1]);

			} 
			P_1_out[i][j] *= a;
		}
	}

	int n_edges = P.elements.size() * 4;

	for(int i = 0; i < 4; i++)
	{
		int id_i = element.edges[i];
		for(int j = 0; j < 4; j++)
		{				
			int id_j = element.nodes[j] + n_edges;
			add_element_to_global_matrix(id_i, id_j, P_1_out[i][j]); 
		}
	}
}
void SLAE::calculate_P_1_out_low(int element_number)
{
	Element element = P.elements[element_number];
	double rho = calculate_rho(element.number_of_area);
	double hx = get_hx(element_number);

	double a =  0.5 * hx * (1 / rho); //якобиан*1/rho

	double n_vec[2] = {0, -1};

	for(int i = 0; i < 4; i++)
	{
		for(int j = 0; j < 4; j++)
		{
			P_1_out[i][j] = 0;
			for(int k = 0; k < 3; k++)
			{
				double p_ksi = gauss_points_1[k];
				P_1_out[i][j] += gauss_weights_1[k]  * psi[j](p_ksi, -1) *
								 (phix[i](p_ksi, -1) * n_vec[0] + phiy[i](p_ksi, -1) * n_vec[1]);

			} 
			P_1_out[i][j] *= a;
		}
	}

	int n_edges = P.elements.size() * 4;

	for(int i = 0; i < 4; i++)
	{
		int id_i = element.edges[i];
		for(int j = 0; j < 4; j++)
		{				
			int id_j = element.nodes[j] + n_edges;
			add_element_to_global_matrix(id_i, id_j, P_1_out[i][j]); 
		}
	}
}
void SLAE::calculate_P_1_out_up(int element_number)
{
	Element element = P.elements[element_number];
	double rho = calculate_rho(element.number_of_area);
	double hx = get_hx(element_number);

	double a =  0.5 * hx * (1 / rho); //якобиан*1/rho

	double n_vec[2] = {0, 1};

	for(int i = 0; i < 4; i++)
	{
		for(int j = 0; j < 4; j++)
		{
			P_1_out[i][j] = 0;
			for(int k = 0; k < 3; k++)
			{
				double p_ksi = gauss_points_1[k];
				P_1_out[i][j] += gauss_weights_1[k]  * psi[j](p_ksi, 1) *
								 (phix[i](p_ksi, 1) * n_vec[0] + phiy[i](p_ksi, 1) * n_vec[1]);

			} 
			P_1_out[i][j] *= a;
		}
	}

	int n_edges = P.elements.size() * 4;

	for(int i = 0; i < 4; i++)
	{
		int id_i = element.edges[i];
		for(int j = 0; j < 4; j++)
		{				
			int id_j = element.nodes[j] + n_edges;
			add_element_to_global_matrix(id_i, id_j, P_1_out[i][j]); 
		}
	}
}

void SLAE::calculate_P_2_out_left(int element_number)
{
	Element element = P.elements[element_number];
	double hy = get_hy(element_number);

	double a =  -0.5 * hy; //-якобиан

	double n_vec[2] = {-1, 0};

	for(int i = 0; i < 4; i++)
	{
		for(int j = 0; j < 4; j++)
		{
			P_2_out[i][j] = 0;
			for(int k = 0; k < 3; k++)
			{
				double p_etta = gauss_points_1[k];
				P_2_out[i][j] += gauss_weights_1[k]  * psi[i](-1, p_etta) *
								 (phix[j](-1, p_etta) * n_vec[0] + phiy[j](-1, p_etta) * n_vec[1]);

			} 
			P_2_out[i][j] *= a;
		}
	}

	int n_edges = P.elements.size() * 4;

	for(int i = 0; i < 4; i++)
	{
		int id_i = element.nodes[i] + n_edges;
		for(int j = 0; j < 4; j++)
		{				
			int id_j = element.edges[j];
			add_element_to_global_matrix(id_i, id_j, P_2_out[i][j]); 
		}
	}
}
void SLAE::calculate_P_2_out_right(int element_number)
{
	Element element = P.elements[element_number];
	double hy = get_hy(element_number);

	double a =  -0.5 * hy; //-якобиан

	double n_vec[2] = {1, 0};

	for(int i = 0; i < 4; i++)
	{
		for(int j = 0; j < 4; j++)
		{
			P_2_out[i][j] = 0;
			for(int k = 0; k < 3; k++)
			{
				double p_etta = gauss_points_1[k];
				P_2_out[i][j] += gauss_weights_1[k]  * psi[i](1, p_etta) *
								 (phix[j](1, p_etta) * n_vec[0] + phiy[j](1, p_etta) * n_vec[1]);

			} 
			P_2_out[i][j] *= a;
		}
	}

	int n_edges = P.elements.size() * 4;

	for(int i = 0; i < 4; i++)
	{
		int id_i = element.nodes[i] + n_edges;
		for(int j = 0; j < 4; j++)
		{				
			int id_j = element.edges[j];
			add_element_to_global_matrix(id_i, id_j, P_2_out[i][j]); 
		}
	}
}
void SLAE::calculate_P_2_out_low(int element_number)
{
	Element element = P.elements[element_number];
	double hx = get_hx(element_number);

	double a =  -0.5 * hx; //-якобиан

	double n_vec[2] = {0, -1};

	for(int i = 0; i < 4; i++)
	{
		for(int j = 0; j < 4; j++)
		{
			P_2_out[i][j] = 0;
			for(int k = 0; k < 3; k++)
			{
				double p_ksi = gauss_points_1[k];
				P_2_out[i][j] += gauss_weights_1[k]  * psi[i](p_ksi, -1) *
								 (phix[j](p_ksi, -1) * n_vec[0] + phiy[j](p_ksi, -1) * n_vec[1]);

			} 
			P_2_out[i][j] *= a;
		}
	}

	int n_edges = P.elements.size() * 4;

	for(int i = 0; i < 4; i++)
	{
		int id_i = element.nodes[i] + n_edges;
		for(int j = 0; j < 4; j++)
		{				
			int id_j = element.edges[j];
			add_element_to_global_matrix(id_i, id_j, P_2_out[i][j]); 
		}
	}
}
void SLAE::calculate_P_2_out_up(int element_number)
{
	Element element = P.elements[element_number];
	double hx = get_hx(element_number);

	double a =  -0.5 * hx; //-якобиан

	double n_vec[2] = {0, 1};

	for(int i = 0; i < 4; i++)
	{
		for(int j = 0; j < 4; j++)
		{
			P_2_out[i][j] = 0;
			for(int k = 0; k < 3; k++)
			{
				double p_ksi = gauss_points_1[k];
				P_2_out[i][j] += gauss_weights_1[k]  * psi[i](p_ksi, 1) *
								 (phix[j](p_ksi, 1) * n_vec[0] + phiy[j](p_ksi, 1) * n_vec[1]);

			} 
			P_2_out[i][j] *= a;
		}
	}

	int n_edges = P.elements.size() * 4;

	for(int i = 0; i < 4; i++)
	{
		int id_i = element.nodes[i] + n_edges;
		for(int j = 0; j < 4; j++)
		{				
			int id_j = element.edges[j];
			add_element_to_global_matrix(id_i, id_j, P_2_out[i][j]); 
		}
	}
}

void SLAE::calculate_SP_out_left(int element_number)
{
	Element element = P.elements[element_number];
	double hy = get_hy(element_number);

	double n_vec[2] = {0, -1};

	double jacobian = 0.5 * hy;
	double st = jacobian * mu2;

	for(int i = 0; i < 4; i++)
	{
		for(int j = 0; j < 4; j++)
		{
			SP_out[i][j] = 0;

			for(int k = 0; k < 3; k++)
			{
				double p_etta = gauss_points_1[k];
				SP_out[i][j] += gauss_weights_1[k] * psi[i](-1, p_etta) * psi[j](-1, p_etta);

			} 
			SP_out[i][j] *= st;
		}
	}

	int n_edges = P.elements.size() * 4;

	for(int i = 0; i < 4; i++)
	{
		int id_i = element.nodes[i] + n_edges;
		for(int j = 0; j < 4; j++)
		{				
			int id_j = element.nodes[j] + n_edges;
			add_element_to_global_matrix(id_i, id_j, SP_out[i][j]); 
		}
	}
}
void SLAE::calculate_SP_out_right(int element_number)
{
	Element element = P.elements[element_number];
	double hy = get_hy(element_number);

	double n_vec[2] = {0, 1};

	double jacobian = 0.5 * hy;
	double st = jacobian * mu2;

	for(int i = 0; i < 4; i++)
	{
		for(int j = 0; j < 4; j++)
		{
			SP_out[i][j] = 0;

			for(int k = 0; k < 3; k++)
			{
				double p_etta = gauss_points_1[k];
				SP_out[i][j] += gauss_weights_1[k] * psi[i](1, p_etta) * psi[j](1, p_etta);

			} 
			SP_out[i][j] *= st;
		}
	}

	int n_edges = P.elements.size() * 4;

	for(int i = 0; i < 4; i++)
	{
		int id_i = element.nodes[i] + n_edges;
		for(int j = 0; j < 4; j++)
		{				
			int id_j = element.nodes[j] + n_edges;
			add_element_to_global_matrix(id_i, id_j, SP_out[i][j]); 
		}
	}
}
void SLAE::calculate_SP_out_low(int element_number)
{
	Element element = P.elements[element_number];
	double hx = get_hx(element_number);

	double n_vec[2] = {-1, 0};

	double jacobian = 0.5 * hx;
	double st = jacobian * mu2;

	for(int i = 0; i < 4; i++)
	{
		for(int j = 0; j < 4; j++)
		{
			SP_out[i][j] = 0;

			for(int k = 0; k < 3; k++)
			{
				double p_ksi = gauss_points_1[k];
				SP_out[i][j] += gauss_weights_1[k] * psi[i](p_ksi, -1) * psi[j](p_ksi, -1);

			} 
			SP_out[i][j] *= st;
		}
	}

	int n_edges = P.elements.size() * 4;

	for(int i = 0; i < 4; i++)
	{
		int id_i = element.nodes[i] + n_edges;
		for(int j = 0; j < 4; j++)
		{				
			int id_j = element.nodes[j] + n_edges;
			add_element_to_global_matrix(id_i, id_j, SP_out[i][j]); 
		}
	}
}
void SLAE::calculate_SP_out_up(int element_number)
{
	Element element = P.elements[element_number];
	double hx = get_hx(element_number);

	double n_vec[2] = {1, 0};

	double jacobian = 0.5 * hx;
	double st = jacobian * mu2;

	for(int i = 0; i < 4; i++)
	{
		for(int j = 0; j < 4; j++)
		{
			SP_out[i][j] = 0;

			for(int k = 0; k < 3; k++)
			{
				double p_ksi = gauss_points_1[k];
				SP_out[i][j] += gauss_weights_1[k] * psi[i](p_ksi, 1) * psi[j](p_ksi, 1);

			} 
			SP_out[i][j] *= st;
		}
	}

	int n_edges = P.elements.size() * 4;

	for(int i = 0; i < 4; i++)
	{
		int id_i = element.nodes[i] + n_edges;
		for(int j = 0; j < 4; j++)
		{				
			int id_j = element.nodes[j] + n_edges;
			add_element_to_global_matrix(id_i, id_j, SP_out[i][j]); 
		}
	}
}

#pragma endregion

#pragma region  краевые условия
void SLAE::input_boundaries1(FILE* f_in)
{
	int count;
	BoundaryCondition1 tmp;

	fscanf(f_in, "%d", &count);	
	boundaries1.reserve(count);

	for(int i = 1; i <= count; i++)		
	{
		fscanf(f_in, "%d", &tmp.elem);
		fscanf(f_in, "%d", &tmp.formula_number);
		fscanf(f_in, "%d", &tmp.edges[0]);
		fscanf(f_in, "%d", &tmp.edges[1]);
		fscanf(f_in, "%d", &tmp.edges[2]);
		fscanf(f_in, "%d", &tmp.edges[3]);
		boundaries1.push_back(tmp);
	}
}

void SLAE::input_boundaries2(FILE* f_in)
{
	int count;
	BoundaryCondition2 tmp;

	fscanf(f_in, "%d", &count);	
	boundaries2.reserve(count);

	for(int i = 1; i <= count; i++)		
	{
		fscanf(f_in, "%d", &tmp.elem);
		fscanf(f_in, "%d", &tmp.formula_number);
		fscanf(f_in, "%d", &tmp.edges[0]);
		fscanf(f_in, "%d", &tmp.edges[1]);
		fscanf(f_in, "%d", &tmp.edges[2]);
		fscanf(f_in, "%d", &tmp.edges[3]);
		boundaries2.push_back(tmp);
	}
}

void SLAE::input_boundaries3(FILE* f_in)
{
	int count;
	BoundaryCondition3 tmp;

	fscanf(f_in, "%d", &count);	
	boundaries3.reserve(count);

	for(int i = 1; i <= count; i++)		
	{
		fscanf(f_in, "%d", &tmp.elem);
		fscanf(f_in, "%d", &tmp.formula_number);
		fscanf(f_in, "%d", &tmp.edges[0]);
		fscanf(f_in, "%d", &tmp.edges[1]);
		fscanf(f_in, "%d", &tmp.edges[2]);
		fscanf(f_in, "%d", &tmp.edges[3]);
		boundaries3.push_back(tmp);
	}
}

void SLAE::calculate_all_boundaries1()
{
	int size_b = boundaries1.size();
	for(int i = 0; i < size_b; i++)
		calculate_boundaries1(i);
}

#pragma region слабые первые краевые условия

void SLAE::calculate_boundaries1_left(int number)
{
	Element element = P.elements[boundaries1[number].elem];
	double hy = get_hy(boundaries1[number].elem);
	double hx = get_hx(boundaries1[number].elem);

	double x0 = P.nodes[element.nodes[0]].x;
	double y0 = P.nodes[element.nodes[0]].y;

	double lambda = calculate_lambda(element.number_of_area);
	double n_vec[2] = {-1, 0};

	double Ug_vector[4];
	double Pg_vector[4];
	double g_x, g_y;
	double jacobian = hy / 2;
	int n_edges = P.elements.size() * 4;

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

void SLAE::calculate_boundaries1_right(int number)
{
	Element element = P.elements[boundaries1[number].elem];
	double hy = get_hy(boundaries1[number].elem);
	double hx = get_hx(boundaries1[number].elem);

	double x0 = P.nodes[element.nodes[0]].x;
	double y0 = P.nodes[element.nodes[0]].y;

	double lambda = calculate_lambda(element.number_of_area);
	double n_vec[2] = {1, 0};

	double Ug_vector[4];
	double Pg_vector[4];
	double g_x, g_y;
	double jacobian = hy / 2;
	int n_edges = P.elements.size() * 4;

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

void SLAE::calculate_boundaries1_low(int number)
{
	Element element = P.elements[boundaries1[number].elem];
	double hy = get_hy(boundaries1[number].elem);
	double hx = get_hx(boundaries1[number].elem);

	double x0 = P.nodes[element.nodes[0]].x;
	double y0 = P.nodes[element.nodes[0]].y;

	double lambda = calculate_lambda(element.number_of_area);
	double n_vec[2] = {0, -1};

	double Ug_vector[4];
	double Pg_vector[4];
	double g_x, g_y;
	double jacobian = hx / 2;
	int n_edges = P.elements.size() * 4;

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

void SLAE::calculate_boundaries1_up(int number)
{
	Element element = P.elements[boundaries1[number].elem];
	double hy = get_hy(boundaries1[number].elem);
	double hx = get_hx(boundaries1[number].elem);

	double x0 = P.nodes[element.nodes[0]].x;
	double y0 = P.nodes[element.nodes[0]].y;

	double lambda = calculate_lambda(element.number_of_area);
	double n_vec[2] = {0, 1};

	double Ug_vector[4];
	double Pg_vector[4];
	double g_x, g_y;
	double jacobian = hx / 2;
	int n_edges = P.elements.size() * 4;

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

#pragma endregion

void SLAE::calculate_boundaries1(int number)
{
	if(boundaries1[number].edges[0] == 1) calculate_boundaries1_left(number);
	if(boundaries1[number].edges[1] == 1) calculate_boundaries1_right(number);
	if(boundaries1[number].edges[2] == 1) calculate_boundaries1_low(number);
	if(boundaries1[number].edges[3] == 1) calculate_boundaries1_up(number);
}

double SLAE::gx(int formula_number, double x, double y)
{
	if(test == 1)
		switch(formula_number)
		{
			case 0: return y; break;
			case 1:	return y; break;
		}
	if(test == 3)
		switch(formula_number)
		{
			case 0: return 20 * x * y * y * y; break;
			case 1:	return 20 * x * y * y * y; break;
		}
}

double SLAE::gy(int formula_number, double x, double y)
{
	if(test == 1)
		switch(formula_number)
		{
			case 0: return x; break;
			case 1:	return x; break;
		}
	if(test == 3)
		switch(formula_number)
		{
			case 0: return 5 * x * x * x * x - 5 * y * y * y * y; break;
			case 1:	return 5 * x * x * x * x - 5 * y * y * y * y; break;
		}
}

#pragma endregion

#pragma endregion

#pragma region выражения для функций и параметров
double SLAE::calculate_fx(int area_number, double x, double y)
{
	if(test == 1)
		switch(area_number)
		{
			case 0: return 1.0 / calculate_rho(area_number); break;
			case 1:	return 1.0 / calculate_rho(area_number); break;
		}
	//if(test == 3)
	//	switch(area_number)
	//	{
	//	case 0: return -calculate_lambda(area_number) * 120 * x * y + 1.0 / calculate_rho(area_number) * 120 * x * y; break;
	//	case 1:	return -calculate_lambda(area_number) * 120 * x * y + 1.0 / calculate_rho(area_number) * 120 * x * y; break;
	//	}
	if(test == 3)
		switch(area_number)
		{
		case 0: return -calculate_lambda(area_number) * 120 * x * y + 
						1.0 / calculate_rho(area_number) * 120 * x * y + 
						400 * x * pow_i(6, y) + 300 * (pow_i(5, x) * y * y
						- x * pow_i(6, y)); break;
		case 1:	return -calculate_lambda(area_number) * 120 * x * y + 
						1.0 / calculate_rho(area_number) * 120 * x * y + 
						400 * x * pow_i(6, y) + 300 * (pow_i(5, x) * y * y
						- x * pow_i(6, y)); break;
		}
}

double SLAE::calculate_fy(int area_number, double x, double y)
{
	if(test == 1)
		switch(area_number)
		{
			case 0: return 1.0 / calculate_rho(area_number); break;
			case 1:	return 1.0 / calculate_rho(area_number); break;
		}
	//if(test == 3)
	//	switch(area_number)
	//	{
	//	case 0: return -calculate_lambda(area_number) * (60 * x * x - 60 * y * y) + 1.0 / calculate_rho(area_number) * (60 * x * x - 60 * y * y); break;
	//	case 1:	return -calculate_lambda(area_number) * (60 * x * x - 60 * y * y) + 1.0 / calculate_rho(area_number) * (60 * x * x - 60 * y * y); break;
	//	}
	if(test == 3)
		switch(area_number)
		{
		case 0: return -calculate_lambda(area_number) * (60 * x * x - 60 * y * y) + 
						1.0 / calculate_rho(area_number) * (60 * x * x - 60 * y * y) + 
						400 * pow_i(4, x)  * pow_i(3, y) + 100 * (pow_i(7, y) 
						- pow_i(4, x)  * pow_i(3, y)); break;
		case 1:	return -calculate_lambda(area_number) * (60 * x * x - 60 * y * y) + 
						1.0 / calculate_rho(area_number) * (60 * x * x - 60 * y * y) + 
						400 * pow_i(4, x)  * pow_i(3, y) + 100 * (pow_i(7, y)
						- pow_i(4, x)  * pow_i(3, y)); break;
		}
}

double SLAE::calculate_ux_analytic(int area_number, double x, double y)
{
	if(test == 1)
		switch(area_number)
		{
			case 0: return y; break;
			case 1:	return y; break;
		}
	if(test == 3)
		switch(area_number)
		{
			case 0: return 20 * x * y * y * y; break;
			case 1:	return 20 * x * y * y * y; break;
		}
}

double SLAE::calculate_uy_analytic(int area_number, double x, double y)
{
	if(test == 1)
		switch(area_number)
		{
			case 0: return x; break;
			case 1:	return x; break;
		}
	if(test == 3)
		switch(area_number)
		{
			case 0: return 5 * x * x * x * x - 5 * y * y * y * y; break;
			case 1:	return 5 * x * x * x * x - 5 * y * y * y * y; break;
		}
}

double SLAE::calculate_uxdx_analytic(int area_number, double x, double y)
{
	if(test == 1)
		switch(area_number)
		{
			case 0: return 0; break;
			case 1:	return 0; break;
		}
	if(test == 3)
		switch(area_number)
		{
		case 0: return 20.0 * y * y * y; break;
		case 1:	return 20.0 * y * y * y; break;
		}
}

double SLAE::calculate_uydy_analytic(int area_number, double x, double y)
{
	if(test == 1)
		switch(area_number)
		{
			case 0: return 0; break;
			case 1:	return 0; break;
		}
	if(test == 3)
		switch(area_number)
		{
		case 0: return -20.0 * y * y * y; break;
		case 1:	return -20.0 * y * y * y; break;
		}
}

double SLAE::calculate_p_analytic(int area_number, double x, double y)
{
	if(test == 1)
		switch(area_number)
		{
			case 0: return x + y - 1; break;
			case 1:	return x + y - 1; break;
		}
	if(test == 3)
		switch(area_number)
		{
			case 0: return 60 * x * x * y - 20 * y * y * y - 5; break;
			case 1:	return 60 * x * x * y - 20 * y * y * y - 5; break;
		}
}

double SLAE::calculate_lambda(int area_number)
{
	if(test == 1)
		switch(area_number)
		{
			case 0: return 1.0; break;
			case 1:	return 1.0; break;
		}
	if(test == 3)
		switch(area_number)
		{
			case 0: return 1.0; break;
			case 1:	return 1.0; break;
		}
}

double SLAE::calculate_rho(int area_number)
{
	if(test == 1)
		switch(area_number)
		{
			case 0: return 1.0; break;
			case 1:	return 1.0; break;
		}
	if(test == 3)
		switch(area_number)
		{
			case 0: return 1.0; break;
			case 1:	return 1.0; break;
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
	int n_edges = P.elements.size() * 4;

	lists = new vector <int>[n];
	int count_elements = P.elements.size();

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
																 P.elements[i].edges, 
																 P.elements[i].edges, 
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
																 P.elements[i].nodes, 
																 P.elements[i].nodes, 
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
																 P.elements[i].nodes, 
																 P.elements[i].edges, 
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
			for(int j = 0; j < lists[i].size(); j++)
			{
				A.jg[k] = lists[i][j];
				k++;
			}
			lists[i].clear();
		}
	}

	delete[] lists;
}
void SLAE::add_element_to_global_matrix(int i, int j, double element)
{
	int id;
	bool flag;

	if(i == j)
		A.di[i] += element;
	else
	{
		if(i < j)
		{	
			flag = false;
			for(id = A.ig[j]; !flag && id <= A.ig[j + 1] - 1; id++)
				if(A.jg[id] == i) flag = true;
			 if(flag) A.ggu[id - 1] += element;
		}
		else
		{
			flag = false;
			for(id = A.ig[i]; !flag && id <= A.ig[i + 1] - 1; id++)
				if(A.jg[id] == j) flag = true;
			if(flag) A.ggl[id - 1] += element;
		}
	}
}

void SLAE::put_element_to_global_matrix(int i, int j, double element)
{
	int id;
	bool flag;

	if(i == j)
		A.di[i] = element;
	else
	{
		if(i < j)
		{	
			flag = false;
			for(id = A.ig[j]; !flag && id <= A.ig[j + 1] - 1; id++)
				if(A.jg[id] == i) flag = true;
			 if(flag) A.ggu[id - 1] = element;
		}
		else
		{
			flag = false;
			for(id = A.ig[i]; !flag && id <= A.ig[i + 1] - 1; id++)
				if(A.jg[id] == j) flag = true;
			if(flag) A.ggl[id - 1] = element;
		}
	}
}
void SLAE::calculate_global_matrix(MyVector q_calc)
{
	int size = P.elements.size();
	int id_i, id_j;

	//локаьные матрицы и вектор правой части
	for(int el_i = 0; el_i < size; el_i++)
		calculate_locals(el_i, q_calc);

	//матрицы межэлементных границ
	for(int el_i = 0; el_i < size; el_i++)
		calculate_internal_boundaries(el_i);

	//расчёт матриц границ области
	//if(out_boundaries)
		for(int el_i = 0; el_i < size; el_i++)
			calculate_outer_boundaries(el_i);

	//учёт первых краевых условий
	calculate_all_boundaries1();
}

#pragma endregion

double SLAE::get_solution_in_point_ux(double x, double y, int element_number, MyVector qi)
{
	int indexes[4];
	double u_in_point, qi_local[4];

	//собираем глобальные номера с элемента
	for(int j = 0; j < 4; j++)
		indexes[j] = P.elements[element_number].edges[j];

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
		indexes[j] = P.elements[element_number].edges[j];

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
		indexes[j] = P.elements[element_number].edges[j];

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
		indexes[j] = P.elements[element_number].edges[j];

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
	int indexes[4], n_edges = P.elements.size() * 4;
	double p_in_point, qi_local[4];

	//собираем глобальные номера с элемента
	for(int j = 0; j < 4; j++)
		indexes[j] = P.elements[element_number].nodes[j] + n_edges;

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
	int size = P.elements.size();

	for(int i = 0; i < size; i++)
	{
		x_left = P.nodes[P.elements[i].nodes[0]].x;
		x_right = P.nodes[P.elements[i].nodes[1]].x;
		y_low = P.nodes[P.elements[i].nodes[0]].y;
		y_up = P.nodes[P.elements[i].nodes[3]].y;
		if(x_left <= x && x <= x_right && y_low <= y && y <= y_up)
			return i;
	}
}

void SLAE::get_vector_solution_in_nodes_ux(MyVector qi, MyVector &solution)
{
	int indexes[4], indexes_nodes[4];
	int size = P.elements.size();
	double u_local[4], qi_local[4];
	double x, y;

	for(int i = 0; i < size; i++)
	{
		//собираем глобальные номера с элемента
		for(int j = 0; j < 4; j++)
			indexes[j] = P.elements[i].edges[j];

		for(int j = 0; j < 4; j++)
			indexes_nodes[j] = P.elements[i].nodes[j];

		//собираем локальный набор весов
		for(int j = 0; j < 4; j++)
			qi_local[j] = qi[indexes[j]];

		//вычисляем в узлах элемента решение
		for(int j = 0; j < 4; j++)
		{
			u_local[j] = 0;
			x = P.nodes[indexes_nodes[j]].x;
			y = P.nodes[indexes_nodes[j]].y;
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
	int size = P.elements.size();
	double u_local[4], qi_local[4];
	double x, y;

	for(int i = 0; i < size; i++)
	{
		//собираем глобальные номера с элемента
		for(int j = 0; j < 4; j++)
			indexes[j] = P.elements[i].edges[j];

		//собираем локальный набор весов
		for(int j = 0; j < 4; j++)
			qi_local[j] = qi[indexes[j]];

		for(int j = 0; j < 4; j++)
			indexes_nodes[j] = P.elements[i].nodes[j];

		//вычисляем в узлах элемента решение
		for(int j = 0; j < 4; j++)
		{
			u_local[j] = 0;
			x = P.nodes[indexes_nodes[j]].x;
			y = P.nodes[indexes_nodes[j]].y;
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
	int indexes[4], n_edges = P.elements.size() * 4;
	int size = P.elements.size();
	double p_local[4], qi_local[4];
	double x, y;

	for(int i = 0; i < size; i++)
	{
		//собираем глобальные номера с элемента
		for(int j = 0; j < 4; j++)
			indexes[j] = P.elements[i].nodes[j];

		//собираем локальный набор весов
		for(int j = 0; j < 4; j++)
			qi_local[j] = qi[indexes[j]+ n_edges];

		//вычисляем в узлах элемента решение
		for(int j = 0; j < 4; j++)
		{
			p_local[j] = 0;
			x = P.nodes[indexes[j]].x;
			y = P.nodes[indexes[j]].y;
			for(int k = 0; k < 4; k++)
				p_local[j] += qi_local[k] * psi_i(k, x, y, i);
		}

		//кладём в результирующий вектор
		for(int j = 0; j < 4; j++)
			solution[indexes[j]] = p_local[j];
	}
}

double SLAE::phix_i(int i, double x, double y, int element_number)
{
	double x_left = P.nodes[P.elements[element_number].nodes[0]].x;
	double x_right = P.nodes[P.elements[element_number].nodes[1]].x;
	double y_low = P.nodes[P.elements[element_number].nodes[0]].y;
	double y_up = P.nodes[P.elements[element_number].nodes[3]].y;
	double hx = x_right - x_left, hy = y_up - y_low;
	double ksi = 2 * (x - (x_left + x_right) / 2) / hx, 
		   etta = 2 * (y - (y_low + y_up) / 2) / hy;
	return phix[i](ksi, etta);
}

double SLAE::phiy_i(int i, double x, double y, int element_number)
{
	double x_left = P.nodes[P.elements[element_number].nodes[0]].x;
	double x_right = P.nodes[P.elements[element_number].nodes[1]].x;
	double y_low = P.nodes[P.elements[element_number].nodes[0]].y;
	double y_up = P.nodes[P.elements[element_number].nodes[3]].y;
	double hx = x_right - x_left, hy = y_up - y_low;
	double ksi = 2 * (x - (x_left + x_right) / 2) / hx, 
		   etta = 2 * (y - (y_low + y_up) / 2) / hy;
	return phiy[i](ksi, etta);
}

double SLAE::phixdx_i(int i, double x, double y, int element_number)
{
	double x_left = P.nodes[P.elements[element_number].nodes[0]].x;
	double x_right = P.nodes[P.elements[element_number].nodes[1]].x;
	double y_low = P.nodes[P.elements[element_number].nodes[0]].y;
	double y_up = P.nodes[P.elements[element_number].nodes[3]].y;
	double hx = x_right - x_left, hy = y_up - y_low;
	double ksi = 2 * (x - (x_left + x_right) / 2) / hx, 
		   etta = 2 * (y - (y_low + y_up) / 2) / hy;
	return dphixksi[i](ksi, etta) / hx;
}

double SLAE::phiydy_i(int i, double x, double y, int element_number)
{
	double x_left = P.nodes[P.elements[element_number].nodes[0]].x;
	double x_right = P.nodes[P.elements[element_number].nodes[1]].x;
	double y_low = P.nodes[P.elements[element_number].nodes[0]].y;
	double y_up = P.nodes[P.elements[element_number].nodes[3]].y;
	double hx = x_right - x_left, hy = y_up - y_low;
	double ksi = 2 * (x - (x_left + x_right) / 2) / hx, 
		   etta = 2 * (y - (y_low + y_up) / 2) / hy;
	return dphiyetta[i](ksi, etta) / hy;
}

double SLAE::psi_i(int i, double x, double y, int element_number)
{
	double x_left = P.nodes[P.elements[element_number].nodes[0]].x;
	double x_right = P.nodes[P.elements[element_number].nodes[1]].x;
	double y_low = P.nodes[P.elements[element_number].nodes[0]].y;
	double y_up = P.nodes[P.elements[element_number].nodes[3]].y;
	double hx = x_right - x_left, hy = y_up - y_low;
	double ksi = 2 * (x - (x_left + x_right) / 2) / hx, 
		   etta = 2 * (y - (y_low + y_up) / 2) / hy;
	return psi[i](ksi, etta);
}

#pragma region предобусловливание

void SLAE::LU()
{
	int i;
	int i0,j0;
	int iend;
	int num;
	int ki,kj;
	double suml,sumu,sumdg;
	int size2 = A.size;

	for(i = 0; i < size2; i++) 
	{
		LU_ggu[i] = A.ggu[i];
		LU_ggl[i] = A.ggl[i];
	}

	for(i = 0; i < n; i++) 
		LU_di[i] = A.di[i];

	for(i = 0; i < n;i++)
	{
		i0 = A.ig[i];
		iend = A.ig[i+1];
		for(num = i0,sumdg = 0; num < iend; num++)
		{
		    j0 = A.ig[A.jg[num]]; //в зависимости от номера фиксируем столбец,какой столбец l,такого столбца  ищем начальный эл у u 
			int jend=A.ig[A.jg[num]+1];
			ki=i0;
			kj=j0;
			for(suml = 0, sumu = 0, ki = i0; ki < num; ki++) //для num учитываются все предыдущие элементы
			{
				for(int m = kj; m < jend; m++)
				if(A.jg[ki]==A.jg[m]) //ищем соответствующие ненулевые элементы для умножения
				{
					suml += LU_ggl[ki] * LU_ggu[m];
					sumu += LU_ggl[m] * LU_ggu[ki];//для симметричного элемента из U
				}
			}
			LU_ggl[num] -= suml;	
			LU_ggu[num] = (LU_ggu[num] - sumu) / LU_di[A.jg[num]];
		sumdg += LU_ggl[num] * LU_ggu[num];//умножаются симметричные элементы	
		}
		LU_di[i]-=sumdg;
	}
}

void SLAE::LYF(MyVector b)
{
	int i, k;
	int i0;//адрес начала строки
	int iend;//адрес конца строки
	double sum;

	yl.make_zero();

	if(use_LU)
	{
		for(i = 0; i < n; i++)
		{
			i0 = A.ig[i]; iend = A.ig[i+1];

			for(k = i0, sum = 0; k < iend; k++)
				sum += yl[A.jg[k]] * LU_ggl[k];

			yl[i] = (b[i] - sum) / LU_di[i];
		}
	}
	else
		yl = b;
}

void SLAE::LYFt(MyVector b)
{
	int i, k;
	int i0;//адрес начала строки
	int iend;//адрес конца строки
	double sum;

	yl.make_zero();
	if(use_LU)
	{
		MyVector bb(n);
		bb = b;
		for(i = n - 1; i >= 0; i--)
		{
			i0 = A.ig[i]; iend = A.ig[i+1];
			yl[i] = bb[i] /LU_di[i];
			for(k = i0, sum = 0; k < iend; k++)
				bb[A.jg[k]] -= yl[i] * LU_ggl[k];
		}
	}
	else
		yl = b;
}

void SLAE::UXY(MyVector b)
{
	int i, k;
	int i0;
	int iend;

	yu.make_zero();
	if(use_LU)
	{
		for(i = n - 1; i >= 0; i--)//проход по столбцам с конца
		{
			yu[i] += b[i];
			i0 = A.ig[i]; iend = A.ig[i+1];

			for(k = iend - 1; k >= i0; k--)//идём по столбцу с конца
				yu[A.jg[k]] -= yu[i] * LU_ggu[k];
		}
	}
	else
		yu = b;
}

void SLAE::UXYt(MyVector b)
{
	int i, k;
	int i0;
	int iend;

	yu.make_zero();
	if(use_LU)
	{
		for(i = n - 1; i >= 0; i--)//проход по столбцам с конца
		{
			yu[i] += b[i];
			i0 = A.ig[i]; iend = A.ig[i+1];

			for(k = iend - 1; k >= i0; k--)//идём по столбцу с конца
				yu[i] -= yu[A.jg[k]] * LU_ggu[k];
		}
	}
	else
		yu = b;
}

void SLAE::convert_to_prof()
{
	int i, j, k;
	vector <double> ggu_new;
	vector <double> ggl_new;
	vector <int>  ig_new;
	for (i = 1; i < n; i++)	//	идём по строкам матрицы
	{
		for (j = A.ig[i]; j < A.ig[i + 1] - 1; j++)	//	идём по элементам строки
		{
			ggu_new.push_back(A.ggu[j]);
			ggl_new.push_back(A.ggl[j]);
			for (k = 1; A.jg[j] + k < A.jg[j + 1]; k++)	//	определяем сколько не хватает нулей для перевода
			{
				ggu_new.push_back(0);
				ggl_new.push_back(0);
			}
		}
		ggu_new.push_back(A.ggu[j]);
		ggl_new.push_back(A.ggl[j]);
		for (k = 1; A.jg[j] + k < i; k++)
		{
			ggu_new.push_back(0);
			ggl_new.push_back(0);
		}
	}
	size_prof = ggu_new.size();

	ig_new.push_back(0);
	ig_new.push_back(0);
	for (i = 1; i < n; i++)	//	идём по строкам матрицы	
		ig_new.push_back(ig_new[i] + i - A.jg[A.ig[i]]);

	LU_ggl2.resize(size_prof);
	LU_ggu2.resize(size_prof);
	LU_di2.resize(n);
	LU_ig2.resize(n + 1);

	LU_ggu2 = ggu_new;
	LU_ggl2 = ggl_new;

	LU_di2 = A.di;
	LU_ig2 = ig_new;
}

void SLAE::LU2()
{
	int i,j;
	int i0,j0;
	int iend, jend;
	int kol_i,kol_j;
	int num;
	int first_i,first_j;
	int ki,kj;
	int dif;
	int d = 0;
	double suml, sumu, sumdg;//суммы для элементов в L,U

	for(i = 0; i < n; i++)
	{
		i0 = LU_ig2[i];//начало i строки (столбца)
		iend = LU_ig2[i + 1];//начало i+1 строки (столбца)
		kol_i = iend - i0;//количество элементов в i строке (столбце)
		first_i = i - kol_i;//номер первого ненулевого элемента i строки (столбца)
		j = first_i;

		for(num = i0, sumdg = 0; j < i; num++, j++)//идём по всем ненулевым элементам строки (столбца)
		{
			j0 = LU_ig2[j];
			jend = LU_ig2[j + 1];
			kol_j = jend - j0;
			first_j = j - kol_j;
			ki = i0;
			kj = j0;
			dif = first_j - first_i;//разница в профилях
			//смещение по профилю для корректного умножения
			if(dif > 0) ki += dif; //ki - номер первого умножаемого элемента в L,если идём по строке
			else kj -= dif;		 //kj - номер первого умножаемого элемента в U,если идём по строке

			for(suml = 0, sumu = 0; ki < num; ki++, kj++) //для num учитываются все предыдущие элементы
			{
				suml += LU_ggl2[ki] * LU_ggu2[kj];
				sumu += LU_ggl2[kj] * LU_ggu2[ki];//считается сумма для симметричного элемента из U
			}

			LU_ggl2[num] -= suml;
			LU_ggu2[num] = (LU_ggu2[num] - sumu) / LU_di2[j];
			sumdg += LU_ggl2[num] * LU_ggu2[num];//умножаются симметричные элементы
		}
		LU_di2[i] -= sumdg; 
	}
}

void SLAE::LYF2(MyVector b)
{
	int i,j,k;
	int i0;//адрес начала строки
	int iend;//адрес конца строки
	int ikol;//количество элементов в строке
	int beg;//номер столбца первого ненулевого элемента строки
	double sum;

	yl.make_zero();

	for(i = 0; i < n; i++)
	{
		i0 = LU_ig2[i]; iend = LU_ig2[i + 1];
		ikol = iend - i0;
		beg = i - ikol;

		for(k = i0, sum = 0, j = beg; j < i; j++, k++)
		{
			sum += yl[j] * LU_ggl2[k];
		}

		yl[i] = (b[i] - sum) / LU_di2[i];
	}
}

void SLAE::UXY2(MyVector b)
{
	int i,j,k;
	int i0;
	int iend;
	int ikol;
	int beg;

	yu.make_zero();

	for(i = n - 1; i >= 0; i--)//проход по столбцам с конца
	{
		yu[i] += yl[i];
		i0 = LU_ig2[i]; iend = LU_ig2[i+1];
		ikol = iend - i0;
		beg = i - ikol;

		for(k = iend-1, j = i-1; j >= beg; j--, k--)//идём по столбцу с конца
		{
			yu[j] -= yu[i] * LU_ggu2[k];
		}
	}
}


MyVector SLAE::Uv(MyVector v)
{
	int i, j, k, kol;
	int iend;
	MyVector new_vector = MyVector(v.ar.size());

	assert(v.ar.size() == n);
	return v;
	for(i = 0; i < n; i++)
	{
		kol = A.ig[i+1] - A.ig[i];//количество ненулевых элементов столбца от первого
								//ненулевого элемента до диагонального элемента (не включая его)
		iend = A.ig[i+1];
		k = A.ig[i]; // адрес первого занятого элемента столбца

		new_vector[i] = v[i];//от главной диагонали (у U на диагонали 1)

		for(; k < iend; k++)//проходим по всем элементам i столбца
		{
			j = A.jg[k];
			new_vector[j] += LU_ggu[k] * v[i];//от верхнего треугольника
		}
	}

	return new_vector;
}

#pragma endregion

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
	int k_iter;

	x = U_begin;
	f = A.b;
	LU();
	r = f - A * x;
	LYF(r); r = yl;
	norm_r = r.norm();
	LYF(f); norm_f = yl.norm();
	x = Uv(U_begin);

	for(int k_iter = 1; k_iter <= max_iter && norm_r / norm_f > eps; k_iter++)
	{
		d.make_zero();
		V[0] = r / norm_r;

		continue_ = true;
		for(int j = 1; j <= m && continue_; j++)
		{
			UXY(V[j - 1]);
			tmp = A * yu;
			LYF(tmp); w = yl;

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
		UXY(x);
		tmp = f - A * yu;
		LYF(tmp); r = yl;
		norm_r = r.norm();
		logger.send_current_information(norm_r / norm_f, k_iter);
		printf("%d\tr=%.10e\n", k_iter, norm_r / norm_f);
	}
	UXY(x); x = yu;
	solution = x;
}

void SLAE::BCGStab(MyVector U_begin, MyVector &solution)
{
	int k_it;
	double  rkr0, ak, gk, bk;
	MyVector r(n), f(n), x(n), r0(n), z(n), p(n), v(n), v1(n), rr2(n);
	double r_norm, f_norm;

	LU();

	x = U_begin;
	f = A.b;
	f_norm = f.norm();
	r0 = f - A * x;
	LYF(r0); r0 = yl;
	r_norm = r0.norm() / f_norm;

	UXY(r0); z = yu;
	r = r0;

	logger.send_current_information(r_norm, 0);

	for(int k_it = 1; k_it <= max_iter && r_norm > eps; k_it++)
	{
		//найдем L^(-1)AU^(-1)zk
		UXY(z); v = yu;// v = U(-1)zk
		v1 = A * v; // v1 = AU^(-1)zk

		LYF(v1); v = yl;// v = L^(-1)AU^(-1)zk

		rkr0 = scal(r, r0);
		ak = rkr0 / scal(v, r0); // ak = (r,r0)/ ( L^(-1)AU^(-1)zk,r0)

		p = r - v * ak; // pk = r - ak*L^(-1)AU^(-1)zk

		//найдем L^(-1)AU^(-1)pk
		UXY(p); v1 = yu; // v1 = U^(-1)pk
		LYF(A * v1); v1 = yl; // v1 = L^(-1)AU^(-1)pk

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
	UXY(x); x = yu;
	solution = x;
}

void SLAE::BCG(MyVector U_begin, MyVector &solution)
{
    MyVector r(n), r_(n), p(n), p_(n), f(n);
    MyVector v1(n), v2(n), v3(n), x(n);
    double alpha, betta, old_r_norm = 1.e+20, sc1, sc2;
	double r_norm, r_norm_;
	int k_it;

	LU();

	f = A.b;

    k_it = 0;
    r_norm = old_r_norm / 10;
    r_norm_ = 0;

	x = U_begin;
	v1 = A * x;
	v1 = f - v1;
    
	LYF(v1); r = yl;
	r_ = r;
	p = r;
	p_ = r_;
	int flag = 0;

    while(flag == 0 && k_it < max_iter)
    {
        sc1 = scal(r, r_);
		UXY(p); v1 = yu;
		v2 = A * v1;
		LYF(v2); v3 = yl;

        sc2 = scal(p_, v3);

        alpha = sc1 / sc2;
		x = x +  v1 * alpha;
		r = r - v3 * alpha;

		LYFt(p_); v1 = yl;

		v2 = A  /  v1;

		UXYt(v2); v3 = yu;

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

	switch(solver)
	{
	case 0:
		{
			convert_to_prof();
			LU2();
			LYF2(A.b);
			UXY2(yl);
			q = yu;
		}
		break;
	case 1:
		{
			use_LU = false;
			BCGStab(U_begin, q);
		}
		break;
	case 2:
		{
			use_LU = false;
			GMRES(U_begin, q);
		}
		break;
	case 3:
		{
			use_LU = false;
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

void SLAE::si_print(FILE *log_f, int iteration_number, double &normL2u, double &normL2p)
{
	FILE *solution_f_out, *info_f_out;
	string f_name_s, f_name_i;

	fprintf(log_f, "---%d---\n", iteration_number);
	f_name_s = string("s_") + to_string(iteration_number) + ".txt";
	f_name_i = string("i_") + to_string(iteration_number) + ".txt";
	solution_f_out = fopen(f_name_s.c_str(), "w");
	info_f_out = fopen(f_name_i.c_str(), "w");
	output(solution_f_out, info_f_out, normL2u, normL2p);
	fclose(solution_f_out);
	fclose(info_f_out);
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
	double w, residual, residual_previous;
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
	int size = P.elements.size();
	double p, function;
	double x0, y0, hx, hy;
	double jacobian;

	diff = 0;
	for(int i = 0; i < size; i++)
	{
		x0 = P.nodes[P.elements[i].nodes[0]].x;
		y0 = P.nodes[P.elements[i].nodes[0]].y;
		hx = get_hx(i);
		hy = get_hy(i);
		jacobian = hx * hy / 4.0;

		diff_local = 0;
		for(int k = 0; k < 9; k++)
		{
			double p_x = hx * gauss_points[0][k] + x0;
			double p_y = hy * gauss_points[1][k] + y0;
			p = get_solution_in_point_p(p_x, p_y, i, q_solution);
			function = calculate_p_analytic(P.elements[i].number_of_area, p_x, p_y);
			function -= p;
			diff_local += gauss_weights[k] * function * function;
		}
		diff_local *= jacobian;

		diff += diff_local;
	}

	return sqrt(diff / P.nodes.size());
}

double SLAE::diff_normL2_u(MyVector q_solution)
{
	double diff_local, diff;
	int size = P.elements.size();
	double function;
	double x0, y0, hx, hy;
	double jacobian;
	double ux, uy, uxdx, uydy;
	double ux_an, uy_an, uxdx_an, uydy_an;

	diff = 0;
	for(int i = 0; i < size; i++)
	{
		x0 = P.nodes[P.elements[i].nodes[0]].x;
		y0 = P.nodes[P.elements[i].nodes[0]].y;
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
			ux_an = calculate_ux_analytic(P.elements[i].number_of_area, p_x, p_y);
			uy_an = calculate_uy_analytic(P.elements[i].number_of_area, p_x, p_y);
			uxdx_an = calculate_uxdx_analytic(P.elements[i].number_of_area, p_x, p_y);
			uydy_an = calculate_uydy_analytic(P.elements[i].number_of_area, p_x, p_y);

			function = (ux - ux_an) * (ux - ux_an) + (uy - uy_an) * (uy - uy_an);
			function += (uxdx - uxdx_an + uydy - uydy_an) * 
						(uxdx - uxdx_an + uydy - uydy_an);

			diff_local += gauss_weights[k] * function;
		}
		diff_local *= jacobian;

		diff += diff_local;
	}

	return sqrt(diff / P.nodes.size());
}

void SLAE::run(FILE *solution_f_out, FILE *info_f_out)
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

void main()
{
	FILE *grid_f_in, *elements_f_in, *l1_f_in;
	FILE *solution_f_out, *info_f_out, *log_f;

	grid_f_in = fopen("grid.txt", "r");
	elements_f_in = fopen("elements.txt", "r");
	l1_f_in = fopen("l1.txt", "r");
	solution_f_out = fopen("solution.txt", "w");
	info_f_out = fopen("info.txt", "w");
	log_f = fopen("log.txt", "w");

	SLAE my_SLAE = SLAE(1000, 100, 1e-12, 30, grid_f_in, elements_f_in, log_f, l1_f_in);
	//my_SLAE.run(solution_f_out, info_f_out);
	my_SLAE.simple_iterations();

	fclose(grid_f_in);
	fclose(elements_f_in);
	fclose(solution_f_out);
	fclose(info_f_out);
	fclose(log_f);
	_getch();
}
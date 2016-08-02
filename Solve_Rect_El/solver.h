#pragma once
#include <stdio.h>
#include <conio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <vector>
#include <assert.h>
#include <algorithm>
#include <functional>
#include <string>
#include <fstream>
#include <iomanip>
using namespace std;

#include "boundary_conditions.h"
#include "element.h"
#include "partition.h"
#include "point.h"
#include "myvector.h"
#include "myfunctions.h"
#include "matrix.h"
#include "testing_parameters.h"
#include "parameters.h"
#include "densematrix.h"
#include "boundaries.h"
#include "logger.h"

using namespace boundary_conditions;
using namespace element;
using namespace partition;
using namespace point;
using namespace myvector;
using namespace matrix;
using namespace densematrix;
using namespace boundaries;
using namespace logger;
using namespace parameters;


struct SLAE : public BoundaryConditionsSupport, public InternalBoundaries, public OuterBoundaries
{
	int n; //размерность СЛАУ
	int m; //глубина метода gmres
	int max_iter; //max количество итераций
	int max_iter_nonlinear;
	double eps; //точность решения СЛАУ
	Matrix A; //матрица
	Logger logger; //логгер для вывода информации о процессе решения СЛАУ


	//S для скорости
	//E для скорости
	//P_ //давление по j
	//P__ //давление по i
	//SP для давления

	//локальные матрицы
	double G[4][4], P1[4][4], P2[4][4], C[4][4]; 
	double F[4]; //локальный вектор правой части

	MyVector Ux_numerical; //численное решение Ux
	MyVector Uy_numerical; //численное решение Uy
	MyVector P_numerical; //численное решение P

	MyVector q_prev; //вектор весов с предыдущей итерации по нелинейности

	vector <double> LU_ggu2; //верхнетреугольные недиагональные элементы U для LU-решателя
	vector <double> LU_ggl2; //нижнетреугольные недиагональные элементы L для LU-решателя
	vector <double> LU_di2; //диагональные элементы L для LU-решателя
	vector <int> LU_ig2; //индексы портрета для LU-решателя
	int size_prof; //размер массивов элементов для профильного формата LU-решателя
    
	SLAE(){};

	SLAE(int max_number_of_iterations, 
		 int max_number_of_iterations_non_lin,
		 double epsilon, 
		 int gmres_m, 
		 ifstream& grid_f_in, 
		 ifstream& elements_f_in,
		 string log_f,
		 ifstream& boundary1)
	{
		initialize(max_number_of_iterations, 
				   max_number_of_iterations_non_lin,
				   epsilon,
				   gmres_m,
				   grid_f_in,
				   elements_f_in,
				   log_f,
				   boundary1);
	}

	~SLAE(){};

	void initialize(int max_number_of_iterations,
					int max_number_of_iterations_non_lin, 
					double epsilon,
					int gmres_m,
					ifstream& grid_f_in,
					ifstream& elements_f_in,
					string log_f,
					ifstream& boundary1);

	void reinitialize();

	int count_unzero_matrix_elements();
	int create_unzero_elements_list(int element_number, 
									vector <int> &list, 
									int dof_num_i, 
									int dof_num_j, 
									int *dof_i, 
									int *dof_j,
									bool dof_j_edge);
	void create_portret();

	void calculate_global_matrix(MyVector q_calc);

	double get_solution_in_point_ux(double x, double y, int element_number, MyVector qi);
	double get_solution_in_point_uy(double x, double y, int element_number, MyVector qi);
	double get_solution_in_point_uxdx(double x, double y, int element_number, MyVector qi);
	double get_solution_in_point_uydy(double x, double y, int element_number, MyVector qi);
	double get_solution_in_point_p(double x, double y, int element_number, MyVector qi);

	double get_solution_in_point2_ux(double x, double y, MyVector qi);
	double get_solution_in_point2_uy(double x, double y, MyVector qi);
	double get_solution_in_point2_p(double x, double y, MyVector qi);

	int search_element(double x, double y);

	void get_vector_solution_in_nodes_ux(MyVector qi, MyVector &solution);
	void get_vector_solution_in_nodes_uy(MyVector qi, MyVector &solution);
	void get_vector_solution_in_nodes_p(MyVector qi, MyVector &solution);

	//локальные матрицы и векторы
	void calculate_locals(int element_number, MyVector q_calc);
	void calculate_G(int element_number);
	void calculate_C(int element_number, MyVector q_calc);
	void calculate_P1(int element_number);
	void calculate_P2(int element_number);
	void calculate_F(int element_number);

	//решатели
	void solve_min_sqr_problem(MyVector d, DenseMatrix H, MyVector &result);
	void GMRES(MyVector U_begin, MyVector &solution);

	void BCGStab(MyVector U_begin, MyVector &solution);
	void BCG(MyVector U_begin, MyVector &solution);

	void convert_to_prof();
	void LU2();
	void LYF2(MyVector b);
	void UXY2(MyVector b);

	void Solve(MyVector U_begin, double &normL2u, double &normL2p);

	void si_print(ofstream& log_f, int iteration_number, double &normL2u, double &normL2p);
	double find_relaxation_parameter(MyVector q_current, MyVector q_previous, double &residual_previous);
	void simple_iterations();

	double SLAE::diff_normL2_p(MyVector q_solution);//погрешность решения в норме L2
	double SLAE::diff_normL2_u(MyVector q_solution);//погрешность решения в норме L2


	//запуск решения
	void run(ofstream& solution_f_out, ofstream& info_f_out);

	void output(ofstream& solution_f_out, ofstream& info_f_out, double normL2u, double normL2p)
	{		
		ofstream res("result.txt");

		int n_nodes = nodes.size();
		for(int i = 0; i < n_nodes; i++)
		{
			if(abs(nodes[i].x - 0.5) < 1e-10)
				res << nodes[i].x << "\t" << nodes[i].y << "\t" << Ux_numerical[i]
				<< "\t" << Uy_numerical[i] << "\t" << P_numerical[i] << "\t"
				<< calculate_p_analytic(0, nodes[i].x, nodes[i].y) << endl;
		}

		res.close();

		solution_f_out << Ux_numerical << endl << endl << endl;
		solution_f_out << Uy_numerical << endl << endl << endl;
		solution_f_out << P_numerical << endl << endl << endl;
		info_f_out << "norm L2 u:|u*-u|=" << scientific << setprecision(4) << normL2u 
			<< endl << "norm L2 p:|p*-p|=" << scientific << setprecision(4) << normL2p
			<< endl << "eps=" << scientific << setprecision(2) << eps << endl;
	};
};


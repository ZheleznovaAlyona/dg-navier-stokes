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
#include "solver.h"

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
using namespace solver;

namespace slae
{
struct SLAE : public BoundaryConditionsSupport, public InternalBoundaries, public OuterBoundaries, public Solver
{
	int n; //размерность СЛАУ
	int m; //глубина метода gmres
	int max_iter; //max количество итераций
	int max_iter_nonlinear;
	double eps; //точность решения СЛАУ
	Matrix A;
	myvector::MyVector b;//вектор правой части
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

	void calculate(MyVector q_calc, 
				   partition::Partition& p,
				   boundaries::InternalBoundaries& internal_bs,
				   boundaries::OuterBoundaries& outer_bs,
				   boundary_conditions::BoundaryConditionsSupport&
				   b_conditions);
							    
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

	double SLAE::diff_normL2_p(MyVector q_solution);//погрешность решения в норме L2
	double SLAE::diff_normL2_u(MyVector q_solution);//погрешность решения в норме L2
};
}


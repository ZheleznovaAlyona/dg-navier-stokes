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
	int n; //����������� ����
	int m; //������� ������ gmres
	int max_iter; //max ���������� ��������
	int max_iter_nonlinear;
	double eps; //�������� ������� ����
	Matrix A;
	Logger logger; //������ ��� ������ ���������� � �������� ������� ����


	//S ��� ��������
	//E ��� ��������
	//P_ //�������� �� j
	//P__ //�������� �� i
	//SP ��� ��������

	//��������� �������
	double G[4][4], P1[4][4], P2[4][4], C[4][4]; 
	double F[4]; //��������� ������ ������ �����

	MyVector Ux_numerical; //��������� ������� Ux
	MyVector Uy_numerical; //��������� ������� Uy
	MyVector P_numerical; //��������� ������� P

	MyVector q_prev; //������ ����� � ���������� �������� �� ������������
    
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

	//��������� ������� � �������
	void calculate_locals(int element_number, MyVector q_calc);
	void calculate_G(int element_number);
	void calculate_C(int element_number, MyVector q_calc);
	void calculate_P1(int element_number);
	void calculate_P2(int element_number);
	void calculate_F(int element_number);

	double SLAE::diff_normL2_p(MyVector q_solution);//����������� ������� � ����� L2
	double SLAE::diff_normL2_u(MyVector q_solution);//����������� ������� � ����� L2


	////������ �������
	//void run(ofstream& solution_f_out, ofstream& info_f_out);

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
}


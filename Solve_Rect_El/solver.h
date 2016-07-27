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
#include "densematrix.h"

using namespace boundary_conditions;
using namespace element;
using namespace partition;
using namespace point;
using namespace myvector;
using namespace matrix;
using namespace densematrix;


struct Logger
{
	ofstream log_f;
	void send_current_information(double r_norm, int iteration_number)
	{
		log_f << iteration_number << "     " << setprecision(20) << r_norm << endl;
	};
};

struct SLAE : public BoundaryConditionsSupport
{
	int n; //����������� ����
	int m; //������� ������ gmres
	int max_iter; //max ���������� ��������
	int max_iter_nonlinear;
	double eps; //�������� ������� ����
	double sigma, mu1, mu2; //������������ ������������
	Matrix A; //�������
	Logger logger; //������ ��� ������ ���������� � �������� ������� ����


	//S ��� ��������
	//E ��� ��������
	//P_ //�������� �� j
	//P__ //�������� �� i
	//SP ��� ��������

	//��������� �������
	double G[4][4], P1[4][4], P2[4][4], C[4][4]; 
	double E[8][8], P_1[8][8], P_2[8][8], SP[8][8]; //��������� �������
	double E_out[4][4], P_1_out[4][4], P_2_out[4][4], SP_out[4][4]; //��������� �������
	double F[4]; //��������� ������ ������ �����

	MyVector Ux_numerical; //��������� ������� Ux
	MyVector Uy_numerical; //��������� ������� Uy
	MyVector P_numerical; //��������� ������� P

	MyVector q_prev; //������ ����� � ���������� �������� �� ������������

	function<double(double, double)> psi[4]; //��������� �� ������� ���������� �������� ������� p � �����
	function<double(double, double)> dpsiksi[4]; //��������� �� ������� ���������� d/dksi �������� ������� p � �����
	function<double(double, double)> dpsietta[4]; //��������� �� ������� ���������� d/detta �������� ������� p � �����
	function<double(double, double)> phix[4]; //��������� �� ������� ���������� �������� ������� ux � �����
	function<double(double, double)> phiy[4]; //��������� �� ������� ���������� �������� ������� uy � �����
	function<double(double, double)> dphixksi[4]; //��������� �� ������� ���������� d/dksi �������� ������� ux � �����
	function<double(double, double)> dphixetta[4]; //��������� �� ������� ���������� d/detta �������� ������� ux � �����
	function<double(double, double)> dphiyksi[4]; //��������� �� ������� ���������� d/dksi �������� ������� uy � �����
	function<double(double, double)> dphiyetta[4]; //��������� �� ������� ���������� d/detta �������� ������� uy � �����

	vector <double> LU_ggu2; //����������������� �������������� �������� U ��� LU-��������
	vector <double> LU_ggl2; //���������������� �������������� �������� L ��� LU-��������
	vector <double> LU_di2; //������������ �������� L ��� LU-��������
	vector <int> LU_ig2; //������� �������� ��� LU-��������
	int size_prof; //������ �������� ��������� ��� ����������� ������� LU-��������
	
	double gauss_points[2][9];//����� ������
	double gauss_weights[9];// ���� ������
	double gauss_points_1[3];//����� ������
	double gauss_weights_1[3];// ���� ������
    
	SLAE(){};

	SLAE(int max_number_of_iterations, 
		 int max_number_of_iterations_non_lin,
		 double epsilon, 
		 int gmres_m, 
		 ifstream& grid_f_in, 
		 ifstream& elements_f_in, 
		 ofstream& log_f, 
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
					ofstream& log_f,
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

	void add_element_to_global_matrix(int i, int j, double element);
	void put_element_to_global_matrix(int i, int j, double element);

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

	double phix_i(int i, double x, double y, int element_number);
	double phiy_i(int i, double x, double y, int element_number);
	double phixdx_i(int i, double x, double y, int element_number);
	double phiydy_i(int i, double x, double y, int element_number);
	double psi_i(int i, double x, double y, int element_number);

	//��������� ������� � �������
	void calculate_locals(int element_number, MyVector q_calc);
	void calculate_G(int element_number);
	void calculate_C(int element_number, MyVector q_calc);
	void calculate_P1(int element_number);
	void calculate_P2(int element_number);
	void calculate_F(int element_number);

	//��������� ������ �� ���������� ��������
	void calculate_internal_boundaries(int element_number);

	void calculate_ES_horizontal(int element_number1, int element_number2);
	void calculate_ES_vertical(int element_number1, int element_number2);

	void calculate_P_1_horizontal(int element_number1, int element_number2);
	void calculate_P_1_vertical(int element_number1, int element_number2);

	void calculate_P_2_horizontal(int element_number1, int element_number2);
	void calculate_P_2_vertical(int element_number1, int element_number2);

	void calculate_SP_horizontal(int element_number1, int element_number2);
	void calculate_SP_vertical(int element_number1, int element_number2);

	void add_ES_to_global(int element_number, int neighbor_element_number);
	void add_P_1_to_global(int element_number, int neighbor_element_number);
	void add_P_2_to_global(int element_number, int neighbor_element_number);
	void add_SP_to_global(int element_number, int neighbor_element_number);


	//��������� ������ �� ������� ��������
	void calculate_outer_boundaries(int element_number);

	void calculate_ES_out_left(int element_number);
	void calculate_ES_out_right(int element_number);
	void calculate_ES_out_low(int element_number);
	void calculate_ES_out_up(int element_number);

	void calculate_P_1_out_left(int element_number);
	void calculate_P_1_out_right(int element_number);
	void calculate_P_1_out_low(int element_number);
	void calculate_P_1_out_up(int element_number);

	void calculate_P_2_out_left(int element_number);
	void calculate_P_2_out_right(int element_number);
	void calculate_P_2_out_low(int element_number);
	void calculate_P_2_out_up(int element_number);

	void calculate_SP_out_left(int element_number);
	void calculate_SP_out_right(int element_number);
	void calculate_SP_out_low(int element_number);
	void calculate_SP_out_up(int element_number);

	double gx(int formula_number, double x, double y);
	double gy(int formula_number, double x, double y);


	//������� � ��������� ��� ������
	double calculate_fx(int area_number, double x, double y);
	double calculate_fy(int area_number, double x, double y);

	double calculate_ux_analytic(int area_number, double x, double y);
	double calculate_uy_analytic(int area_number, double x, double y);
	double calculate_uxdx_analytic(int area_number, double x, double y);
	double calculate_uydy_analytic(int area_number, double x, double y);

	double calculate_p_analytic(int area_number, double x, double y);

	double calculate_lambda(int area_number);
	double calculate_rho(int area_number);


	//��������
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

	double SLAE::diff_normL2_p(MyVector q_solution);//����������� ������� � ����� L2
	double SLAE::diff_normL2_u(MyVector q_solution);//����������� ������� � ����� L2


	//������ �������
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


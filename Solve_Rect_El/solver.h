#include <stdio.h>
#include <conio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <vector>
#include <assert.h>
#include <algorithm>
#include <functional>
#include <string.h>
#include <string>
using namespace std;

bool use_LU;
int test = 3;
int solver = 2;

struct Point
{
	double x;
	double y;
	bool operator==(Point point)
	{
		if(point.x == x && point.y == y)
			return true;
		else
			return false;
	}

};

struct Element
{
	int nodes[4];
	int edges[4];
	int number_of_area;
	int neighbors[4]; //�����, ������, ������, �������

	Element& operator=(Element element)
	{
		for(int i = 0; i < 4; i++)
			nodes[i] = element.nodes[i];
		for(int i = 0; i < 4; i++)
			edges[i] = element.edges[i];
		number_of_area = element.number_of_area;
		for(int i = 0; i < 4; i++)
		neighbors[i] = element.neighbors[i];

		return *this;
	}
};

struct Partition
{
	vector <Element> elements;
	vector <Point> nodes;
	void input(FILE *grid_f_in, FILE *elements_f_in);
};

struct MyVector
{
	vector <double> ar;

	MyVector(){};

	MyVector(int size)
	{
		ar.reserve(size);
		for(int i = 0; i < size; i++)
			ar.push_back(0.0);
	}

	~MyVector(){};

	double& operator[](int j) 
    {
        return ar[j];
    }

	
	MyVector operator+(MyVector a) 
	{
		MyVector new_vector = MyVector(ar.size());
		if(a.ar.size() == ar.size())
		for(int i = 0; i < ar.size(); i++)
			 new_vector.ar[i] = ar[i] + a.ar[i];
		return new_vector;
	}

	MyVector operator-(MyVector a) 
	{
		MyVector new_vector = MyVector(ar.size());
		assert(a.ar.size() == ar.size());
		for(int i = 0; i < ar.size(); i++)
			 new_vector.ar[i] = ar[i] - a.ar[i];
		return new_vector;
	}

	MyVector& operator=(MyVector vec)
	{
		assert(ar.size() == vec.ar.size());
		ar = vec.ar;
		return *this;
	}

	MyVector operator*(double a) 
	{
		MyVector new_vector = MyVector(ar.size());
		for(int i = 0; i < ar.size(); i++)
			 new_vector.ar[i] = ar[i] * a;
		return new_vector;
	}

	MyVector operator/(double a) 
	{
		MyVector new_vector = MyVector(ar.size());
		for(int i = 0; i < ar.size(); i++)
			 new_vector.ar[i] = ar[i] / a;
		return new_vector;
	}

	void initialize(int size)
	{
		ar.reserve(size);
		for(int i = 0; i < size; i++)
			ar.push_back(0.0);
	}

	void make_zero()
	{
		for(int i = 0; i < ar.size(); i++)
			ar[i] = 0.0;
	}

	double norm()
	{
		double sum = 0;
		for(int i = 0; i < ar.size(); i++)
			sum += ar[i] * ar[i];

		return sqrt(sum);
	}

	void output(FILE *f_out)
	{
		for(int i = 0; i < ar.size(); i++)
			fprintf(f_out, "%.20lf\n", ar[i]);
	}

};

double scal(MyVector v1, MyVector v2)
{
	double sum = 0;
	if(v1.ar.size() == v2.ar.size())
	for(int i = 0; i < v1.ar.size(); i++)
		sum += v1[i] * v2[i];
	return sum;
}

struct Matrix
{
	//����������� �������� ������
	int n;//����������� �������
	int size;//����������� ��������,��� �������� �������������� �������� 

	vector <double> ggl;//������ � ����������������� ��������������� ����������
	vector <double> ggu;//������ � ������������������ ��������������� ����������
	vector <double> di;//���������
	MyVector b;//������ ������ �����
	vector <int> ig;//��������� ������ �����(��������)
	vector <int> jg;//������ ��������(�����) ��������������� ���������

	
	Matrix(){};

	Matrix(int size1, int size2)
	{
		initialize(size1, size2);
	}

	void initialize(int size1, int size2);
	void reinitialize();

	Matrix& operator=(Matrix matrix)
	{
		n = matrix.n;
		size = matrix.size;
		ggl = matrix.ggl;
		ggu = matrix.ggu;
		di = matrix.di;
		b = matrix.b;
		ig = matrix.ig;
		jg = matrix.jg;
		return *this;
	}

	//��������� �� ������
	MyVector operator*(MyVector a) 
	{
		int i, j, k, kol;
		int iend;
		MyVector new_vector = MyVector(a.ar.size());

		assert(a.ar.size() == n);
		for(i = 0; i < n; i++)
		{
			kol = ig[i+1] - ig[i];//���������� ��������� ��������� ������ (�������) �� �������
								  //���������� �������� �� ������������� �������� (�� ������� ���)
			iend = ig[i+1];
			k = ig[i]; // ����� ������� �������� �������� ������ (�������) 

			new_vector[i] = di[i] * a[i];//�� ������� ���������

			for(; k < iend; k++)//�������� �� ���� ��������� i ������ (�������)
			{
				j = jg[k];
				new_vector[i] += ggl[k] * a[j];//�� ������� ������������
				new_vector[j] += ggu[k] * a[i];//�� �������� ������������
			}
		}

		return new_vector;
	}
	MyVector operator/(MyVector a) 
	{
		int i, j, k, kol;
		int iend;
		MyVector new_vector = MyVector(a.ar.size());

		assert(a.ar.size() == n);
		for(i = 0; i < n; i++)
		{
			kol = ig[i+1] - ig[i];//���������� ��������� ��������� ������ (�������) �� �������
								  //���������� �������� �� ������������� �������� (�� ������� ���)
			iend = ig[i+1];
			k = ig[i]; // ����� ������� �������� �������� ������ (�������) 

			new_vector[i] = di[i] * a[i];//�� ������� ���������

			for(; k < iend; k++)//�������� �� ���� ��������� i ������ (�������)
			{
				j = jg[k];
				new_vector[i] += ggu[k] * a[j];//�� ������� ������������
				new_vector[j] += ggl[k] * a[i];//�� �������� ������������
			}
		}

		return new_vector;
	}

	~Matrix(){};

	MyVector Uv(MyVector v);

};

struct DenseMatrix
{
	int n_lines, n_columns;
	vector <MyVector> ar;
	
	DenseMatrix(){};

	DenseMatrix(int size1, int size2)
	{
		n_lines = size1; n_columns = size2;

		ar.reserve(n_columns);
		for(int i = 0; i < n_columns; i++)
			ar.push_back(MyVector (n_lines));
	}

	DenseMatrix& operator=(DenseMatrix matrix)
	{
		n_lines = matrix.n_lines, n_columns = matrix.n_columns;
		for(int i = 0; i < n_columns; i++)
			ar[i] = matrix.ar[i];
		
		return *this;
	}

	MyVector& operator[](int j) 
    {
        return ar[j];
    }

	//��������� �� ������
	MyVector operator*(MyVector a) 
	{
		MyVector new_vector = n_lines;

		assert(a.ar.size() == n_columns);
		for(int j = 0; j < n_columns; j++)
			for(int i = 0; i < n_lines; i++)
				new_vector[i] += ar[j][i] * a[j];

		return new_vector;
	}

	~DenseMatrix(){};

	void initialize(int size1, int size2)
	{
		n_lines = size1; n_columns = size2;

		ar.reserve(n_columns);
		for(int i = 0; i < n_columns; i++)
			ar.push_back(MyVector (n_lines));
	}
};

struct Logger
{
	FILE *log_f;
	void send_current_information(double r_norm, int iteration_number)
	{
		fprintf(log_f, "%d     %.20lf\n", iteration_number, r_norm);
	};
};

struct BoundaryCondition1
{
	int elem;
	int edges[4]; //�����,������,������, �������: 1 - ����, 0 - ���
	int formula_number;
};

struct BoundaryCondition2
{
	int elem;
	int edges[4]; //�����,������,������, �������: 1 - ����, 0 - ���
	int formula_number;
};

struct BoundaryCondition3
{
	int elem;
	int edges[4]; //�����,������,������, �������: 1 - ����, 0 - ���
	int formula_number;
};

struct SLAE
{
	int n; //����������� ����
	int m; //������� ������ gmres
	int max_iter; //max ���������� ��������
	int max_iter_nonlinear;
	double eps; //�������� ������� ����
	double sigma, mu1, mu2; //������������ ������������
	Matrix A; //�������
	Partition P; //��������� �������
	Logger logger; //������ ��� ������ ���������� � �������� ������� ����
	vector <BoundaryCondition1> boundaries1; //������ ������� �������
	vector <BoundaryCondition2> boundaries2; //������ -//-
	vector <BoundaryCondition3> boundaries3; //������ -//-

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

	vector <double> LU_ggu; //����������������� �������������� �������� U
	vector <double> LU_ggl; //���������������� �������������� �������� L
	vector <double> LU_di; //������������ �������� L
	vector <double> LU_ggu2; //����������������� �������������� �������� U ��� LU-��������
	vector <double> LU_ggl2; //���������������� �������������� �������� L ��� LU-��������
	vector <double> LU_di2; //������������ �������� L ��� LU-��������
	vector <int> LU_ig2; //������� �������� ��� LU-��������
	int size_prof; //������ �������� ��������� ��� ����������� ������� LU-��������
	MyVector yl; //������� ������� Lyl=F
	MyVector yu; //������� ������� Uyu=F
	
	double gauss_points[2][9];//����� ������
	double gauss_weights[9];// ���� ������
	double gauss_points_1[3];//����� ������
	double gauss_weights_1[3];// ���� ������
    
	SLAE(){};

	SLAE(int max_number_of_iterations, 
		 int max_number_of_iterations_non_lin,
		 double epsilon, 
		 int gmres_m, 
		 FILE *grid_f_in, 
		 FILE *elements_f_in, 
		 FILE *log_f, 
		 FILE *boundary1)
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
					FILE *grid_f_in,
					FILE *elements_f_in,
					FILE *log_f,
					FILE *boundary1);

	void reinitialize();

	double get_hx(int element_number);
	double get_hy(int element_number);

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


	//������� �������
	void input_boundaries1(FILE* f_in);
	void input_boundaries2(FILE* f_in);
	void input_boundaries3(FILE* f_in);

	void calculate_all_boundaries1();
	void calculate_boundaries1(int number);

	void calculate_boundaries1_left(int number);
	void calculate_boundaries1_right(int number);
	void calculate_boundaries1_low(int number);
	void calculate_boundaries1_up(int number);

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

	void si_print(FILE *log_f, int iteration_number, double &normL2u, double &normL2p);
	double find_relaxation_parameter(MyVector q_current, MyVector q_previous, double &residual_previous);
	void simple_iterations();

	double SLAE::diff_normL2_p(MyVector q_solution);//����������� ������� � ����� L2
	double SLAE::diff_normL2_u(MyVector q_solution);//����������� ������� � ����� L2


	//�������������������
	void LU();
	void LYF(MyVector b);
	void UXY(MyVector b);
	void LYFt(MyVector b);
	void UXYt(MyVector b);
	MyVector Uv(MyVector v);

	//������ �������
	void run(FILE *solution_f_out, FILE *info_f_out);

	void output(FILE *solution_f_out, FILE *info_f_out, double normL2u, double normL2p)
	{		
		FILE *res = fopen("result.txt", "w");
		int n_nodes = P.nodes.size();
		for(int i = 0; i < n_nodes; i++)
		{
			if(abs(P.nodes[i].x - 0.5) < 1e-10)
			fprintf(res, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", 
					P.nodes[i].x, 
					P.nodes[i].y, 
					Ux_numerical[i], 
					Uy_numerical[i], 
					P_numerical[i],
					calculate_p_analytic(0,P.nodes[i].x,P.nodes[i].y));
		}
		fclose(res);
		Ux_numerical.output(solution_f_out);
		fprintf(solution_f_out, "\n\n\n");
		Uy_numerical.output(solution_f_out);
		fprintf(solution_f_out, "\n\n\n");
		P_numerical.output(solution_f_out);
		fprintf(solution_f_out, "\n\n\n");
		fprintf(info_f_out, "norm L2 u:|u*-u|=%.4e\nnorm L2 p:|p*-p|=%.4e\neps=%.2e\n", 
				normL2u, normL2p, eps);
	};
};


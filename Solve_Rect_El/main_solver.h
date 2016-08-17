#pragma once
#include <string>
#include <fstream>

#include "boundary_conditions.h"
#include "myvector.h"
#include "boundaries.h"
#include "logger.h"
#include "solver.h"
#include "SLAE.h"

namespace mainsolver
{
	class MainSolver : public boundary_conditions::BoundaryConditionsSupport, public boundaries::InternalBoundaries, public boundaries::OuterBoundaries
	{
	public:
		slae::SLAE my_slae;
		logger::Logger logger; //логгер для вывода информации о процессе решения СЛАУ

		//локальные матрицы
		double G[4][4], P1[4][4], P2[4][4], C[4][4]; 
		double F[4]; //локальный вектор правой части

		myvector::MyVector Ux_numerical; //численное решение Ux
		myvector::MyVector Uy_numerical; //численное решение Uy
		myvector::MyVector P_numerical; //численное решение P

		myvector::MyVector q_prev; //вектор весов с предыдущей итерации по нелинейности
    
		MainSolver(){};

		MainSolver(std::ifstream& grid_f_in, 
				   std::ifstream& elements_f_in,
				   std::string log_f,
				   std::ifstream& boundary1,
				   std::ifstream& boundary2,
				   std::ifstream& boundary3)
		{
			initialize(grid_f_in,
					   elements_f_in,
					   log_f,
					   boundary1,
					   boundary2,
					   boundary3);
		}

		~MainSolver(){};

		void initialize(std::ifstream& grid_f_in,
						std::ifstream& elements_f_in,
						std::string log_f,
						std::ifstream& boundary1,
						std::ifstream& boundary2,
						std::ifstream& boundary3);

		void reinitialize();

		void build_slae(myvector::MyVector q_calc);

		double find_relaxation_parameter(myvector::MyVector q_current, 
										 myvector::MyVector q_previous, 
										 double &residual_previous);
		void simple_iterations(solver::Solver& s);
		void linear(std::ofstream& solution_f_out,
					std::ofstream& info_f_out,
					solver::Solver& s);
		void solve();
		void solve(std::ofstream& solution_f_out, std::ofstream& info_f_out);
							    
		double get_solution_in_point_ux(double x, double y, int element_number, myvector::MyVector qi);
		double get_solution_in_point_uy(double x, double y, int element_number, myvector::MyVector qi);
		double get_solution_in_point_uxdx(double x, double y, int element_number, myvector::MyVector qi);
		double get_solution_in_point_uydy(double x, double y, int element_number, myvector::MyVector qi);
		double get_solution_in_point_p(double x, double y, int element_number, myvector::MyVector qi);

		double get_solution_in_point2_ux(double x, double y, myvector::MyVector qi);
		double get_solution_in_point2_uy(double x, double y, myvector::MyVector qi);
		double get_solution_in_point2_p(double x, double y, myvector::MyVector qi);

		void get_vector_solution_in_nodes_ux(myvector::MyVector qi, myvector::MyVector &solution);
		void get_vector_solution_in_nodes_uy(myvector::MyVector qi, myvector::MyVector &solution);
		void get_vector_solution_in_nodes_p(myvector::MyVector qi, myvector::MyVector &solution);

		//локальные матрицы и векторы
		void calculate_locals(int element_number, myvector::MyVector q_calc);
		void calculate_G(int element_number);
		void calculate_C(int element_number, myvector::MyVector q_calc);
		void calculate_P1(int element_number);
		void calculate_P2(int element_number);
		void calculate_F(int element_number);

		double diff_normL2_p(myvector::MyVector q_solution);//погрешность решения в норме L2
		double diff_normL2_u(myvector::MyVector q_solution);//погрешность решения в норме L2
	};
}
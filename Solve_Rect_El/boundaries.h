#pragma once
#include "partition.h"
#include "parameters.h"
#include "integration.h"
#include "basis.h"
#include "matrix.h"

namespace boundaries
{
	//численные потоки по внутренним границам
	class InternalBoundaries : virtual public partition::Partition,
							   virtual public parameters::Parameters,
							   virtual public integration::Gauss_integration,
							   virtual public basis::Basis
	{
		//S для скорости
		//E для скорости
		//P_1 //давление по j
		//P_2 //давление по i
		//SP для давления

		double sigma, mu2; //коэффициенты стабилизации

		void calculate_ES_horizontal(int element_number1, int element_number2, matrix::Matrix& A);
		void calculate_ES_vertical(int element_number1, int element_number2, matrix::Matrix& A);

		void calculate_P_1_horizontal(int element_number1, int element_number2, matrix::Matrix& A);
		void calculate_P_1_vertical(int element_number1, int element_number2, matrix::Matrix& A);

		void calculate_P_2_horizontal(int element_number1, int element_number2, matrix::Matrix& A);
		void calculate_P_2_vertical(int element_number1, int element_number2, matrix::Matrix& A);

		void calculate_SP_horizontal(int element_number1, int element_number2, matrix::Matrix& A);
		void calculate_SP_vertical(int element_number1, int element_number2, matrix::Matrix& A);

		void add_ES_to_global(int element_number, int neighbor_element_number, matrix::Matrix& A, std::vector <std::vector<double>>& ES);
		void add_P_1_to_global(int element_number, int neighbor_element_number, matrix::Matrix& A, std::vector <std::vector<double>>& P_1);
		void add_P_2_to_global(int element_number, int neighbor_element_number, matrix::Matrix& A, std::vector <std::vector<double>>& P_2);
		void add_SP_to_global(int element_number, int neighbor_element_number, matrix::Matrix& A, std::vector <std::vector<double>>& SP);

	protected:

		void initialize_penalty_parameters();
		void calculate_internal_boundaries(int element_number, matrix::Matrix& A);
	};

	//численные потоки по внешним границам
	class OuterBoundaries : virtual public partition::Partition,
							virtual public parameters::Parameters,
							virtual public integration::Gauss_integration,
							virtual public basis::Basis
	{

		double sigma, mu2; //коэффициенты стабилизации	
		double E_out[4][4], P_1_out[4][4], P_2_out[4][4], SP_out[4][4]; //локальные матрицы

		void calculate_ES_out_left(int element_number, matrix::Matrix& A);
		void calculate_ES_out_right(int element_number, matrix::Matrix& A);
		void calculate_ES_out_low(int element_number, matrix::Matrix& A);
		void calculate_ES_out_up(int element_number, matrix::Matrix& A);

		void calculate_P_1_out_left(int element_number, matrix::Matrix& A);
		void calculate_P_1_out_right(int element_number, matrix::Matrix& A);
		void calculate_P_1_out_low(int element_number, matrix::Matrix& A);
		void calculate_P_1_out_up(int element_number, matrix::Matrix& A);

		void calculate_P_2_out_left(int element_number, matrix::Matrix& A);
		void calculate_P_2_out_right(int element_number, matrix::Matrix& A);
		void calculate_P_2_out_low(int element_number, matrix::Matrix& A);
		void calculate_P_2_out_up(int element_number, matrix::Matrix& A);

		void calculate_SP_out_left(int element_number, matrix::Matrix& A);
		void calculate_SP_out_right(int element_number, matrix::Matrix& A);
		void calculate_SP_out_low(int element_number, matrix::Matrix& A);
		void calculate_SP_out_up(int element_number, matrix::Matrix& A);

	protected:

		void initialize_penalty_parameters();
		void calculate_outer_boundaries(int element_number, matrix::Matrix& A);
	};

}
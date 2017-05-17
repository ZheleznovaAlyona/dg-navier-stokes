#pragma once
#include "partition.h"
#include "parameters.h"
#include "integration.h"
#include "basis.h"
#include "matrix.h"
#include "EdgeLocation.h"
#include <string>

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

		double gamma, sigma; //коэффициенты стабилизации

		void calculate_ES(int element_number1, int element_number2, matrix::Matrix& A, EdgeOrient orient);
		void calculate_P_1(int element_number1, int element_number2, matrix::Matrix& A, EdgeOrient orient);
		void calculate_P_2(int element_number1, int element_number2, matrix::Matrix& A, EdgeOrient orient);
		void calculate_SP(int element_number1, int element_number2, matrix::Matrix& A, EdgeOrient orient);

		void add_ES_to_global(int element_number, int neighbor_element_number, matrix::Matrix& A, std::vector <std::vector<double>>& ES);
		void add_P_1_to_global(int element_number, int neighbor_element_number, matrix::Matrix& A, std::vector <std::vector<double>>& P_1);
		void add_P_2_to_global(int element_number, int neighbor_element_number, matrix::Matrix& A, std::vector <std::vector<double>>& P_2);
		void add_SP_to_global(int element_number, int neighbor_element_number, matrix::Matrix& A, std::vector <std::vector<double>>& SP);

	protected:

		void initialize_penalty_parameters(std::string fileName);
		void calculate_internal_boundaries(int element_number, matrix::Matrix& A);
	};

	//численные потоки по внешним границам
	class OuterBoundaries : virtual public partition::Partition,
							virtual public parameters::Parameters,
							virtual public integration::Gauss_integration,
							virtual public basis::Basis
	{

		double gamma, sigma; //коэффициенты стабилизации	

		void calculate_ES_out(int element_number, matrix::Matrix& A, EdgeSide edgeSide);
		void calculate_P_1_out(int element_number, matrix::Matrix& A, EdgeSide edgeSide);
		void calculate_P_2_out(int element_number, matrix::Matrix& A, EdgeSide edgeSide);
		void calculate_SP_out(int element_number, matrix::Matrix& A, EdgeSide edgeSide);

	protected:

		void initialize_penalty_parameters(std::string fileName);
		void calculate_outer_boundaries(int element_number, matrix::Matrix& A);
	};

}
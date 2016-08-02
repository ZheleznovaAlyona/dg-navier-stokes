#pragma once
#include "partition.h"
#include "parameters.h"
#include "integration.h"
#include "basis.h"
#include "matrix.h"

namespace boundaries
{
	//��������� ������ �� ���������� ��������
	class InternalBoundaries : virtual public partition::Partition,
							   virtual public parameters::Parameters,
							   virtual public integration::Gauss_integration,
							   virtual public basis::Basis
	{
	public:
		double sigma, mu2; //������������ ������������
		double E[8][8], P_1[8][8], P_2[8][8], SP[8][8]; //��������� �������

		void initialize_penalty_parameters();

		void calculate_internal_boundaries(int element_number, matrix::Matrix& A);

		void calculate_ES_horizontal(int element_number1, int element_number2);
		void calculate_ES_vertical(int element_number1, int element_number2);

		void calculate_P_1_horizontal(int element_number1, int element_number2);
		void calculate_P_1_vertical(int element_number1, int element_number2);

		void calculate_P_2_horizontal(int element_number1, int element_number2);
		void calculate_P_2_vertical(int element_number1, int element_number2);

		void calculate_SP_horizontal(int element_number1, int element_number2);
		void calculate_SP_vertical(int element_number1, int element_number2);

		void add_ES_to_global(int element_number, int neighbor_element_number, matrix::Matrix& A);
		void add_P_1_to_global(int element_number, int neighbor_element_number, matrix::Matrix& A);
		void add_P_2_to_global(int element_number, int neighbor_element_number, matrix::Matrix& A);
		void add_SP_to_global(int element_number, int neighbor_element_number, matrix::Matrix& A);
	};

	//��������� ������ �� ������� ��������
	class OuterBoundaries : virtual public partition::Partition,
							virtual public parameters::Parameters,
							virtual public integration::Gauss_integration,
							virtual public basis::Basis
	{
	public:
		double sigma, mu2; //������������ ������������	
		double E_out[4][4], P_1_out[4][4], P_2_out[4][4], SP_out[4][4]; //��������� �������

		void initialize_penalty_parameters();

		void calculate_outer_boundaries(int element_number, matrix::Matrix& A);

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
	};

}
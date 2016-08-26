#pragma once
#include <fstream>
#include <string>
#include "myvector.h"
#include "partition.h"
#include "parameters.h"

namespace logger
{
	class Logger
	{
		std::ofstream log_f;

		enum 
		{
			TO_FILE,
			TO_STDO
		} p;

	public:
		void open(std::string file_name);
		Logger();
		~Logger();

		void send_current_information(double r_norm, int iteration_number);
		void send_current_information_to_screen(double r_norm, int iteration_number);
		void send_current_information_to_screen_si(int si_iteration_number);
		void send_message_build_slae();
		void send_message_create_portret();
		void send_message_solution();
		void send_message_Ux();
		void send_message_Uy();
		void send_message_P();
		void send_message_norms();
		void send_inf_UL2_norm(double normL2u);
		void send_inf_PL2_norm(double normL2p);

		void si_print(int iteration_number,
					  double &normL2u,
					  double &normL2p,
					  myvector::MyVector& Ux_num,
					  myvector::MyVector& Uy_num,
					  myvector::MyVector& P_num,
					  partition::Partition& P);

		void output(std::ofstream& solution_f_out, 
					std::ofstream& info_f_out, 
					double normL2u, 
					double normL2p, 
					myvector::MyVector& Ux_num,
					myvector::MyVector& Uy_num,
					myvector::MyVector& P_num,
					partition::Partition& P);

	};
}
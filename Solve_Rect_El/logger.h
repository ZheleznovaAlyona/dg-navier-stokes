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
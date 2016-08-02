#pragma once
#include <fstream>
#include <string>

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
		void send_current_information(double r_norm, int iteration_number);
		void open(std::string file_name);
		Logger();
		~Logger();
	};
}
#include "logger.h"
#include <iomanip>
#include <iostream>

using namespace std;

namespace logger
{
	void Logger::send_current_information(double r_norm, int iteration_number)
	{
		ostream* os;

		if(p == TO_FILE)
			os = &log_f;
		else
			os = &cout;

		*os << iteration_number << "     " << setprecision(20) << r_norm << endl;
	}

	Logger::Logger()
	{
		p = TO_STDO;
	}

	void Logger::open(string file_name)
	{
		log_f.open(file_name);
		p = TO_FILE;
	}

	Logger::~Logger()
	{
		if(p == TO_FILE)
			log_f.close();
	}
}
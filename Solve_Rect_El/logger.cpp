#include "logger.h"
#include <iomanip>
#include <iostream>

using namespace std;
using namespace myvector;
using namespace partition;
using namespace parameters;

namespace logger
{
	void Logger::send_current_information(double r_norm, int iteration_number)
	{
		ostream* os;

		if(p == TO_FILE)
			os = &log_f;
		else
			os = &cout;

		*os << iteration_number << "     " << scientific << setprecision(20) << r_norm << endl;
	}

	void Logger::send_current_information_to_screen_si(int si_iteration_number)
	{
		cout << endl << si_iteration_number << " iteration is in process" << endl;
	}

	void Logger::send_current_information_to_screen(double r_norm, int iteration_number)
	{
		cout << iteration_number << "\tr=" << scientific << setprecision(10) << r_norm << endl;
	}

	void Logger::send_message_build_slae()
	{
		cout << "Building SLAE..." << endl;
	}

	void Logger::send_message_create_portret()
	{
		cout << "Creating matrix-portret..." << endl;
	}

	void Logger::send_message_solution()
	{
		cout << "Getting solution..." << endl;
	}

	void Logger::send_message_Ux()
	{
		cout << "Getting Ux..." << endl;
	}

	void Logger::send_message_Uy()
	{
		cout << "Getting Uy..." << endl;
	}

	void Logger::send_message_P()
	{
		cout << "Getting P..." << endl;
	}

	void Logger::send_message_norms()
	{
		cout << "Calculating L2-norms..." << endl;
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

	void Logger::send_inf_UL2_norm(double normL2u)
	{
		cout << "|U|(L2) = " << scientific << setprecision(4) << normL2u << endl;
	}

	void Logger::send_inf_PL2_norm(double normL2p)
	{
		cout << "|P|(L2) = " << scientific << setprecision(4) << normL2p << endl;
	}

	void Logger::si_print(int iteration_number,
					double &normL2u,
					double &normL2p,
					MyVector& Ux_num,
					MyVector& Uy_num,
					MyVector& P_num,
					Partition& P)
	{
		string f_name_s, f_name_i;

		log_f << "---" << iteration_number << "---" << endl;
		f_name_s = string("s_") + to_string(iteration_number) + ".txt";
		f_name_i = string("i_") + to_string(iteration_number) + ".txt";
		ofstream solution_f_out(f_name_s), info_f_out(f_name_i);

		output(solution_f_out, info_f_out, normL2u, normL2p, Ux_num, Uy_num, P_num, P);
		solution_f_out.close();
		info_f_out.close();
	}

	void Logger::output(ofstream& solution_f_out, 
				ofstream& info_f_out, 
				double normL2u, 
				double normL2p, 
				MyVector& Ux_num,
				MyVector& Uy_num,
				MyVector& P_num,
				Partition& P)
	{		
		ofstream res("result.txt");

		int n_nodes = P.nodes.size();
		for(int i = 0; i < n_nodes; i++)
		{
			if(abs(P.nodes[i].x - 0.5) < 1e-10)
				res << P.nodes[i].x << "\t" << P.nodes[i].y << "\t" << Ux_num[i]
				<< "\t" << Uy_num[i] << "\t" << P_num[i] << "\t"
				<< Parameters::calculate_p_analytic(0, P.nodes[i].x, P.nodes[i].y) << endl;
		}

		res.close();

		solution_f_out << Ux_num << endl << endl << endl;
		solution_f_out << Uy_num << endl << endl << endl;
		solution_f_out << P_num << endl << endl << endl;

		info_f_out << "norm L2 u:|u*-u|=" << scientific << setprecision(4) << normL2u 
			<< endl << "norm L2 p:|p*-p|=" << scientific << setprecision(4) << normL2p
			<< endl;
	}
}
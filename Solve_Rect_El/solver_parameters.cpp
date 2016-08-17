#include "solver_parameters.h"
#include "rapidjson/document.h"

using namespace std;
using namespace rapidjson;


namespace solverparameters
{
	void SolverParameters::initialize(std::string file_name)
	{
		ifstream file_in(file_name);
		string json;

		while (!file_in.eof())
		{
			string line;
			file_in >> line;
			json += line;
		}

		file_in.close();
		
		Document d_in;
		d_in.Parse(json.c_str());

		auto& j_solver = d_in["solver"];
		this->epsilon = j_solver["epsilon"].GetDouble();
		this->gmres_m = j_solver["gmres_m"].GetInt();
		this->max_number_of_iterations = j_solver["max_number_of_iterations"].GetInt();
		this->max_number_of_iterations_non_lin = j_solver["max_number_of_iterations_non_linear"].GetInt();
	}
}
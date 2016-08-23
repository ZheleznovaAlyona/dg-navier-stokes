#include "testing_parameters.h"
#include "rapidjson/document.h"
#include <fstream>

using namespace std;
using namespace rapidjson;

namespace testingparameters
{
	void Testing_parameters::initialize(std::string file_name)
	{
		ifstream file_in(file_name);
		string json;
		string solver_name;

		while (!file_in.eof())
		{
			string line;
			file_in >> line;
			json += line;
		}

		file_in.close();

		Document d_in;
		d_in.Parse(json.c_str());

		auto& j_testp = d_in["testp"];
		this->use_LU = j_testp["useLU"].GetBool();
		this->test = j_testp["test"].GetInt();
		solver_name = j_testp["solver"].GetString();

		if (solver_name == "BiCGStab") this->solver = 1;
		else
			if (solver_name == "GMRES") this->solver = 2;
			else
				if(solver_name == "BCGandGMRES") this->solver = 3;
				else this->solver = 2;
	}

}
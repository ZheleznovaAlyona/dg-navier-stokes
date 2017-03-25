#include "testing_parameters.h"
#include "rapidjson/document.h"
#include <fstream>
#include <map>

using namespace std;
using namespace rapidjson;

namespace testingparameters
{
	bool Testing_parameters::use_LU = false;
	int Testing_parameters::test = 0;
	int Testing_parameters::solver = 0;

	Testing_parameters::Testing_parameters() {}
	Testing_parameters::~Testing_parameters() {}
	Testing_parameters::Testing_parameters(Testing_parameters& parameters) {}

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
		use_LU = j_testp["useLU"].GetBool();
		test = j_testp["test"].GetInt();
		solver_name = j_testp["solver"].GetString();

		map<string, int> _map;

		_map["BiCGStab"] = 1;
		_map["GMRES"] = 2;
		_map["BCGandGMRES"] = 3;
		_map["BCG"] = 4;

		if (_map.find(solver_name) == _map.end())
			solver = 2;
		else
			solver = _map[solver_name];
	}

}
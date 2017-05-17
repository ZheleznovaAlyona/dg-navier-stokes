#include <fstream>
#include "main_solver.h"
#include "testing_parameters.h"

using namespace std;
using namespace mainsolver;
using namespace solver;
using namespace testingparameters;

void main()
{
	ifstream l1_in("l1.dt"), l2_in("l2.dt"), l3_in("l3.dt"), grid_in("grid.dt"), elements_in("elements.dt");
	Testing_parameters::initialize("testing_parameters.json");

	MainSolver problem_solver(grid_in, elements_in, "log.txt", l1_in, l2_in, l3_in, "PenaltyParameters.json");	
	//problem_solver.solve();
	ofstream solution_out("solution.txt"), info_out("info.txt");
	problem_solver.solve(solution_out, info_out);
	solution_out.close();
	info_out.close();

	l1_in.close();
	l2_in.close();
	l3_in.close();
	grid_in.close();
	elements_in.close();

	system("pause");
}
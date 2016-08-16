#include <fstream>
#include "main_solver.h"

using namespace std;
using namespace mainsolver;
using namespace solver;

void main()
{
	ifstream l1_in("l1.txt"), l2_in("l2.txt"), l3_in("l3txt"), grid_in("grid.txt"), elements_in("elements.txt");
	ofstream solution_out("solution.txt"), info_out("info.txt");

	MainSolver problem_solver(grid_in, elements_in, "log.txt", l1_in, l2_in, l3_in);	
	problem_solver.solve();

	l1_in.close();
	grid_in.close();
	elements_in.close();
	solution_out.close();
	info_out.close();

	system("pause");
}
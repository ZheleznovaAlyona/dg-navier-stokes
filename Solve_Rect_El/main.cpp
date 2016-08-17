#include <fstream>
#include "main_solver.h"

using namespace std;
using namespace mainsolver;
using namespace solver;

void main()
{
	ifstream l1_in("l1.txt"), l2_in("l2.txt"), l3_in("l3.txt"), grid_in("grid.txt"), elements_in("elements.txt");

	MainSolver problem_solver(grid_in, elements_in, "log.txt", l1_in, l2_in, l3_in);	
	problem_solver.solve();

	l1_in.close();
	l2_in.close();
	l3_in.close();
	grid_in.close();
	elements_in.close();

	system("pause");
}
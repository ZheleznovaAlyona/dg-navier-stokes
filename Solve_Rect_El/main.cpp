#include <fstream>
#include "main_solver.h"

using namespace std;
using namespace mainsolver;

void main()
{
	ifstream l1_in("l1.txt"), grid_in("grid.txt"), elements_in("elements.txt");
	ofstream solution_out("solution.txt"), info_out("info.txt");

	SLAE my_SLAE(grid_in, elements_in, l1_in);
	//my_SLAE.run(solution_f_out, info_f_out);
	my_SLAE.simple_iterations();

	l1_in.close();
	grid_in.close();
	elements_in.close();
	solution_out.close();
	info_out.close();

	system("pause");
}
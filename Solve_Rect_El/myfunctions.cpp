#include "myfunctions.h"

using namespace std;
using namespace myvector;

double pow_i(int i, double a)
{
	double res = 1;
	for(int j = 1; j <= i; j++)
		res *= a;
	return res;
}

void initialize_vector(vector <double> &v, int size)
{
	v.resize(size);
	memset(&v[0], 0, size * sizeof(double)); //обнуляем
}

void initialize_vector(vector <int> &v, int size)
{
	v.resize(size);
	memset(&v[0], 0, size * sizeof(int)); //обнуляем
}

double scal(MyVector v1, MyVector v2)
{
	double sum = 0;
	if(v1.ar.size() == v2.ar.size())
	for(unsigned int i = 0; i < v1.ar.size(); i++)
		sum += v1[i] * v2[i];
	return sum;
}
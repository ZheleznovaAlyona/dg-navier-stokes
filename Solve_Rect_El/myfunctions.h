#pragma once
#include <vector>
#include "myvector.h"

double pow_i(int i, double a);
void initialize_vector(std::vector <double> &v, int size);
void initialize_vector(std::vector <int> &v, int size);
double scal(myvector::MyVector v1, myvector::MyVector v2);
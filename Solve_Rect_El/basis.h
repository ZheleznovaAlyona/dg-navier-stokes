#pragma once
#include <functional>
#include "partition.h"

namespace basis
{
	class Basis
	{
	public:
		int n_func_p; //число базисных функций p
		int n_func_u; //число базисных функций u
		std::function<double(double, double)> psi[4]; //указатели на функции вычисления базисных функций p в точке
		std::function<double(double, double)> dpsiksi[4]; //указатели на функции вычисления d/dksi базисных функций p в точке
		std::function<double(double, double)> dpsietta[4]; //указатели на функции вычисления d/detta базисных функций p в точке
		std::function<double(double, double)> phix[4]; //указатели на функции вычисления базисных функций ux в точке
		std::function<double(double, double)> phiy[4]; //указатели на функции вычисления базисных функций uy в точке
		std::function<double(double, double)> dphixksi[4]; //указатели на функции вычисления d/dksi базисных функций ux в точке
		std::function<double(double, double)> dphixetta[4]; //указатели на функции вычисления d/detta базисных функций ux в точке
		std::function<double(double, double)> dphiyksi[4]; //указатели на функции вычисления d/dksi базисных функций uy в точке
		std::function<double(double, double)> dphiyetta[4]; //указатели на функции вычисления d/detta базисных функций uy в точке

		void initialize();

		double phix_i(int i, double x, double y, int element_number, partition::Partition& p);
		double phiy_i(int i, double x, double y, int element_number, partition::Partition& p);
		double phixdx_i(int i, double x, double y, int element_number, partition::Partition& p);
		double phiydy_i(int i, double x, double y, int element_number, partition::Partition& p);
		double psi_i(int i, double x, double y, int element_number, partition::Partition& p);
	};
}
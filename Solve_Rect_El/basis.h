#pragma once
#include <functional>
#include "partition.h"

namespace basis
{
	const int n_func_p = 4; //число базисных функций p
	const int n_func_u = 4; //число базисных функций u

	class Basis
	{

	protected:

		std::function<double(double, double)> psi[n_func_p]; //указатели на функции вычисления базисных функций p в точке
		std::function<double(double, double)> dpsiksi[n_func_p]; //указатели на функции вычисления d/dksi базисных функций p в точке
		std::function<double(double, double)> dpsietta[n_func_p]; //указатели на функции вычисления d/detta базисных функций p в точке
		std::function<double(double, double)> phix[n_func_u]; //указатели на функции вычисления базисных функций ux в точке
		std::function<double(double, double)> phiy[n_func_u]; //указатели на функции вычисления базисных функций uy в точке
		std::function<double(double, double)> dphixksi[n_func_u]; //указатели на функции вычисления d/dksi базисных функций ux в точке
		std::function<double(double, double)> dphixetta[n_func_u]; //указатели на функции вычисления d/detta базисных функций ux в точке
		std::function<double(double, double)> dphiyksi[n_func_u]; //указатели на функции вычисления d/dksi базисных функций uy в точке
		std::function<double(double, double)> dphiyetta[n_func_u]; //указатели на функции вычисления d/detta базисных функций uy в точке

		void initialize();

		double phix_i(int i, double x, double y, int element_number, partition::Partition& p);
		double phiy_i(int i, double x, double y, int element_number, partition::Partition& p);
		double phixdx_i(int i, double x, double y, int element_number, partition::Partition& p);
		double phiydy_i(int i, double x, double y, int element_number, partition::Partition& p);
		double psi_i(int i, double x, double y, int element_number, partition::Partition& p);
	};
}
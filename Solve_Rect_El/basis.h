#pragma once
#include <functional>
#include "partition.h"

namespace basis
{
	const int n_func_p = 4; //����� �������� ������� p
	const int n_func_u = 4; //����� �������� ������� u

	class Basis
	{

	protected:

		std::function<double(double, double)> psi[n_func_p]; //��������� �� ������� ���������� �������� ������� p � �����
		std::function<double(double, double)> dpsiksi[n_func_p]; //��������� �� ������� ���������� d/dksi �������� ������� p � �����
		std::function<double(double, double)> dpsietta[n_func_p]; //��������� �� ������� ���������� d/detta �������� ������� p � �����
		std::function<double(double, double)> phix[n_func_u]; //��������� �� ������� ���������� �������� ������� ux � �����
		std::function<double(double, double)> phiy[n_func_u]; //��������� �� ������� ���������� �������� ������� uy � �����
		std::function<double(double, double)> dphixksi[n_func_u]; //��������� �� ������� ���������� d/dksi �������� ������� ux � �����
		std::function<double(double, double)> dphixetta[n_func_u]; //��������� �� ������� ���������� d/detta �������� ������� ux � �����
		std::function<double(double, double)> dphiyksi[n_func_u]; //��������� �� ������� ���������� d/dksi �������� ������� uy � �����
		std::function<double(double, double)> dphiyetta[n_func_u]; //��������� �� ������� ���������� d/detta �������� ������� uy � �����

		void initialize();

		double phix_i(int i, double x, double y, int element_number, partition::Partition& p);
		double phiy_i(int i, double x, double y, int element_number, partition::Partition& p);
		double phixdx_i(int i, double x, double y, int element_number, partition::Partition& p);
		double phiydy_i(int i, double x, double y, int element_number, partition::Partition& p);
		double psi_i(int i, double x, double y, int element_number, partition::Partition& p);
	};
}
#pragma once
#include <functional>

namespace basis
{
	class Basis
	{
	public:
		int n_func; //����� �������� �������
		std::function<double(double, double)> psi[4]; //��������� �� ������� ���������� �������� ������� p � �����
		std::function<double(double, double)> dpsiksi[4]; //��������� �� ������� ���������� d/dksi �������� ������� p � �����
		std::function<double(double, double)> dpsietta[4]; //��������� �� ������� ���������� d/detta �������� ������� p � �����
		std::function<double(double, double)> phix[4]; //��������� �� ������� ���������� �������� ������� ux � �����
		std::function<double(double, double)> phiy[4]; //��������� �� ������� ���������� �������� ������� uy � �����
		std::function<double(double, double)> dphixksi[4]; //��������� �� ������� ���������� d/dksi �������� ������� ux � �����
		std::function<double(double, double)> dphixetta[4]; //��������� �� ������� ���������� d/detta �������� ������� ux � �����
		std::function<double(double, double)> dphiyksi[4]; //��������� �� ������� ���������� d/dksi �������� ������� uy � �����
		std::function<double(double, double)> dphiyetta[4]; //��������� �� ������� ���������� d/detta �������� ������� uy � �����

		void initialize();
	};
}
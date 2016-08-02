#pragma once
namespace parameters
{

	class Parameters
	{
	public:
		
		//������ ������� �������
		static double gx(int formula_number, double x, double y);
		static double gy(int formula_number, double x, double y);

		//������� ������ �����
		static double calculate_fx(int area_number, double x, double y);
		static double calculate_fy(int area_number, double x, double y);

		//������������� ������� ��� ��������
		static double calculate_ux_analytic(int area_number, double x, double y);
		static double calculate_uy_analytic(int area_number, double x, double y);

		//����������� �������������� ������� ��� ��������
		static double calculate_uxdx_analytic(int area_number, double x, double y);
		static double calculate_uydy_analytic(int area_number, double x, double y);

		//������������� ������� ��� ��������
		static double calculate_p_analytic(int area_number, double x, double y);

		//������������ ������
		static double calculate_lambda(int area_number);
		static double calculate_rho(int area_number);

	};

}	

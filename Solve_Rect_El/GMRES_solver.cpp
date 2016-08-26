#include "solver.h"
#include "myfunctions.h"
#include <iostream>
#include <iomanip>

using namespace myvector;
using namespace densematrix;
using namespace logger;
using namespace std;

namespace solver
{
	MyVector GMRES::solve(MyVector U_begin, double &normL2u, double &normL2p, slae::SLAE& slae_in, Logger& my_logger)
	{
		MyVector r(slae_in.n), x(slae_in.n), d(s_parameters.gmres_m + 1), z(s_parameters.gmres_m), w(slae_in.n), tmp(slae_in.n), f(slae_in.n), rr2(slae_in.n);
		DenseMatrix V(slae_in.n, s_parameters.gmres_m), H(s_parameters.gmres_m + 1, s_parameters.gmres_m);

		double norm_r, norm_f;
		bool continue_;

		x = U_begin;
		f = slae_in.b;
		r = f - slae_in.A * x;
		slae_in.A.LYF(r); r = slae_in.A.yl;
		norm_r = r.norm();
		slae_in.A.LYF(f); norm_f = slae_in.A.yl.norm();
		x = slae_in.A.Uv(U_begin);

		for(int k_iter = 1; k_iter <= s_parameters.max_number_of_iterations && norm_r / norm_f > s_parameters.epsilon; k_iter++)
		{
			d.make_zero();
			V[0] = r / norm_r;

			continue_ = true;
			for(int j = 1; j <= s_parameters.gmres_m && continue_; j++)
			{
				slae_in.A.UXY(V[j - 1]);
				tmp = slae_in.A * slae_in.A.yu;
				slae_in.A.LYF(tmp); w = slae_in.A.yl;

				for(int l = 1; l <= j; l++)
				{
					H[j - 1][l - 1] = scal(V[l - 1], w);
					w = w - V[l - 1] * H[j - 1][l - 1];
				}

				H[j - 1][j] = w.norm();
				if(abs(H[j - 1][j]) < 1E-14)
				{
					s_parameters.gmres_m = j;
					continue_ = false;
				}
				else
				{
					if(s_parameters.gmres_m != j)
						V[j] = w / H[j - 1][j];
				}
			}

			d[0] = norm_r;
			solve_min_sqr_problem(d, H, z);
			x = x + V * z;
			slae_in.A.UXY(x);
			tmp = f - slae_in.A * slae_in.A.yu;
			slae_in.A.LYF(tmp); r = slae_in.A.yl;
			norm_r = r.norm();
			my_logger.send_current_information(norm_r / norm_f, k_iter);
			my_logger.send_current_information_to_screen(norm_r, k_iter);
		}
		slae_in.A.UXY(x); x = slae_in.A.yu;
		return x;
	}

	void GMRES::solve_min_sqr_problem(MyVector d, DenseMatrix H, MyVector &result)
	{
		int m2 = H.n_columns;
		DenseMatrix H_previous(m2 + 1, m2), H2(s_parameters.gmres_m, s_parameters.gmres_m);
		MyVector d1(m2 + 1), d2(s_parameters.gmres_m);
		double ci, si, tmp;

		double *tmp1, *tmp2;
		tmp1 = new double[s_parameters.gmres_m];
		tmp2 = new double[s_parameters.gmres_m];
		double tmp11, tmp22;

		H_previous = H;
		d1 = d;
		for(int i = 0; i < s_parameters.gmres_m; i++)
		{
			tmp = sqrt(H_previous[i][i] * H_previous[i][i] +
										H[i][i + 1] * H[i][i + 1]);
			ci = H_previous[i][i] / tmp;
			si = H[i][i + 1] / tmp;

			#pragma region H_prev=R*H_prev
			//расчитываем заранее элементы строк, где блок синусов-косинусов
			for(int l = 0; l < s_parameters.gmres_m; l++)
			{
				tmp1[l] = H_previous[l][i] * ci + H_previous[l][i + 1] * si;
				tmp2[l] = -H_previous[l][i] * si + H_previous[l][i + 1] * ci;
			}

			//заполняем строки,где блок синусов-косинусов
			for(int l = 0; l < s_parameters.gmres_m; l++)
			{

				H_previous[l][i] = tmp1[l];
				H_previous[l][i + 1] = tmp2[l];
			}

			#pragma endregion

			#pragma region d1=R*d1
			tmp11 = 0;
			//рассчитываем элементы вектора, где блок синусов-косинусов
			tmp11 = d1[i] * ci + d1[i + 1] * si;
			tmp22 = -d1[i] * si + d1[i + 1] * ci;

			//заполняем элементы вектора, где блок синусов-косинусов
			//остальные не изменяются
			d1[i] = tmp11;
			d1[i + 1] = tmp22;
			#pragma endregion
		}

		//исключаем m+1-ю строку и m+1-й элемент
		for(int j = 0; j < s_parameters.gmres_m; j++)
		{
			for(int i = 0; i < s_parameters.gmres_m; i++)
				H2[j][i] = H_previous[j][i];
			d2[j] = d1[j];
		}

		//находим неизвестный вектор из СЛАУ H2*result=d2
		for(int i = s_parameters.gmres_m - 1; i >= 0; i--)
		{
			result[i] = d2[i];
			for(int j = i + 1; j < s_parameters.gmres_m; j++)
			{		
				result[i] -= result[j] * H2[j][i];
			}
			result[i] = result[i] / H2[i][i];
		}
	}
}
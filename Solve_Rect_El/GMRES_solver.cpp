#include "solver.h"

using namespace myvector;
using namespace densematrix;

namespace solver
{
	MyVector GMRES::Solve(MyVector U_begin, double &normL2u, double &normL2p, MyVector& solution)
	{
		MyVector r(n), x(n), d(m + 1), z(m), w(n), tmp(n), f(n), rr2(n);
		DenseMatrix V(n, m), H(m + 1, m);

		double norm_r, norm_f;
		bool continue_;

		x = U_begin;
		f = A.b;
		A.LU();
		r = f - A * x;
		A.LYF(r); r = A.yl;
		norm_r = r.norm();
		A.LYF(f); norm_f = A.yl.norm();
		x = A.Uv(U_begin);

		for(int k_iter = 1; k_iter <= max_iter && norm_r / norm_f > eps; k_iter++)
		{
			d.make_zero();
			V[0] = r / norm_r;

			continue_ = true;
			for(int j = 1; j <= m && continue_; j++)
			{
				A.UXY(V[j - 1]);
				tmp = A * A.yu;
				A.LYF(tmp); w = A.yl;

				for(int l = 1; l <= j; l++)
				{
					H[j - 1][l - 1] = scal(V[l - 1], w);
					w = w - V[l - 1] * H[j - 1][l - 1];
				}

				H[j - 1][j] = w.norm();
				if(abs(H[j - 1][j]) < 1E-14)
				{
					m = j;
					continue_ = false;
				}
				else
				{
					if(m != j)
						V[j] = w / H[j - 1][j];
				}
			}

			d[0] = norm_r;
			solve_min_sqr_problem(d, H, z);
			x = x + V * z;
			A.UXY(x);
			tmp = f - A * A.yu;
			A.LYF(tmp); r = A.yl;
			norm_r = r.norm();
			logger.send_current_information(norm_r / norm_f, k_iter);
			printf("%d\tr=%.10e\n", k_iter, norm_r / norm_f);
		}
		A.UXY(x); x = A.yu;
		solution = x;
	}

	void GMRES::solve_min_sqr_problem(MyVector d, DenseMatrix H, MyVector &result)
	{
		int m2 = H.n_columns;
		DenseMatrix H_previous(m2 + 1, m2), H2(m, m);
		MyVector d1(m2 + 1), d2(m);
		double ci, si, tmp;

		double *tmp1, *tmp2;
		tmp1 = new double[m];
		tmp2 = new double[m];
		double tmp11, tmp22;

		H_previous = H;
		d1 = d;
		for(int i = 0; i < m; i++)
		{
			tmp = sqrt(H_previous[i][i] * H_previous[i][i] +
										H[i][i + 1] * H[i][i + 1]);
			ci = H_previous[i][i] / tmp;
			si = H[i][i + 1] / tmp;

			#pragma region H_prev=R*H_prev
			//расчитываем заранее элементы строк, где блок синусов-косинусов
			for(int l = 0; l < m; l++)
			{
				tmp1[l] = H_previous[l][i] * ci + H_previous[l][i + 1] * si;
				tmp2[l] = -H_previous[l][i] * si + H_previous[l][i + 1] * ci;
			}

			//заполняем строки,где блок синусов-косинусов
			for(int l = 0; l < m; l++)
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
		for(int j = 0; j < m; j++)
		{
			for(int i = 0; i < m; i++)
				H2[j][i] = H_previous[j][i];
			d2[j] = d1[j];
		}

		//находим неизвестный вектор из СЛАУ H2*result=d2
		for(int i = m-1; i >= 0; i--)
		{
			result[i] = d2[i];
			for(int j = i + 1; j < m; j++)
			{		
				result[i] -= result[j] * H2[j][i];
			}
			result[i] = result[i] / H2[i][i];
		}
	}
}
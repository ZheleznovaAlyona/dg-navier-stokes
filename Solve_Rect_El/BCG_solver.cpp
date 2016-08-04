#include "solver.h"

using namespace myvector;

namespace solver
{
	MyVector BCG::Solve(MyVector U_begin, double &normL2u, double &normL2p, MyVector& solution)
	{
		MyVector r(n), r_(n), p(n), p_(n), f(n);
		MyVector v1(n), v2(n), v3(n), x(n);
		double alpha, betta, old_r_norm = 1.e+20, sc1, sc2;
		double r_norm, r_norm_;
		int k_it;

		A.LU();

		f = A.b;

		k_it = 0;
		r_norm = old_r_norm / 10;
		r_norm_ = 0;

		x = U_begin;
		v1 = A * x;
		v1 = f - v1;
    
		A.LYF(v1); r = A.yl;
		r_ = r;
		p = r;
		p_ = r_;
		int flag = 0;

		while(flag == 0 && k_it < max_iter)
		{
			sc1 = scal(r, r_);
			A.UXY(p); v1 = A.yu;
			v2 = A * v1;
			A.LYF(v2); v3 = A.yl;

			sc2 = scal(p_, v3);

			alpha = sc1 / sc2;
			x = x +  v1 * alpha;
			r = r - v3 * alpha;

			A.LYFt(p_); v1 = A.yl;

			v2 = A  /  v1;

			A.UXYt(v2); v3 = A.yu;

			r_ = r_ - v3 * alpha;
			sc2 = scal(r, r_);

			betta = sc2 / sc1;

			old_r_norm = r_norm;
			r_norm = r.norm();
			printf("%d\tr=%.10e\n", k_it, r_norm);

			logger.send_current_information(r_norm, k_it);
			if( r_norm < eps) flag = 1;

			if(flag == 0)
			{
				p = r + p * betta;
				p_ = r_ + p_ * betta;

				k_it++;
			}
		}

		printf ("\nk_iterations: %ld	r:%.6e", k_it, r_norm);
		if(flag == 0) printf ("\nexit r\n");
		if(flag == 1) printf ("\nsolution end\n");

		if(flag == 2) printf ("\nsolution not end!\n- change vector R_\n");
		solution = x;
	}
}
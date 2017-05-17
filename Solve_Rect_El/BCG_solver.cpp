#include "solver.h"
#include "myfunctions.h"
#include <iostream>
#include <iomanip>

using namespace myvector;
using namespace logger;
using namespace std;

namespace solver
{
	MyVector BCG::solve(MyVector& U_begin, double &normL2u, double &normL2p, slae::SLAE& slae_in, Logger& my_logger)
	{
		MyVector r(slae_in.n), r_(slae_in.n), p(slae_in.n), p_(slae_in.n), f(slae_in.n);
		MyVector v1(slae_in.n), v2(slae_in.n), v3(slae_in.n), x(slae_in.n);
		double alpha, betta, old_r_norm = 1.e+20, sc1, sc2;
		double r_norm, r_norm_;
		int k_it;

		f = slae_in.b;

		k_it = 0;
		r_norm = old_r_norm / 10;
		r_norm_ = 0;

		x = U_begin;
		v1 = slae_in.A * x;
		v1 = f - v1;
    
		slae_in.A.LYF(v1); r = slae_in.A.yl;
		r_ = r;
		p = r;
		p_ = r_;
		int flag = 0;

		while(flag == 0 && k_it < s_parameters.max_number_of_iterations)
		{
			sc1 = scal(r, r_);
			slae_in.A.UXY(p); v1 = slae_in.A.yu;
			v2 = slae_in.A * v1;
			slae_in.A.LYF(v2); v3 = slae_in.A.yl;

			sc2 = scal(p_, v3);

			alpha = sc1 / sc2;
			x = x +  v1 * alpha;
			r = r - v3 * alpha;

			slae_in.A.LYFt(p_); v1 = slae_in.A.yl;

			v2 = slae_in.A  /  v1;

			slae_in.A.UXYt(v2); v3 = slae_in.A.yu;

			r_ = r_ - v3 * alpha;
			sc2 = scal(r, r_);

			betta = sc2 / sc1;

			old_r_norm = r_norm;
			r_norm = r.norm();
			
			my_logger.send_current_information_to_screen(r_norm, k_it);
			my_logger.send_current_information(r_norm, k_it);

			if( r_norm < s_parameters.epsilon) flag = 1;

			if(flag == 0)
			{
				p = r + p * betta;
				p_ = r_ + p_ * betta;

				k_it++;
			}
		}

		return x;
	}
}
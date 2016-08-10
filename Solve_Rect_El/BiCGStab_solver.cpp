#include "solver.h"
#include "myfunctions.h"

using namespace myvector;
using namespace logger;

namespace solver
{
	MyVector BiCGStab::solve(MyVector U_begin, double &normL2u, double &normL2p, slae::SLAE& slae_in, Logger& my_logger)
	{
		double  rkr0, ak, gk, bk;
		MyVector r(slae_in.n), f(slae_in.n), x(slae_in.n), r0(slae_in.n), z(slae_in.n), p(slae_in.n), v(slae_in.n), v1(slae_in.n), rr2(slae_in.n);
		double r_norm, f_norm;

		slae_in.A.LU();

		x = U_begin;
		f = slae_in.b;
		f_norm = f.norm();
		r0 = f - slae_in.A * x;
		slae_in.A.LYF(r0); r0 = slae_in.A.yl;
		r_norm = r0.norm() / f_norm;

		slae_in.A.UXY(r0); z = slae_in.A.yu;
		r = r0;

		my_logger.send_current_information(r_norm, 0);

		for(int k_it = 1; k_it <= s_parameters.max_number_of_iterations && r_norm > s_parameters.epsilon; k_it++)
		{
			//найдем L^(-1)AU^(-1)zk
			slae_in.A.UXY(z); v = slae_in.A.yu;// v = U(-1)zk
			v1 = slae_in.A * v; // v1 = AU^(-1)zk

			slae_in.A.LYF(v1); v = slae_in.A.yl;// v = L^(-1)AU^(-1)zk

			rkr0 = scal(r, r0);
			ak = rkr0 / scal(v, r0); // ak = (r,r0)/ ( L^(-1)AU^(-1)zk,r0)

			p = r - v * ak; // pk = r - ak*L^(-1)AU^(-1)zk

			//найдем L^(-1)AU^(-1)pk
			slae_in.A.UXY(p); v1 = slae_in.A.yu; // v1 = U^(-1)pk
			slae_in.A.LYF(slae_in.A * v1); v1 = slae_in.A.yl; // v1 = L^(-1)AU^(-1)pk

			gk = scal(v1, p) / scal(v1, v1); // gk = (L^(-1)AU^(-1)pk,pk)/ ( L^(-1)AU^(-1)pk,L^(-1)AU^(-1)pk)

			//x = x + ak * zk + gk * pk;
			x = x + z * ak + p * gk;

			r = p - v1 * gk; //r = pk - gk * v1;

			bk = scal(r, r0) / rkr0  * (ak / gk);

			//zk = rk + bk * zk - gk * bk * v;
			z = r + z * bk - v * (gk * bk);

			r_norm = r.norm() / f_norm;
			printf("%d\tr=%.10e\n", k_it, r_norm);
			my_logger.send_current_information(r_norm, k_it);
		}
		slae_in.A.UXY(x); x = slae_in.A.yu;
		return x;
	}
}
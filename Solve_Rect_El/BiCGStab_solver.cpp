#include "solver.h"

using namespace myvector;

namespace solver
{
	MyVector BiCGStab::Solve(MyVector U_begin, double &normL2u, double &normL2p, MyVector& solution)
	{
		double  rkr0, ak, gk, bk;
		MyVector r(n), f(n), x(n), r0(n), z(n), p(n), v(n), v1(n), rr2(n);
		double r_norm, f_norm;

		A.LU();

		x = U_begin;
		f = A.b;
		f_norm = f.norm();
		r0 = f - A * x;
		A.LYF(r0); r0 = A.yl;
		r_norm = r0.norm() / f_norm;

		A.UXY(r0); z = A.yu;
		r = r0;

		logger.send_current_information(r_norm, 0);

		for(int k_it = 1; k_it <= max_iter && r_norm > eps; k_it++)
		{
			//найдем L^(-1)AU^(-1)zk
			A.UXY(z); v = A.yu;// v = U(-1)zk
			v1 = A * v; // v1 = AU^(-1)zk

			A.LYF(v1); v = A.yl;// v = L^(-1)AU^(-1)zk

			rkr0 = scal(r, r0);
			ak = rkr0 / scal(v, r0); // ak = (r,r0)/ ( L^(-1)AU^(-1)zk,r0)

			p = r - v * ak; // pk = r - ak*L^(-1)AU^(-1)zk

			//найдем L^(-1)AU^(-1)pk
			A.UXY(p); v1 = A.yu; // v1 = U^(-1)pk
			A.LYF(A * v1); v1 = A.yl; // v1 = L^(-1)AU^(-1)pk

			gk = scal(v1, p) / scal(v1, v1); // gk = (L^(-1)AU^(-1)pk,pk)/ ( L^(-1)AU^(-1)pk,L^(-1)AU^(-1)pk)

			//x = x + ak * zk + gk * pk;
			x = x + z * ak + p * gk;

			r = p - v1 * gk; //r = pk - gk * v1;

			bk = scal(r, r0) / rkr0  * (ak / gk);

			//zk = rk + bk * zk - gk * bk * v;
			z = r + z * bk - v * (gk * bk);

			r_norm = r.norm() / f_norm;
			printf("%d\tr=%.10e\n", k_it, r_norm);
			logger.send_current_information(r_norm, k_it);
		}
		A.UXY(x); x = A.yu;
		solution = x;
	}
}
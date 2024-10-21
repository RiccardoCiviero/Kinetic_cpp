#include "bertoli.h"

using Eigen::ArrayXd;

ArrayXd Bertoli(Profile p, int i)
{
	i++;
	double delta_res = 5e5;
	double delta = 5e8; //arbitrary big
	ArrayXd nu, f, b, h, eps, rho, eps_t, rho_t, rho_tt, Y;

	Y = p.Y(Eigen::seqN(0, i));
	b = p.b(Eigen::seqN(0, i));
	h = p.h(Eigen::seqN(0, i));
	f = p.f(Eigen::seqN(0, i));
	nu = p.nu(Eigen::seqN(0, i));

	ArrayXd dh(Y.size());
	dh.setZero();
	dh(dh.size() - 1) = 1e-9; // Small just not to make the log diverge 

	eps = f; 
	rho.setZero(i);
	auto test = cumsum(h);
	ArrayXd K = Y * b.square() * (1+nu) * (1- nu*1/4)/ (4*EIGEN_PI) * (log((h.sum() - (cumsum(h) - dh) )/b) +1)  ; // TODO capire moduli
	double E[2]; 


	do // Decided to hardcode the loop because of LAM
	{
		// Layers 0 to N-2
		for (auto j = 0; j < i - 1; j++)
		{
			rho_t = rho;
			rho_tt = rho;
			eps_t = eps;
			E[1] = 1e20; // Abitrary big

			// E--
			rho_t(j) = rho_tt(j) - delta;
			rho_t(j + 1) = rho_tt(j) - delta;
			eps_t = f + cumsum(h * b * rho_t);
			E[0] = (eps.square() * Y * h + K * rho * h).sum();
			if (E[0] < E[1]) {
				E[1] = E[0];
				rho = rho_t;
				eps = eps_t;
			}

			// E-=
			rho_t(j) = rho_tt(j) - delta;
			rho_t(j + 1) = rho_tt(j);
			eps_t = f + cumsum(h * b * rho_t);
			E[0] = (eps.square() * Y * h + K * rho * h).sum();
			if (E[0] < E[1]) {
				E[1] = E[0];
				rho = rho_t;
				eps = eps_t;
			}

			// E-+
			rho_t(j) = rho_tt(j) - delta;
			rho_t(j + 1) = rho_tt(j) + delta;
			eps_t = f + cumsum(h * b * rho_t);
			E[0] = (eps.square() * Y * h + K * rho * h).sum();
			if (E[0] < E[1]) {
				E[1] = E[0];
				rho = rho_t;
				eps = eps_t;
			}

			// E=-
			rho_t(j) = rho_tt(j);
			rho_t(j + 1) = rho_tt(j) - delta;
			eps_t = f + cumsum(h * b * rho_t);
			E[0] = (eps.square() * Y * h + K * rho * h).sum();
			if (E[0] < E[1]) {
				E[1] = E[0];
				rho = rho_t;
				eps = eps_t;
			}

			// E==
			rho_t(j) = rho_tt(j);
			rho_t(j + 1) = rho_tt(j);
			eps_t = f + cumsum(h * b * rho_t);
			E[0] = (eps.square() * Y * h + K * rho * h).sum();
			if (E[0] < E[1]) {
				E[1] = E[0];
				rho = rho_t;
				eps = eps_t;
			}

			// E=+
			rho_t(j) = rho_tt(j);
			rho_t(j + 1) = rho_tt(j) + delta;
			eps_t = f + cumsum(h * b * rho_t);
			E[0] = (eps.square() * Y * h + K * rho * h).sum();
			if (E[0] < E[1]) {
				E[1] = E[0];
				rho = rho_t;
				eps = eps_t;
			}

			// E+-
			rho_t(j) = rho_tt(j) + delta;
			rho_t(j + 1) = rho_tt(j) - delta;
			eps_t = f + cumsum(h * b * rho_t);
			E[0] = (eps.square() * Y * h + K * rho * h).sum();
			if (E[0] < E[1]) {
				E[1] = E[0];
				rho = rho_t;
				eps = eps_t;
			}

			// E+=
			rho_t(j) = rho_tt(j) + delta;
			rho_t(j + 1) = rho_tt(j);
			eps_t = f + cumsum(h * b * rho_t);
			E[0] = (eps.square() * Y * h + K * rho * h).sum();
			if (E[0] < E[1]) {
				E[1] = E[0];
				rho = rho_t;
				eps = eps_t;
			}
			// E++
			rho_t(j) = rho_tt(j) + delta;
			rho_t(j + 1) = rho_tt(j) + delta;
			eps_t = f + cumsum(h * b * rho_t);
			E[0] = (eps.square() * Y * h + K * rho * h).sum();
			if (E[0] < E[1]) {
				E[1] = E[0];
				rho = rho_t;
				eps = eps_t;
			}
		}

		//Layer N-1 (there is no j+1 layer to consider for LAM
		rho_t = rho;
		rho_tt = rho;
		eps_t = eps;
		E[1] = 1e20; // Arbitrary big
		auto j = i ? i - 1 : 0;
		
		// E-
		rho_t(j) = rho_tt(j) - delta;
		eps_t = f + cumsum(h * b * rho_t);
		E[0] = (eps.square() * Y * h + K * rho * h).sum();
		if (E[0] < E[1]) {
			E[1] = E[0];
			rho = rho_t;
			eps = eps_t;
		}

		// E=
		rho_t(j) = rho_tt(j);
		eps_t = f + cumsum(h * b * rho_t);
		E[0] = (eps.square() * Y * h + K * rho * h).sum();
		if (E[0] < E[1]) {
			E[1] = E[0];
			rho = rho_t;
			eps = eps_t;
		}

		// E+
		rho_t(j) = rho_tt(j) + delta;
		eps_t = f + cumsum(h * b * rho_t);
		E[0] = (eps.square() * Y * h + K * rho * h).sum();
		if (E[0] < E[1]) {
			E[1] = E[0];
			rho = rho_t;
			eps = eps_t;
		}

		delta *= 0.9;
	} while (delta > delta_res);
	//std::cout << eps(i-1) << std::endl;
    return eps;
}

#pragma once
#include "kinetic_model.h"


double constexpr schmidt = 0.4082482905; // 1/sqrt(6) see Fitzgerald Dislocations... 1991
double constexpr sinsin = 0.5; // sin(alpha)sin(phi) see Ayers Heteroepitaxy... 1991
double constexpr kb = 8.617333262e-5; // eV/K


using Eigen::ArrayXd;

void kinetic_calculation(Profile p, Results r)
{
	// Working variables
	auto size = p.t.rows();
	ArrayXd tau_eff(size);
	ArrayXd eps(size);
	ArrayXd eps_eq(size);
	ArrayXd dgamma(size);
	ArrayXd rho(size);

	// Placeholder for physical parameters I don't know yet
	double U = 0;
	double rho_0 = 0;
	double K = 0;
	double B = 0;

	// Layer 0
	eps_eq = Bertoli(p, 0);
	eps(0) = p.f(0);
	tau_eff(0) = p.Y(0) * schmidt * (eps(0) - eps_eq(0));
	dgamma(0) = sgn(eps_eq(0) - eps(0)) * K * B * p.b(0) * sinsin * tau_eff(0) * tau_eff(0) * exp(-U / (kb * p.T(0))) * rho_0 * p.h(0) * p.t(0);
	eps(0) += dgamma(0) ;
	rho(0) = (eps(0) - p.f(0)) * p.h(0) / (p.b(0) * sinsin);
	 
	for (auto j = 1; j < size - 1; j++) {

		eps(j) = eps(j - 1) + p.f(j) - p.f(j-1);
		
		eps += p.th_exp_coeff_diff * (p.T(j) - p.T(j-1)); // Should be integrated between t(j-1) and t(j) but it would require explicit T dependence f coeff, here considered constant
		eps_eq = Bertoli(p, j);

		for (auto n = 0; n < j; n++) {
			tau_eff(n) = schmidt * p.h.segment(n,j-n).sum() * (p.Y(n) * (eps(n)-eps_eq(n))) / p.h.segment(n,j-n).sum(); //TODO
		}
		dgamma = sgn(eps_eq(0) - eps(0)) * K * B * p.b * sinsin * tau_eff.square() * exp(-U / (kb * p.T))* p.t * cumsum((rho + rho_0)*p.h);
		eps = eps + dgamma;
		
		rho(0) = (eps(0) - p.f(0))/(p.b(0) * p.h(0) * sinsin);
		for (auto i = 1; i < j; i++)
			rho(i) = (eps(0) - eps(i-1)  - p.f(i) - p.f(i-1) )/ (p.b(0) * p.h(0) * sinsin);
	}

	ArrayXd dLmd = 2 * B * tau_eff * exp(-U / (kb * p.T));

	r.Lmd = cumsum(dLmd);
	r.misfit_density = rho;
	r.strain_profile = eps;

}

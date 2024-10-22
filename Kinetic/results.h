#pragma once
#include <Eigen>

class Results
{
public:
	void save();
	Results(int i) : t(i), a(i), x(i), T(i), strain_profile(i), misfit_density(i), Lmd(i), D(i) {
	}
	Eigen::ArrayXd t, a, x, T, strain_profile, misfit_density, Lmd, D, strain_eq;
};


#pragma once
#include <Eigen>

class Results
{
public:
	void save();
	
	Eigen::ArrayXd t, a, x, T, strain_profile, misfit_density, Lmd, D;
};


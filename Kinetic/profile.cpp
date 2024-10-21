#include "profile.h"
#include<string>
#include<fstream>
#include<iostream>
#include<sstream>
#include<vector>

void Profile::loadCSV(const std::string& filePath)
{
	std::vector<double> p_t[10];

	std::ifstream file(filePath);
	std::string line;

	while (file >> line) {
		std::stringstream ss(line);
		std::string cell;
		int i = 0;
		while (std::getline(ss, cell, ',')) {
			p_t[i].push_back(std::stod(cell));
			i++;
		}

		t = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(p_t[0].data(), p_t[0].size());
		a = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(p_t[1].data(), p_t[1].size());
		x = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(p_t[2].data(), p_t[2].size());
		Y = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(p_t[3].data(), p_t[3].size());
		nu = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(p_t[4].data(), p_t[4].size());
		f = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(p_t[5].data(), p_t[5].size());
		b = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(p_t[6].data(), p_t[6].size());
		h = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(p_t[7].data(), p_t[7].size());
		T = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(p_t[8].data(), p_t[8].size());
		th_exp_coeff_diff = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(p_t[9].data(), p_t[9].size());

	}

}
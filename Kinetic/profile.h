#pragma once
#include <string>
#include <Eigen>

// Profile is a class initialized by reading a csv file with the growth profile already divided in steps, which must be containing all the 
// required properties of the material for each step: a, x, elastic moduli.  

class Profile
{
public:
	void loadCSV(const std::string& filePath);
    Profile(int i) : t(i), a(i), x(i), Y(i), nu(i), f(i), b(i), h(i), T(i), th_exp_coeff_diff(i) {
    }
	Eigen::ArrayXd t, a, x, Y, nu, f, b, h, T, th_exp_coeff_diff;
};

#include <fstream>

inline int countlines(std::string path) {
    std::ifstream file(path);
    int count = 0;
    std::string line;
    while (std::getline(file, line)) {
        count++;
    }
    return count;
}
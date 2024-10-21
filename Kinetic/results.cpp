#include "results.h"
#include <fstream>
#include <iostream>

void Results::save()
{
	std::string path = "C:\\dev.ricca\\Kinetic_cpp\\results.csv";
    std::ofstream file(path);
    if (file.is_open())
    {
        file << "t,a,x,T,strain_profile,misfit_density,Lmd,D\n";
        for (int i = 0; i < a.size(); i++)
        {
            file << t(i) << "," << a(i) << "," << x(i) << "," << T(i) << "," << strain_profile(i) << "," << misfit_density(i) << "," << Lmd(i) << "," << D(i) << "\n";
        }
        file.close();
		std::cout << "Results saved to" << path << std::endl;
    }
}

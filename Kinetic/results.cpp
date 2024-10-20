#include "results.h"
#include <fstream>

void Results::save()
{
    std::ofstream file("results.csv");
    if (file.is_open())
    {
        file << "t,a,x,T,strain_profile,misfit_density,Lmd,D\n";
        for (int i = 0; i < a.size(); i++)
        {
            file << t[i] << a[i] << "," << x[i] << "," << T[i] << "," << strain_profile[i] << "," << misfit_density[i] << "," << Lmd[i] << "," << D[i] << "\n";
        }
        file.close();
    }
}

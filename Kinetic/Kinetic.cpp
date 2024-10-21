#include <iostream>
#include<string>
#include<fstream>
#include "profile.h"
#include "results.h"
#include "kinetic_model.h"
#include "tdd_calculation.h"
#include <Eigen>

int main(int argc, char* argv[])
{
	
	//if (argv[1]) {
	//	auto size = countlines(argv[1]);
	//	Profile p(size);
	//	p.loadCSV(argv[1]);
	//}
	//else {
		std::string path = "C:\\dev\\Kinetic_cpp\\profile.csv";
		auto size = countlines(path);
		Profile p(size);
		p.loadCSV(path);
	//}


	Results r(size);

	r.a << p.a;
	r.x << p.x;
	r.T << p.T;
	r.t << p.t;

	kinetic_calculation(p, r);

	//TDD(p,r);

	r.save();

}

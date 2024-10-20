#include <iostream>
#include<string>
#include "profile.h"
#include "results.h"
#include "kinetic_model.h"
#include "tdd_calculation.h"

int main(int argc, char* argv[])
{
	
	//if (argv[1]) {
	//	auto size = countlines(argv[1]);
	//	Profile p(size);
	//	p.loadCSV(argv[1]);
	//}
	//else {
		std::string path = "C:\\dev.ricca\\Kinetic_cpp\\profile.csv";
		auto size = countlines(path);
		Profile p(size);
		p.loadCSV(path);
	//}

	std::cout << "Loaded profile: " << p.h.size() << std::endl;

	Results r;

	kinetic_calculation(p, r);

	TDD(p,r);

	r.save();

}
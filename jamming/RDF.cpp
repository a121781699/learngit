#include <iostream>
#include"FIRE.h"
#include"CG.h"
#include<iomanip>
using namespace std;

void boxtest(FIRE& b, read_input& file)
{
	b.i = 1;
	b.FIREinitial(file.initial_configuration, file.ignore_line);
	b.NVEdsys();
}


int main(int argc, char* argv[])
{
	read_input file;
	int error = file.read(argc, argv);
	std::cout << std::setprecision(16) << std::endl;
	std::cout << "=====Input FIRE Parameter=====" << std::endl;
	int error0 = file.read_FIRE(argc, argv);

	FIRE b(file.initial_configuration, file.ignore_line, file.N, file.box_length, file.phi, file.Nmin);
	//std::cout << std::setprecision(16) << std::endl;


	double fi = 0;
	boxtest(b,file);

	//b.FIREinitial(file.initial_configuration, file.ignore_line);

	//b.NVEdsys();
}


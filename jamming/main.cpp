//#include "pch.h"
#include <iostream>
#include<fstream>
#include"box.h"
using namespace std;

void boxtest(box* Box, read_input* file)
{
	//ofstream filename;
	//filename.open("Information.dat");

	Box->read_positions(file->initial_configuration, file->ignore_line);
	Box->initialize();
	Box->inflate_volume(0.58);
	Box->FIRErelax();
	//filename << Box->MOVE.fire.GetCur() << ' ' << Box->MONITOR.Uenergy << Box->SHEAR.GetGamma() << endl;
	Box->InitShear(0.02);
	//for (double i = 0.02; i <= 0.4; i += 0.02)
	{
		Box->XzShearStrain();
		Box->FIRErelax();
		//filename << Box->MOVE.fire.GetCur() << ' ' << Box->MONITOR.Uenergy << Box->SHEAR.GetGamma() << endl;
	}
	//filename.close();
}

int main(int argc,char* argv[])
{
	read_input *file = new read_input;
	int error = file->read(argc, argv);
	if (error == 1)
	{
		std::cerr << "There are some errors" << std::endl;
		exit(0);
	}
	box *Box = new box(file->N, file->monitor, file->MonitorFile);

	boxtest(Box, file);

	delete file;
	delete Box;
	
	
	//CG b(file.initial_configuration, file.ignore_line, file.N, file.box_length, file.phi);
	//double fi = 0;
	//while (b.i <= 100)
	//{
	//	b.CGinitial(file.initial_configuration, file.ignore_line);
	//	b.CGrelax();
	//	b.record_positions(file.final_configuration);
	//	if (b.potential_energy / file.N > 1e-16)
	//		++fi;
	//	++b.i;
	//}
	//fi = fi / 100;
	//cout << fi << endl;

	return 0;
}

//#include "pch.h"
#include "read_input.h"
#include<fstream>
#include<iostream>

int read_input::read(int argc, char* argv[])
{
	int error = 0;
	std::ifstream infile;
	infile.open(argv[1]);
	if (!infile)
	{
		std::cerr << "can't open this file" << std::endl;
		error = 1;
		return error;
	}

	char buf[500], c;
	infile.get(buf, 500, '='); infile.get(c); infile >> N;
	infile.get(buf, 500, '='); infile.get(c); infile >> ignore_line;
	infile.get(buf, 500, '='); infile.get(c); infile.width(NAME_LENGTH - 1); infile >> initial_configuration;
	infile.get(buf, 500, '='); infile.get(c); infile.width(NAME_LENGTH - 1); infile >> final_configuration;
	infile.get(buf, 500, '='); infile.get(c); infile >> monitor;
	infile.get(buf, 500, '='); infile.get(c); infile.width(NAME_LENGTH - 1); infile >> MonitorFile;

	std::cout << "N:  " << N << std::endl;
	std::cout << "ignore_line:  " << ignore_line << std::endl;
	std::cout << "initial_configuration:  " << initial_configuration << std::endl;
	std::cout << "final_configuration:  " << final_configuration << std::endl;
	std::cout << "monitor: " << monitor << std::endl;
	std::cout << "MonitorFile: " << MonitorFile << std::endl;

	infile.close();
	return error;
}

int read_input::read_FIRE(int argc, char* argv[])
{
	int error = 0;
	std::ifstream infile;
	if (argc < 3)
	{
		std::cerr << "not find FIREinput script" << std::endl;
		error = 1;
		return error;
	}
	infile.open(argv[2]);
	if (!infile)
	{
		std::cerr << "can't open this file" << std::endl;
		error = 2;
		return error;
	}

	char buf[500], c;
	infile.get(buf, 500, '='); infile.get(c); infile >> phi;
	infile.get(buf, 500, '='); infile.get(c); infile >> Nmin;
	//infile.get(buf, 500, '='); infile.get(c); infile >> ignore_line;
	//infile.get(buf, 500, '='); infile.get(c); infile.width(NAME_LENGTH - 1); infile >> initial_configuration;
	//infile.get(buf, 500, '='); infile.get(c); infile.width(NAME_LENGTH - 1); infile >> final_configuration;

	std::cout << "phi :  " << phi << std::endl;
	std::cout << "Nmin:  " << Nmin << std::endl;
	//std::cout << "ignore_line:  " << ignore_line << std::endl;
	//std::cout << "initial_configuration:  " << initial_configuration << std::endl;
	//std::cout << "final_configuration:  " << final_configuration << std::endl;

	infile.close();
	return error;
}

read_input::~read_input()
{
}

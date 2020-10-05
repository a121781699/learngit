#pragma once
#include "box.h"
#include<vector>
using std::vector;

class CG :
	public box
{
public:
	CG(const char* file, int ignore_line);
	CG(const char* file, int ignore_line, int N_i, double box_length,double phi_i);

	void CGinitial(char* file, int ignore_line);
	void CGrelax();
	~CG();
private:
	double phi;
	vector<double> F;
};


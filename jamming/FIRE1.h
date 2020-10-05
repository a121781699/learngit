#pragma once
#include<vector>
#include"Tool.h"
using std::vector;
class FIRE1
{
public:
	FIRE1();

	void GetFIREval(const double& dt, const int& N);
	void InitFIREval(double& dt, const int& N);
	void update_para(double& dt, const int& N, particle* s);
	int GetCur() { return CurrentStep; }

	~FIRE1();
private:
	double P_value(const int& N, particle* s);
	int CurrentStep;
	int NegativeStep;

	double P;
	vector<double> F_norm;
	vector<double> V;

	double cdt0;
	double dt0, dtmax0;
	double dtmax; 
	double finc, falpha, fdec;
	double alpha, alpha0;
	int Nmin;

	
};


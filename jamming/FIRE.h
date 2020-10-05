#pragma once
#include "box.h"
#include<vector>
using std::vector;

class FIRE :
	public box
{
public:
	FIRE(const char* file,int ignore_line);
	FIRE(const char* file, int ignore_line,int N_i, double box_length);
	FIRE(const char* file, int ignore_line,int N_i, double box_length, double phi_i, int N_min);

	void FIREinitial(char* file, int ignore_line);
	inline void upCurrstep() { ++CurrentStep; }
	inline void upNegstep() { NegativeStep = CurrentStep; }
	inline int Differstep() { int differ = CurrentStep - NegativeStep; return differ; }
	double P_value();
	void update_para();

	void FIRErelax(double phi_i);
	void compress(double dphi_i, double dphi_f, double phi_factor);
	void FIREcompress(double phi_i);
	void test();
	
	~FIRE();
protected:
	double phi;

private:
	//FIRE parameter
	int CurrentStep;
	int NegativeStep;

	double P;
	vector<double> F_norm;
	vector<double> V;
	
	double dtmax;
	double finc, falpha, fdec;
	double alpha, alpha0;
	int Nmin;
};



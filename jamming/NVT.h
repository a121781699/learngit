#pragma once
#include"particle.h"
#define XINTEGRATE 1
#define VINTEGRATE 2
class NVT
{
public:
	NVT(double Temperature, int N);

	double Random();
	double Gauss(double sigma, double l0 = 0, double l = 0);
	void Integrate(int SWITCH, double& temp, int Npart, particle* s, double dt);

    //Nos¨¦-Hoover chain with M=2
	void Chain(particle* s, double& Kenergy, double dt, int N);    //calculate force after Chain function
	void Pos_Vel(particle* s, double& Kenergy, double dt, int N);

	~NVT() {}
private:
	double nu;

	double T, L;
	double Q1, Q2, G1, G2;
	double vxi1, vxi2;
	double xi1, xi2;
	double l;
};


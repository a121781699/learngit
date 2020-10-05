#pragma once
#include"particle.h"
#define tmax 200          //maximum correlation time steps
#define it0 10            //intervel time steps between two original time steps
#define t0max 10           //max number of original time steps

class Diffusion
{
public:
	Diffusion(particle* p, double DT, int npart);

	void dif(int SWITCH, int nsamp, particle* p);
	//-4.4.1 Diffusion
	//-Algorithm 8
	double dt;
	int N;
	int ntel;
	int nsamp; //the number of time steps in one sample
	double dtime; //time between two samples
	int ntime[tmax];
	int time0[t0max];
	double vacf[tmax];
	double r2t[tmax];
	int t0;

	void show();

	~Diffusion();
private:
	particle* s;
	vec<DIM,double>** x0;
	vec<DIM,double>** vx0;
};


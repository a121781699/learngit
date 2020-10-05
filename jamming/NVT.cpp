#include "NVT.h"
#include<math.h>
#include<random>

NVT::NVT(double Temperature, int N)
{
	T = Temperature;
	L = 3 * double(N);

	nu = 1.0;

	Q1 = 1000.0;
	Q2 = 1000.0;
	vxi1 = vxi2 = 0;
	xi1 = xi2 = 0;
}

double NVT::Random()
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dist(0, 1);

	double t = dist(gen);
	return t;
}
double NVT::Gauss(double sigma, double l0, double l)
{
	double r = 2;
	double v1, v2;
	while (r >= 1) {
		v1 = 2 * Random() - 1;
		v2 = 2 * Random() - 1;
		r = v1 * v1 + v2 * v2;
	}
	l = v1 * sqrt(-2 * log(r) / r);
	l = l0 + sigma * l;
	return l;
}
//-6.1.1 The Andersen Thermotat
//-Algorithm 15
//========================================================================
void NVT::Integrate(int SWITCH, double& temp, int Npart, particle *s, double dt)
{
	vec<DIM, double> vt, ft2, ft;
	double v2;
	switch (SWITCH)
	{
	case XINTEGRATE:
		for (int i = 0; i < Npart; ++i) {
			vecMulti(vt, dt, s[i].v); 
			vecMulti(ft2, 0.5*dt*dt, s[i].f);
			Add(vt, vt, ft2);
			Add(s[i].x, s[i].x, vt);
			vecMulti(ft, 0.5*dt, s[i].f);
			Add(s[i].v, s[i].v, ft);
		}
		break;
	case VINTEGRATE:
		temp = 0;
		for (int i = 0; i < Npart; ++i) {
			vecMulti(ft, 0.5*dt, s[i].f);
			Add(s[i].v, s[i].v, ft);
			vecDot(v2, s[i].v, s[i].v);
			temp = temp + v2;
		}
		temp = temp * (1.0 / (3 * Npart));
		double CONST_T = 2.0;
		double sigma = sqrt(CONST_T);
		for (int i = 0; i < Npart; ++i) {
			if (Random() < nu*dt) {
				s[i].v[0] = Gauss(sigma);
				s[i].v[1] = Gauss(sigma);
				s[i].v[2] = Gauss(sigma);
			}
		}
	}
}

void NVT::Chain(particle* s,double& Kenergy,double dt,int N)
{
	G2 = (Q1*vxi1*vxi1 - T);
	vxi2 = vxi2 + G2 * (dt / 4.0);
	vxi1 = vxi1 * exp(-vxi2 * dt / 8.0);
	G1 = (2 * Kenergy - L * T) / Q1;
	vxi1 = vxi1 + G1 * dt / 4.0;
	vxi1 = vxi1 * exp(-vxi2 * dt / 8.0);
	xi1 = xi1 + vxi1 * dt / 2.0;
	xi2 = xi2 + vxi2 * dt / 2.0;
	l = exp(-vxi1 * dt / 2.0);
	for (int i = 0; i < N; ++i)
		vecMulti(s[i].v, l, s[i].v);
	Kenergy = Kenergy * l*l;
	vxi1 = vxi1 * exp(-vxi2 * dt / 8.0);
	G1 = (2 * Kenergy - L * T) / Q1;
	vxi1 = vxi1 + G1 * dt / 4.0;
	vxi1 = vxi1 * exp(-vxi2 * dt / 8.0);
	G2 = (Q1*vxi1*vxi1 - T) / Q2;
	vxi2 = vxi2 + G2 * dt / 4.0;

	vec<DIM, double> vt;
	for (int i = 0; i < N; ++i) {
		vecMulti(vt, 0.5*dt, s[i].v);
		Add(s[i].x, s[i].x, vt);
	}
}

void NVT::Pos_Vel(particle* s, double& Kenergy, double dt, int N)
{
	double e = 0;
	vec<DIM, double> ft_m, vt;
	double v2;
	for (int i = 0; i < N; ++i)
	{
		vecMulti(ft_m, dt / s[i].m, s[i].f);
		Add(s[i].v, s[i].v, ft_m);
		Add(s[i].x, s[i].x, vt);
		vecDot(v2, s[i].v, s[i].v);
		e = e + v2*s[i].m / 2.0;
	}
	Kenergy = e;
}
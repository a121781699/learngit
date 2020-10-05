#pragma once
#include"particle.h"
#include"Tool.h"
class Monitor
{
	friend class box;
public:
	Monitor() { Uenergy = Energy = Kenergy = Temperature = aveCoordNumWithRattler = aveCoordNumNoRattler = 0; SphInfoInit = false; }
	
	double Uenergy;
	double Energy;
	double Kenergy;
	double Temperature;
	double volFrac;
	double aveCoordNumWithRattler;
	double aveCoordNumNoRattler;

	double Uenergyi(particle* s, int i, const int& N, const double& epsilon, double* Length, double* invLength);
	double UenergySys(particle* s, const int& N, const double& epsilon, double* Length, double* invLength);

	double Kenergyi(particle* s, int i, const int& N, const double& epsilon, double* Length);
	double KenergySys(particle* s, const int& N, const double& epsilon, double* Length);
	double TemperaSys(particle* s, const int& N, const double& epsilon, double* Length);

	double EnergySys(particle* s, const int& N, const double& epsilon, double* Length, double* invLength);

	void updateSclInfo(particle* s, const int& N, const double& boxVolume, Allparticle& Ap);

	void computerCoordination(particle* s, const int& N, const double& epsilon, double* Length, double* invLength);
	~Monitor() {}
private:
	bool SphInfoInit;
};


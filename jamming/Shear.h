#pragma once
#include"particle.h"
class Shear
{
	friend class box;
public:
	Shear() { delgamma = 0; gamma = 0; }
	double shearPhi;
	bool ShearFLAG;
	void InitShear(double delg);
	void XzShearStrain(particle* s, int N,double* Length);

	double GetGamma() { return gamma; }
	~Shear() {}
private:
	double delgamma;
	double gamma;
};


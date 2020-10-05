#include "Shear.h"


void Shear::InitShear(double delg)
{
	delgamma = delg;
	ShearFLAG = false;
}

void Shear::XzShearStrain(particle* s,int N,double* Length)
{
	gamma += delgamma;
	for (int i = 0; i < N; ++i)
		s[i].x[0] = s[i].x[0] + s[i].x[2] * delgamma;
	Length[4] += delgamma * Length[2];
	ShearFLAG = true;
}




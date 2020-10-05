#include "Allparticle.h"
#include<cmath>


Allparticle::Allparticle()
{
	meanDiameter = 0.0;
	maxDiameterScale = 0.0;
	minDiameterScale = 0.0;
}

double Allparticle::mean_Diameter(particle* s, const int& N)
{
	double sum = 0;
	double mean_d;
	for (int i = 0; i < N; ++i)
		sum += s[i].DiameterScale;
	mean_d = sum / N;
	meanDiameter = mean_d;
	return mean_d;
}
double Allparticle::max_Diameter(particle* s, const int& N)
{
	double maxScale = 0.0;
	for (int iatom = 0; iatom < N; iatom++) 
		maxScale = (maxScale > s[iatom].DiameterScale ? maxScale : s[iatom].DiameterScale);
	maxDiameterScale = maxScale;
	return maxScale;
}
double Allparticle::min_Diameter(particle* s, const int& N)
{
	double minScale = 1E10;
	for (int iatom = 0; iatom < N; iatom++) {
		double tmp = s[iatom].DiameterScale;
		minScale = (minScale < tmp ? minScale : tmp);
	}
	minDiameterScale = minScale;
	return minScale;
}

void Allparticle::updateSphInfo(particle* s, const int& N)
{
	double factor = 0;
	double minScale = 1E10, maxScale = 0.0;
	for (int iatom = 0; iatom < N; iatom++) {
		factor += pow(s[iatom].DiameterScale, 3);
		double tmp = s[iatom].DiameterScale;
		minScale = (minScale < tmp ? minScale : tmp);
		maxScale = (maxScale > tmp ? maxScale : tmp);
	}
	minDiameterScale = minScale;
	maxDiameterScale = maxScale;

	allSphFrac = VolUnitSphere * factor * pow(meanDiameter * 0.5, 3);
}

#pragma once
#include"particle.h"
#define PI 3.1415926535897932
#ifndef VolUnitSphere
#define VolUnitSphere (pow(PI, ((double)(DIM)) / 2.0) / exp(lgamma(1 + ((double)(DIM)) / 2.0)))  // volume pre-factor for sphere
#endif
class Allparticle
{
	friend class box;
	friend class Monitor;
	friend class NebrList;
public:
	Allparticle();
	double mean_Diameter(particle* s, const int& N);
	double max_Diameter(particle* s, const int& N);
	double min_Diameter(particle* s, const int& N);
	
	void updateSphInfo(particle* s, const int& N);
	~Allparticle() {};
private:
	double meanDiameter;
	double maxDiameterScale;
	double minDiameterScale;
	double allSphFrac;
};


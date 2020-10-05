#pragma once
#include"vec.h"
class particle
{
public:
	particle();
	particle(particle& p);
	particle& operator=(const particle& s);

	//properties of particle
	vec<DIM,double> x;   //position
	vec<DIM,double> v;   //velocity
	vec<DIM,double> f;   //force
	double m;            //mass
	double r;            //radius
	double e;            //energy
	double p;            //pressure
	double z;          //number of contact

	double DiameterScale;    //particles' diameter scale in multibody system.  
	                         //diameter_scale=diameter/mean_diameter
	//========================================================================
	//====   Understanding Molecular Simulation                 ==============
	//====                             --- Daan Frenkel         ==============
	//========================================================================
	//-4.2.3 Integrating the Equations of Motion
	//-Algorithm 6
	vec<DIM,double> xm;          //verlet algorithm      xm is used to store position previous time

	~particle();
};


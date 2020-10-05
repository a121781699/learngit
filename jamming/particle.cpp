//#include"pch.h"
#include "particle.h"

particle::particle()
{
	e = 0.0;
	z = 0.0;
	m = 1.0;
	r = 0.06;
	p = 0.0;
	DiameterScale = 1.0;

	for (int i = 0; i < DIM; ++i)
		v[i] = 0;
	for (int i = 0; i < DIM; ++i)
		f[i] = 0;
}

particle& particle::operator=(const particle& s)
{
	x = s.x;
	v = s.v;
	f = s.f;
	m = s.m;
	r = s.r;
	e = s.e;
	z = s.z;
	p = s.p;

	xm = s.xm;
	return *this;
}

particle::particle(particle& s)
{
	x = s.x;
	v = s.v;
	f = s.f;
	m = s.m;
	r = s.r;
	e = s.e;
	z = s.z;
	p = s.p;

	xm = s.xm;
}


particle::~particle()
{
}

#pragma once
#include"FIRE1.h"
class Move
{
public:
	Move() {}

	FIRE1 fire;
	//MD algorithms
	//1.
	void first_integration(particle* s, const int& N, const double& dt);
	void second_integration(particle* s, const int& N, const double& dt);
	//2.
	void Euler_integration(particle* s, const int& N, const double& dt);

	~Move() {}
};


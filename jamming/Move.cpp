#include "Move.h"


//==================================================================================
//                         Molecular dynamics Algorithms
//==================================================================================

//======================twice integration of velocity algorithm=====================
void Move::first_integration(particle* s, const int& N,const double& dt)
{
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < DIM; ++j)
		{
			s[i].v[j] += 0.5 * dt*s[i].f[j] / s[i].m;           //update volecity
			s[i].x[j] += s[i].v[j] * dt;                        //update position
		}

	}
}
void Move::second_integration(particle* s, const int& N, const double& dt)
{
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < DIM; ++j)
			s[i].v[j] += 0.5 * dt * s[i].f[j] / s[i].m;           //update volecity
	}
}
//==================================================================================

//============================Euler Algorithm=======================================
//update xn using vn
void Move::Euler_integration(particle* s, const int& N, const double& dt)
{
	vec<DIM, double> vt, ft;
	for (int i = 0; i < N; ++i) {
		vecMulti(vt, dt, s[i].v);
		Add(s[i].x, s[i].x, vt);
		vecMulti(ft, dt, s[i].f);
		Add(s[i].v, s[i].v, ft);
	}
}
//===================================================================================

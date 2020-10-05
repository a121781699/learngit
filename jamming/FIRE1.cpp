#include "FIRE1.h"

FIRE1::FIRE1()
{
}

void FIRE1::GetFIREval(const double& dt, const int& N)
{
	V.clear();
	F_norm.clear();

	CurrentStep = 0;
	NegativeStep = 0;

	P = 0;

	cdt0 = dt;
	dt0 = dt; dtmax0 = 100 * dt0;
	dtmax = 100 * dt; 
	finc = 1.1;
	fdec = 0.5;
	alpha0 = 0.1;
	falpha = 0.99;
	alpha = alpha0;
	Nmin = 5;

	V.resize(DIM*N, 0);
	F_norm.resize(DIM*N, 0);
}
void FIRE1::InitFIREval(double& dt, const int& N)
{
	CurrentStep = 0;
	NegativeStep = 0;

	P = 0;

	dt = cdt0;
	dt0 = dt; dtmax0 = 100 * dt0;
	dtmax = 100 * dt;
	finc = 1.1;
	fdec = 0.5;
	alpha0 = 0.1;
	falpha = 0.99;
	alpha = alpha0;
}


double FIRE1::P_value(const int& N,particle* s)
{
	//´ý²âÊÔ
	for (int i = 0; i < N; ++i)
		for (int j = 0; j < DIM; ++j)
		{
			V[i*DIM + j] = s[i].v[j];
			F_norm[i*DIM + j] = s[i].f[j];
		}
	double p0 = 0;
	p0 = scalar_product(V, F_norm);
	P = p0;
	return p0;
}

void FIRE1::update_para(double& dt, const int& N,particle* s)
{
	++CurrentStep;
	P = P_value(N, s);
	if (P > 0)
	{
		if (CurrentStep - NegativeStep > Nmin)
		{
			dt *= finc;
			dt = (dt < dtmax ? dt : dtmax);
			alpha *= falpha;
		}

		double VV = 0, FF = 0;
		VV = scalar_product(V, V);
		FF = scalar_product(F_norm, F_norm);

		double VF = sqrt(VV / FF);
		for (int i = 0; i < N; ++i)
			for (int j = 0; j < DIM; ++j)
				s[i].v[j] = (1 - alpha)*V[i*DIM + j] + alpha * F_norm[i*DIM + j] * VF;
	}
	else
	{
		int deltaStep = CurrentStep - NegativeStep;
		if (deltaStep < log(dtmax / dt0) * 10.5 + Nmin) {
			double a = log(dtmax / dt0);
			dtmax = dtmax * 0.9;
			dtmax = (dtmax > (5.0 * cdt0) ? dtmax : (5.0*cdt0));
		}
		else if (deltaStep > 3.0 * log(dtmax / dt0) * 10.5 + Nmin) {
			double a = log(dtmax / dt0);
			dtmax = dtmax * 1.1;
			dtmax = (dtmax < dtmax0 ? dtmax : dtmax0);
		}

		dt *= fdec;
		dt0 = dt;
		NegativeStep = CurrentStep;
		alpha = alpha0;

		for (int i = 0; i < N; ++i)
			for (int j = 0; j < DIM; ++j)
				s[i].v[j] = 0;
	}
}

FIRE1::~FIRE1()
{
}

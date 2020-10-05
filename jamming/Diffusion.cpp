#include "pch.h"
#include<iostream>
#include "Diffusion.h"

Diffusion::Diffusion(particle* p, double DT, int npart)
{
	s = p; dt = DT; N = npart;
	x0 = new vec<DIM,double>*[N];
	for (int i = 0; i < N; ++i)
		x0[i] = new vec<DIM,double>[t0max];
	vx0 = new vec<DIM,double>*[N];
	for (int i = 0; i < N; ++i)
		vx0[i] = new vec<DIM,double>[t0max];
}

//-4.4.1 Diffusion
//-Algorithm 8
//========================================================================
void Diffusion::dif(int SWITCH, int nsamp, particle* p)
{
	if (SWITCH == 0)
	{
		ntel = 0;
		t0 = 0;
		dtime = dt * nsamp;
		for (int i = 0; i < tmax + 1; ++i) {
			ntime[i] = 0;
			vacf[i] = 0;
			r2t[i] = 0;
		}
	}
	else if (SWITCH == 1) {
		s = p;
		++ntel;
		if (ntel%it0 == 0) {
			t0++;
			int tt0 = (t0 - 1) % t0max;
			time0[tt0] = ntel;
			for (int i = 0; i < N; ++i) {
				x0[i][tt0] = s[i].x;
				vx0[i][tt0] = s[i].v;
			}
		}
		for (int i = 0; i < ((t0 < t0max) ? t0 : t0max); ++i) {
			int delt = ntel - time0[i];
			if (delt < tmax) {
				ntime[delt] = ntime[delt] + 1;
				for (int j = 0; j < N; ++j) {
					double f;
					vecDot(f, s[j].v, vx0[j][i]);
					vacf[delt] += f;

					vec<DIM,double> temp;
					Minus(temp, s[j].x, x0[j][i]);
					double x2 = 0;
					for (int k = 0; k < DIM; ++k)
						x2 += temp[k] * temp[k];

					r2t[delt] += x2;
				}
			}
		}
	}
	else if (SWITCH == 2) {
		for (int i = 0; i < tmax; ++i) {
			double time = dtime * (i + 0.5);
			vacf[i] = vacf[i] / (N*ntime[i]);
			r2t[i] = r2t[i] / (N*ntime[i]);
		}
	}
}

void Diffusion::show() {
	for (int i = 0; i < tmax; ++i)
		std::cout << "time: " << (dtime*(i + 0.5)) << ' ' << "velocity autocorr: " << vacf[i] << ' '
		<< "mean-squared dislacement: " << r2t[i] << std::endl;
}
Diffusion::~Diffusion()
{
	for (int i = 0; i < N; ++i)
		delete[]x0[i];
	delete[] x0;
	for (int i = 0; i < N; ++i)
		delete[]vx0[i];
	delete[] vx0;
}

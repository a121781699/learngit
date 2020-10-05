#include "CG.h"




CG::CG(const char* file, int ignore_line):
	box()
{
	initialize();


	for (int i = 0; i < N*DIM; ++i)
		F.push_back(0);
}

CG::CG(const char* file, int ignore_line, int N_i, double box_length,double phi_i):
	box(N_i, box_length),
	phi(phi_i)
{
	initialize();

	for (int i = 0; i < N*DIM; ++i)
		F.push_back(0);
}

void CG::CGinitial(char* file, int ignore_line)
{
	read_positions(file, ignore_line);
}

void CG::CGrelax()
{
	bool stop = false;
	double old_e = -INFINITY;
	inflate_volume(phi);

	for (int i = 0; i < N; ++i)
		s[i].f = force_i(i);
	for (int i = 0; i < N; ++i)
		for (int j = 0; j < DIM; ++j)
			F[i*DIM + j] = s[i].f[j];
	double F0 = scalar_product(F, F);
	for (int i = 0; i < N; ++i)
		s[i].v = force_i(i)*(1.0/sqrt(F0));

	while (!stop) {
		system_property();
		if (potential_energy / N < 1e-16 || fabs(potential_energy - old_e) / N < 1e-15)
			stop = true;
		old_e = potential_energy;
		for (int i = 0; i < N; ++i)
			for (int j = 0; j < DIM; ++j)
				F[i*DIM + j] = s[i].f[j];
		double Fn_old = scalar_product(F, F);
		for (int i = 0; i < N; ++i)
			s[i].x = s[i].x + s[i].v;

		for (int i = 0; i < N; ++i)
			s[i].f = force_i(i);
		for (int i = 0; i < N; ++i)
			for (int j = 0; j < DIM; ++j)
				F[i*DIM + j] = s[i].f[j];
		double Fn = scalar_product(F, F);
		for (int i = 0; i < N; ++i)
			s[i].v = s[i].f + s[i].v*(Fn*(1.0 / Fn_old));
	}
	system_property();
}

CG::~CG()
{
}

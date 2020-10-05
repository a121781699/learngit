#include "Monitor.h"

double Monitor::Uenergyi(particle* s,int i, const int& N,const double& epsilon,double* Length, double* invLength)
{
	double factor = 0.5*epsilon / 2.0;
	double e = 0;
	vec<DIM,double> vec_rij;                      //vec_rij = xi - xj
	double rij, sigmaij;                             //rij = |xi - xj|; sigma_ij = r_i + r_j;
	s[i].z = 0;                                         //coordination number
	for (int j = 0; j < N; j++) {
		if (i == j) continue;
		vec_rij = distance_vector(s[i].x, s[j].x, Length, invLength);
		rij = scalarization(vec_rij);
		sigmaij = s[i].r + s[j].r;
		if (rij < sigmaij) {
			e = e + (1.0 - rij / sigmaij)*(1.0 - rij / sigmaij);
			++s[i].z;              //increase contact number
		}
	}
	e = factor * e;
	s[i].e = e;
	return e;
}
double Monitor::UenergySys(particle* s, const int& N, const double& epsilon, double* Length, double* invLength)
{
	double e_i, e0 = 0;
	for (int i = 0; i < N; i++)
	{
		e_i = Uenergyi(s, i, N, epsilon, Length, invLength);
		e0 += e_i;
	}
	Uenergy = e0;
	return e0;
}


double Monitor::Kenergyi(particle* s, int i, const int& N, const double& epsilon, double* Length)
{
	double k, fac;
	vecDot(fac, s[i].v, s[i].v);
	k = 0.5*s[i].m*fac;
	return k;
}
double Monitor::KenergySys(particle* s, const int& N, const double& epsilon, double* Length)
{
	double e = 0;
	for (int i = 0; i < N; ++i)
		e += Kenergyi(s, i, N, epsilon, Length);
	Kenergy = e;
	return e;
}
double Monitor::TemperaSys(particle* s, const int& N, const double& epsilon, double* Length)
{
	double k, e2 = 0;
	for (int i = 0; i < N; ++i)
	{
		vecDot(k, s[i].v, s[i].v);
		e2 += k;
	}
	e2 = e2 / (3 * N);
	Temperature = e2;
	return e2;
}

double Monitor::EnergySys(particle* s, const int& N, const double& epsilon, double* Length, double* invLength)
{
	double e0 = 0;
	KenergySys(s, N, epsilon, Length);
	UenergySys(s, N, epsilon, Length, invLength);

	e0 = Kenergy + Uenergy;
	Energy = e0;
	return e0;
}

void Monitor::updateSclInfo(particle* s, const int& N, const double& boxVolume, Allparticle& Ap)
{
	if (SphInfoInit == true) return;
	Ap.updateSphInfo(s, N);
	volFrac = Ap.allSphFrac / boxVolume;
	SphInfoInit = true;
}

void Monitor::computerCoordination(particle* s, const int& N, const double& epsilon, double* Length, double* invLength)
{
	int *nCoordNumWithRattler = new int[N];
	int *nCoordNumNoRattler = new int[N];
	
	for (int iatom = 0; iatom < N; ++iatom)
	{
		nCoordNumNoRattler[iatom] = 0;
		nCoordNumWithRattler[iatom] = 0;
	}

	for (int iatom = 0; iatom < N; ++iatom)
	{
		int nwr = nCoordNumWithRattler[iatom];
		for (int jatom = iatom + 1; jatom < N; ++jatom)
		{
			vec<DIM,double> vec = distance_vector(s[iatom].x, s[jatom].x, Length, invLength);
			double rij = scalarization(vec);
			double sRc = s[iatom].r + s[jatom].r;
			if (rij < sRc)
			{
				++nwr;
				++nCoordNumWithRattler[jatom];
			}
		}
		nCoordNumWithRattler[iatom] = nwr;
	}
	for (int iatom = 0; iatom < N; ++iatom)
		nCoordNumNoRattler[iatom] = nCoordNumWithRattler[iatom];

	int cordNumWR = 0, cordNumNR = 0, nRattler = 0;
	for (int iatom = 0; iatom < N; iatom++) {
		if (nCoordNumNoRattler[iatom] < (DIM + 1)) {
			nCoordNumNoRattler[iatom] = 0;
			nRattler++;
		}
		cordNumWR += nCoordNumWithRattler[iatom];
		cordNumNR += nCoordNumNoRattler[iatom];
	}
	aveCoordNumWithRattler = (double)cordNumWR / (double)N;
	if (nRattler != N) {
		aveCoordNumNoRattler =
			(double)cordNumNR / (double)(N - nRattler);
	}
	else {
		aveCoordNumNoRattler = 0;
	}

	delete[] nCoordNumWithRattler;
	delete[] nCoordNumNoRattler;
}

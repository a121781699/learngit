#include "Force.h"
#include<time.h>

vec<DIM,double> Force::HarmonicForce_ij(particle* s, int i, int j, bool& zeroforce, double* Length, double* invLength,const double& epsilon)
{
	vec<DIM,double> fij;
	vec<DIM,double> vec_ij, norm_ij;
	double rij, rij2, sigmaij;
	vec_ij = distance_vector(s[i].x, s[j].x, Length, invLength);
	rij2 = square(vec_ij);
	sigmaij = s[i].r + s[j].r;
	if (rij2 < sigmaij*sigmaij)
	{
		rij = sqrt(rij2);
		vecDevi(norm_ij, vec_ij, rij);
		double factor = epsilon / sigmaij * (1.0 - rij / sigmaij);
		vecMulti(fij, factor, norm_ij);
		zeroforce = false;
	}
	else
		zeroforce = true;
	return fij;
}

void Force::HarmonicForce(particle* s, const double& epsilon, const int& N, double* Length,double* invLength, Monitor& MONITOR)
{
	zeroforceratio = 0.0;
	for (int i = 0; i < N; ++i)
		for (int j = 0; j < DIM; ++j)
			s[i].f[j] = 0;
	MONITOR.Uenergy = 0;
	for (int i = 0; i < N; ++i) {
		vec<DIM,double> iforce = s[i].f;
		vec<DIM,double> iPos = s[i].x;
		double iRc = s[i].r;
		for (int j = i + 1; j < N; j++) {
			if (i >= j) continue;

			vec<DIM,double> vec_ij;
			double sigmaij = iRc + s[j].r;
			vec_ij = distance_vector(s[i].x, s[j].x, Length, invLength);
			double rij2 = square(vec_ij);
			if (rij2 > sigmaij*sigmaij) continue;
			double rij = sqrt(rij2);
			double fpair = epsilon / sigmaij * (1.0 - rij / sigmaij) / rij;
			vec<DIM,double> H;
			vecMulti(H, fpair, vec_ij);
			Add(iforce, iforce, H);
			vecMulti(H, -fpair, vec_ij);
			Add(s[j].f, s[j].f, H);

			MONITOR.Uenergy += 0.5*epsilon * (1.0 - rij / sigmaij)*(1.0 - rij / sigmaij);
		}
		s[i].f = iforce;
		if (s[i].f[0] < 1e-14&&s[i].f[1] < 1e-14&&s[i].f[2] < 1e-14)
			++zeroforceratio;
	}
	zeroforceratio /= N;
	double AVEFORCE = 0;
	for (int i = 0; i < N; ++i)
		AVEFORCE += scalarization(s[i].f);
	aveforce = AVEFORCE / N;
}

//update--cell list
vec<DIM,double> Force::HarmonicForce_iwithCelllist(particle* s, int i, bool& zeroforce, double* Length, double* invLength, const double& epsilon, NebrList* nebrlist)
{
	int jcell = 0;
	vec<DIM,double> f, fij;
	int celli = (int)((s[i].x[0] + 0.5*Length[0]) / nebrlist->cLen[0]);
	int ci = celli % nebrlist->Lcnum[0];
	if (s[i].x[0] + 0.5*Length[0] < 0)
		celli = ci + nebrlist->Lcnum[0] - 1;
	else celli = ci;
	int cellj = (int)((s[i].x[1] + 0.5*Length[0]) / nebrlist->cLen[1]);
	int cj = cellj % nebrlist->Lcnum[1];
	if (s[i].x[1] + 0.5*Length[0] < 0)
		cellj = cj + nebrlist->Lcnum[1] - 1;
	else cellj = cj;
	int cellk = (int)((s[i].x[2] + 0.5*Length[0]) / nebrlist->cLen[2]);
	int ck = cellk % nebrlist->Lcnum[2];
	if (s[i].x[2] + 0.5 *Length[0] < 0)
		cellk = ck + nebrlist->Lcnum[2] - 1;
	else cellk = ck;
	int icell = celli + cellj * nebrlist->Lcnum[0] + cellk * nebrlist->Lcnum[0] * nebrlist->Lcnum[1];

	int* ptemp = nebrlist->Celllist->First[icell + 1];
	while (ptemp != nebrlist->Celllist->Last[icell + 1])
	{
		jcell = *ptemp;
		int j = 0;
		for (j = nebrlist->Head[jcell]; j != -1; j = nebrlist->List[j])
		{
			if (i == j) continue;
			fij = HarmonicForce_ij(s, i, j, zeroforce, Length, invLength, epsilon);
			if (!zeroforce)
				Add(f, f, fij);
		}
		ptemp++;
	}
	if (f[0] < 1e-14&&f[1] < 1e-14&&f[2] < 1e-14)
		++zeroforceratio;
	return f;
}
void Force::HarmonicForcewithCelllist(particle* s, bool& zeroforce, const double& epsilon, const Allparticle& Ap, const int& N, double* Length, double* invLength, NebrList* nebrlist)
{
	nebrlist->ParticleCell(s, Length, invLength, Ap, N);
	zeroforceratio = 0;
	for (int i = 0; i < N; ++i)
		s[i].f = HarmonicForce_iwithCelllist(s, i, zeroforce, Length, invLength, epsilon, nebrlist);
	zeroforceratio /= N;
}

void Force::HarmonicForcewithNebrlist(particle* s, const double& epsilon, const Allparticle& Ap, const int& N, double* Length, double* invLength, NebrList* nebrlist,Monitor& MONITOR)
{
	//if (nebrlist->IsNebrlistInvalid(s, Length, invLength, sigma_mean, N))
		//nebrlist->BuildNebr(s,Length, invLength, sigma_mean,N);
	bool IsInvalid = nebrlist->IsNebrlistInvalid(s, Length, invLength, Ap, N);
	//std::cout << "IsInvalid:" << IsInvalid << std::endl;
	if (IsInvalid) {
		if (cntForce < 10) {
			nebrlist->skinSet *= 1.1;
			nebrlist->skinSet = ((nebrlist->skinSet < 0.5) ? nebrlist->skinSet : 0.5);
		}
		else if (cntForce > 20) {
			nebrlist->skinSet *= 0.9;
			nebrlist->skinSet = ((nebrlist->skinSet > 5E-3) ? nebrlist->skinSet : 5E-3);
		}
	}
	else {
		if (nebrlist->skinSet > 5E-3 && cntForce > 20) {
			nebrlist->skinSet *= 0.9;
			nebrlist->skinSet = ((nebrlist->skinSet > 5E-3) ? nebrlist->skinSet : 5E-3);
			IsInvalid = true;
		}
	}
	if (IsInvalid) {
		nebrlist->BuildNebr(s, Length, invLength, Ap, N, cntForce);
	}
	
	zeroforceratio = 0.0;
	for (int i = 0; i < N; ++i)
		for (int j = 0; j < DIM; ++j)
			s[i].f[j] = 0;
	MONITOR.Uenergy = 0;

	clock_t start = clock();
	for (int i = 0; i < N; ++i)
	{
		vec<DIM,double> iforce = s[i].f;
		vec<DIM,double> iPos = s[i].x;
		double iRc = s[i].r;
		for (int k=1; k <= nebrlist->nlist[i]; ++k)
		{
			int j = nebrlist->Nnebrlist[INCREMENT*i + k];
			if (j <= i) continue;
			
			double sigmaij = iRc + s[j].r;
			vec<DIM, double> vec_ij;
			Minus(vec_ij, s[i].x, s[j].x);
			periodicVec(vec_ij, Length);
			//vec_ij = distance_vector(s[i].x, s[j].x, Length, invLength);
			double rij2 = square(vec_ij);
			if (rij2 > sigmaij*sigmaij) continue;
			double rij = sqrt(rij2);
			double fpair = epsilon / sigmaij * (1.0 - rij / sigmaij) / rij;
			vec<DIM,double> H;
			vecMulti(H, fpair, vec_ij);
			Add(iforce, iforce, H);
			vecMulti(H, -fpair, vec_ij);
			Add(s[j].f, s[j].f, H);
			
			MONITOR.Uenergy += 0.5*epsilon * (1.0 - rij / sigmaij)*(1.0 - rij / sigmaij);
		}
		s[i].f = iforce;
	}
	clock_t end = clock();
	double time = (double)(end - start) / CLK_TCK;
	std::cout << "time: "<< time << std::endl;
	/*zeroforceratio = 0.0;
	for (int i = 0; i < N; ++i)
		if (s[i].f[0] < 1e-14&&s[i].f[1] < 1e-14&&s[i].f[2] < 1e-14)
			++zeroforceratio;
	zeroforceratio /= N;*/
	++cntForce;
	double AVEFORCE = 0;
	for (int i = 0; i < N; ++i)
		AVEFORCE += scalarization(s[i].f);
	aveforce = AVEFORCE / N;
}


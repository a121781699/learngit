#pragma once
#include"NebrList.h"
#include"Monitor.h"
class Force
{
	friend class box;
public:
	Force() { zeroforce = false; zeroforceratio = 0.0; aveforce = 0.0; cntForce = 0; }

	vec<DIM,double> HarmonicForce_ij(particle* s, int i, int j, bool& zeroforce, double* Length, double* invLength,const double& epsilon);
	vec<DIM,double> HarmonicForce_iwithCelllist(particle* s, int i, bool& zeroforce, double* Length, double* invLength, const double& epsilon, NebrList* nebrlist);
	void HarmonicForcewithCelllist(particle* s, bool& zeroforce, const double& epsilon, const Allparticle& Ap, const int& N, double* Length, double* invLength, NebrList* nebrlist);
    void HarmonicForce(particle* s, const double& epsilon, const int& N, double* Length,double* invLength, Monitor& MONITOR);
	void HarmonicForcewithNebrlist(particle* s, const double& epsilon, const Allparticle& Ap, const int& N, double* Length, double* invLength, NebrList* nebrlist,Monitor& MONITOR);

	~Force() {}
private:
	bool zeroforce;
	double zeroforceratio;
	double aveforce;
	int cntForce;
};


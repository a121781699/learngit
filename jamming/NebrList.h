#pragma once
#include"Tool.h"
#include"particle.h"

#include"ListArray.h"
#define INCREMENT 50


class NebrList {
	friend class box;
	friend class Force;
public:
	NebrList() { List = 0; Head = 0; sHold = 0; nlist = 0; IsList = false; };
	~NebrList();
	void mallocCell(particle* s, double* Length, const Allparticle& Ap, const int& N);
	void InitCell(particle* s, double* Length, double* invLength, const Allparticle& Ap, const int& N);
	void ParticleCell(particle* s, double* Length, double* invLength, const Allparticle& Ap, const int& N);

	void mallocNebr(particle* s, double* Length, const Allparticle& Ap, const int& N);
	void InitNebr(particle* s, double* Length, double* invLength, const Allparticle& Ap, const int& N);
	void BuildNebr(particle* s, double* Length, double* invLength, const Allparticle& Ap, const int& N, int& cntForce);
	bool IsNebrlistInvalid(particle* s, double* Length, double* invLength, const Allparticle& Ap, const int& N);
private:
	int *List;
	int *Head;
	int Lcnum[DIM];        //the number of cells toward any direction in box
	double cLen[DIM];      //the length of one cell in three directions
	int cnum;
	ListArray<int> *Celllist;

	bool IsList;
	bool IsValid;
	double dskin;
	double skinSet;
	particle* sHold;              //hold positions of particles after updating Nebrlist;
	int* nlist;                   //nlist[i] hold the number of particles in the Verlet list of particle i
	int* Nnebrlist;        //Nnebrlist[i] is Verlet list of particle i
};
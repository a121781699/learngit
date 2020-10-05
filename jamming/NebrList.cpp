#include "pch.h"
#include "NebrList.h"

//allocate space for serial number of particles and serial number of cells
void NebrList::mallocCell(particle* s, double* Length, const Allparticle& Ap, const int& N)
{
	//cell list
	double DistPlane[DIM];
	// L(a,b): distance between plane (a,b); {vol / |a cross b|} = Length[2];
	DistPlane[2] = Length[2];
	// L(c,a): distance between plane (c,a); {vol / |c cross a|}
	double sca = sqrt(pow(Length[0] * Length[2], 2) +
				      pow(Length[0] * Length[3], 2));
	DistPlane[1] = Length[0] * Length[1] * Length[2] / sca;
	// L(b,c): distance between plane (b,c); {vol / |b cross c|}
	double sbc = sqrt(pow(Length[1] * Length[2], 2) +
					  pow(Length[5] * Length[2], 2) +
					  pow(Length[5] * Length[3] - 
						  Length[1] * Length[4], 2));
	DistPlane[0] = Length[0] * Length[1] * Length[2] / sbc;

	dskin = Ap.minDiameterScale*Ap.meanDiameter * skinSet;
	double sysRcs = Ap.maxDiameterScale*Ap.meanDiameter + dskin;
	
	Lcnum[0] = (int)floor(DistPlane[0] / sysRcs);
	Lcnum[1] = (int)floor(DistPlane[1] / sysRcs);
	Lcnum[2] = (int)floor(DistPlane[2] / sysRcs);
	Lcnum[0] = (Lcnum[0] <= 0 ? 1 : Lcnum[0]);
	Lcnum[1] = (Lcnum[1] <= 0 ? 1 : Lcnum[1]);
	Lcnum[2] = (Lcnum[2] <= 0 ? 1 : Lcnum[2]);
	cnum = Lcnum[0] * Lcnum[1] * Lcnum[2];
	if (cnum > N) {
		int num0 = (int)floor(pow(
			N * DistPlane[1] * DistPlane[2] / DistPlane[0] / DistPlane[0],
			1.0 / 3.0));
		int num1 = (int)floor(DistPlane[1] / DistPlane[0] * num0);
		int num2 = (int)floor(DistPlane[2] / DistPlane[0] * num0);

		Lcnum[0] = num0;
		Lcnum[1] = num1;
		Lcnum[2] = num2;
		cnum = Lcnum[0] * Lcnum[1] * Lcnum[2];
	}
	for (int i = 0; i < DIM; ++i)
		cLen[i] = Length[i] / Lcnum[i];
	if(List==NULL)
		List = new int[N];
	if(Head==NULL)
		Head = new int[cnum];
	if(Celllist==NULL)
		Celllist = new ListArray<int>;
}

//set neighboor cells for each cell and reserve them in ListArray
void NebrList::InitCell(particle* s, double* Length, double* invLength, const Allparticle& Ap, const int& N)
{
	mallocCell(s, Length, Ap, N);
	/*vec<DIM,double>* x = new vec<DIM,double>[cnum];
	for (int i = 0; i < Lcnum[2]; ++i)
	{
		for (int j = 0; j < Lcnum[1]; ++j)
		{
			for (int k = 0; k < Lcnum[0]; ++k)
			{
				int inum = i * Lcnum[0] * Lcnum[1] + j * Lcnum[0] + k;
				x[inum][0] = k * cLen[0];
				x[inum][1] = j * cLen[1];
				x[inum][2] = i * cLen[2];
			}
		}
	}

	vec<DIM,double> vec;
	double Rdia2 = cLen[0] * cLen[0] + cLen[1] * cLen[1] + cLen[2] * cLen[2];
	double Rdia = sqrt(Rdia2);
	for (int inum = 0; inum < cnum; ++inum)
	{
		Celllist->Insert(inum + 1, 0, inum);
		for (int jnum = 0; jnum < cnum; ++jnum)
		{
			if (inum == jnum) continue;
			vec = OrthoDistance(x[inum], x[jnum], Length);
			double rij = vec.scalarization();
			rij = (double)(((int)(rij * 10000)) / 10000.0);
			if (rij <= Rdia)
			{
				Celllist->Insert(inum + 1, 0, jnum);
			}
		}
	}
	delete[] x;*/

	/*for (int kbin = 0; kbin < Lcnum[2]; kbin++) {
		for (int jbin = 0; jbin < Lcnum[1]; jbin++) {
			for (int ibin = 0; ibin < Lcnum[0]; ibin++) {
				int binIdx = kbin * Lcnum[1] * Lcnum[0] +
					jbin * Lcnum[0] + ibin;
				Celllist->Insert(binIdx + 1, 0, binIdx);
				{
					double iRcut = s[0].r * MeanSigma;
					double mRcutRskin = iRcut + MeanSigma*MeanSigma + dskin;

					int zNum = (int)ceil(mRcutRskin / cLen[2]);
					int zstart = -zNum, zstop = zNum;
					if (2 * zNum + 1 > Lcnum[2]) {
						zstart = 0;
						zstop = Lcnum[2] - 1;
					}
					for (int bindz = zstart; bindz <= zstop; bindz++) {
						int tbinz = (kbin + bindz + Lcnum[2]) % Lcnum[2];

						int yNum = (int)ceil(mRcutRskin / cLen[1]);
						int ystart = -yNum, ystop = yNum;
						if (2 * yNum + 1 > Lcnum[1]) {
							ystart = 0;
							ystop = Lcnum[1] - 1;
						}
						for (int bindy = ystart; bindy <= ystop; bindy++) {
							int tbiny = (jbin + bindy + Lcnum[1]) % Lcnum[1];

							int xNum = (int)ceil(mRcutRskin / cLen[0]);
							int xstart = -xNum, xstop = xNum;
							if (2 * xNum + 1 > Lcnum[0]) {
								xstart = 0;
								xstop = Lcnum[0] - 1;
							}

							for (int bindx = xstart; bindx <= xstop; bindx++) {
								int tbinx =
									(ibin + bindx + Lcnum[0]) % Lcnum[0];
								int tbinIdx = tbinz * Lcnum[1] * Lcnum[0] +
									tbiny * Lcnum[0] + tbinx;
								Celllist->Insert(tbinIdx + 1, 0, tbinIdx);
							}
						}
					}
				}
			}
		}
	}*/
}

//find out all particles in one cell and link those particles one by one,reserve the biggest 
//serial number of particle in this cell.
void NebrList::ParticleCell(particle *s, double* Length, double* invLength, const Allparticle& Ap, const int& N)
{
	//Initialize Head
	memset(Head, -1, cnum * sizeof(int));
	//set List(put index of particles into List)
	//set Head(reserve the index of last particle of List)
	for (int iatom = 0; iatom < N; iatom++) {
		particle *posXyz = &s[iatom];

		int cIdxx, cIdxy, cIdxz;
		double lamda_x = invLength[0] * posXyz->x[0] +
						 invLength[5] * posXyz->x[1] +
						 invLength[4] * posXyz->x[2] + 0.5;
		double lamda_y = invLength[1] * posXyz->x[1] +
						 invLength[3] * posXyz->x[2] + 0.5;
		double lamda_z = invLength[2] * posXyz->x[2] + 0.5;
		cIdxx = ((int)floor(lamda_x * Lcnum[0])) % Lcnum[0];
		cIdxy = ((int)floor(lamda_y * Lcnum[1])) % Lcnum[1];
		cIdxz = ((int)floor(lamda_z * Lcnum[2])) % Lcnum[2];

		cIdxx = (cIdxx < 0 ? Lcnum[0] - 1 : cIdxx);
		cIdxx = (cIdxx >= Lcnum[0] ? 0 : cIdxx);
		cIdxy = (cIdxy < 0 ? Lcnum[1] - 1 : cIdxy);
		cIdxy = (cIdxy >= Lcnum[1] ? 0 : cIdxy);
		cIdxz = (cIdxz < 0 ? Lcnum[2] - 1 : cIdxz);
		cIdxz = (cIdxz >= Lcnum[2] ? 0 : cIdxz);

		int inumber =
			cIdxz * Lcnum[1] * Lcnum[0] + cIdxy * Lcnum[0] + cIdxx;
		List[iatom] = Head[inumber];
		Head[inumber] = iatom;
	}
}

void NebrList::mallocNebr(particle* s, double* Length, const Allparticle& Ap, const int& N)
{
	if(sHold==NULL)
		sHold = new particle[N];
	if(nlist==NULL)
		nlist = new int[N];
	if(Nnebrlist==NULL)
		Nnebrlist = new int[INCREMENT*N];
}

void NebrList::InitNebr(particle* s, double* Length, double* invLength, const Allparticle& Ap, const int& N)
{
	InitCell(s, Length, invLength, Ap, N);
	mallocNebr(s, Length, Ap, N);
	memset(Nnebrlist, 0, N*INCREMENT * sizeof(int));
	IsList = true;
	IsValid = false;
}

void NebrList::BuildNebr(particle* s, double* Length, double* invLength, const Allparticle& Ap, const int& N, int& cntForce)
{
	InitNebr(s, Length, invLength, Ap, N);
	ParticleCell(s, Length, invLength, Ap, N);
	memset(nlist, 0, N * sizeof(int));
	memcpy(sHold, s, N*sizeof(particle));
	double meanRadius = 0.5 * Ap.meanDiameter;
	double maxRadius = 0.5 * Ap.maxDiameterScale * Ap.meanDiameter;
	for (int kbin = 0; kbin < Lcnum[2]; ++kbin) {        //========
		for (int jbin = 0; jbin < Lcnum[1]; ++jbin) {    //========
			for (int ibin = 0; ibin < Lcnum[0]; ++ibin) {     //=========
				int icell = kbin * Lcnum[0] * Lcnum[1] + jbin * Lcnum[0] + ibin;   //========

				int i = Head[icell];
				while (i != -1) {
					double iRcut = s[i].DiameterScale*meanRadius;
					double mRcutRskin = maxRadius + iRcut + dskin;


					int zNum = (int)ceil(mRcutRskin / cLen[2]);
					int zstart = -zNum, zstop = zNum;
					if (2 * zNum + 1 > Lcnum[2]) {
						zstart = 0;
						zstop = Lcnum[2] - 1;
					}
					for (int bindz = zstart; bindz <= zstop; bindz++) {
						int tbinz = (kbin + bindz + Lcnum[2]) % Lcnum[2];

						int yNum = (int)ceil(mRcutRskin / cLen[1]);
						int ystart = -yNum, ystop = yNum;
						if (2 * yNum + 1 > Lcnum[1]) {
							ystart = 0;
							ystop = Lcnum[1] - 1;
						}
						for (int bindy = ystart; bindy <= ystop; bindy++) {
							int tbiny = (jbin + bindy + Lcnum[1]) % Lcnum[1];

							int xNum = (int)ceil(mRcutRskin / cLen[0]);
							int xstart = -xNum, xstop = xNum;
							if (2 * xNum + 1 > Lcnum[0]) {
								xstart = 0;
								xstop = Lcnum[0] - 1;
							}

							for (int bindx = xstart; bindx <= xstop; bindx++) {
								int tbinx = (ibin + bindx + Lcnum[0]) % Lcnum[0];
								int tbinIdx = tbinz * Lcnum[1] * Lcnum[0] +
									tbiny * Lcnum[0] + tbinx;

								int jcell = tbinIdx;
								for (int j = Head[jcell]; j != -1; j = List[j])  //loop over all particles in jcell
								{
									vec<DIM,double> vec_ij;
									double rij2, rcutrskin2;
									if (i >= j) continue;
									Minus(vec_ij, s[i].x, s[j].x);
									periodicVec(vec_ij, Length);
									//vec_ij = distance_vector(s[i].x, s[j].x, Length, invLength);
									rij2 = square(vec_ij);
									rcutrskin2 = pow(s[i].r + s[j].r + dskin, 2);
									if (rij2 < rcutrskin2)
									{
										++nlist[i];
										++nlist[j];
										Nnebrlist[INCREMENT*i + nlist[i]] = j;
										Nnebrlist[INCREMENT*j + nlist[j]] = i;
									}
								}
							}
						}
					}
					i = List[i];
				}
			}  //kbin end 
		}      //jbin end
	}          //ibin end 

	for (int i = 0; i < N; ++i) 
		Insertsort(&Nnebrlist[INCREMENT*i + 1], nlist[i]);
		
	IsValid = true;
	cntForce = 0;
}

bool NebrList::IsNebrlistInvalid(particle* s, double* Length, double* invLength, const Allparticle& Ap, const int& N)
{
	if (!IsValid) return true;
	for (int i = 0; i < N; ++i)
	{
		vec<DIM,double> vec_ij;
		Minus(vec_ij, s[i].x, sHold[i].x);
		double ri0;
		ri0 = scalarization(vec_ij);
		double dsig = s[i].r - sHold[i].r;

		if (ri0 + dsig >= dskin * 0.5) return true;
	}
	return false;
}

NebrList::~NebrList()
{
	delete[] List;
	delete[] Head;
	delete Celllist;
	if (IsList)
	{
		delete[] sHold;
		delete[] nlist;
		delete[] Nnebrlist;
	}
}

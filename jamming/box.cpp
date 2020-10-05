//#include"pch.h"
#include "box.h"
#include<iostream>
#include<fstream>
#include<iomanip>
#include<stdio.h>
#include <algorithm> 
#include<stdexcept>

#include<time.h>
using namespace std;

box::box(int input_N, bool monitor_i,const char* monitor_file_i) :
	N(input_N), 
	monitor(monitor_i)
{
	s = new particle[N];
	nebrlist = new NebrList;
}

void box::initialize()
{
	Length[0] = 1.0; Length[1] = 1.0; Length[2] = 1.0;
	Length[3] = 0.0; Length[4] = 0.0; Length[5] = 0.0;
	setBoxPara();
	epsilon = 1.0;

	MONITOR.updateSclInfo(s, N, boxVolume, Ap);
	for (int i = 0; i < N; ++i)
		s[i].r = 0.5*s[i].DiameterScale*Ap.meanDiameter;

	dt = Ap.meanDiameter / epsilon * 0.005;
	nebrlist->skinSet = 0.1 / 1.1;
	nebrlist->InitNebr(s, Length, invLength, Ap, N);
	MOVE.fire.GetFIREval(dt, N);
	MONITOR.updateSclInfo(s, N, boxVolume, Ap);
}


void box::inflate_volume(double target_volFrac)
{
	MONITOR.updateSclInfo(s, N, boxVolume, Ap);
	double factor = pow(target_volFrac / MONITOR.volFrac, 1.0 / DIM);
	Ap.meanDiameter = Ap.meanDiameter*factor;
	MONITOR.SphInfoInit = false;
	for (int i = 0; i < N; ++i)
		s[i].r = 0.5*s[i].DiameterScale*Ap.meanDiameter;
	MONITOR.updateSclInfo(s, N, boxVolume, Ap);
	dt = Ap.meanDiameter / epsilon * 0.005;
	nebrlist->InitNebr(s, Length, invLength, Ap, N);
	MOVE.fire.GetFIREval(dt, N);
}
void box::InitShear(double delg)
{
	SHEAR.InitShear(delg);
}
void box::XzShearStrain()
{
	SHEAR.XzShearStrain(s, N, Length);
	MOVE.fire.InitFIREval(dt, N);
	setBoxPara();
	nebrlist->InitNebr(s, Length, invLength, Ap, N);
}

void box::read_positions(const char* file,int ignore_line)
{
	std::ifstream infile;
	infile.open(file);
	if (!infile)
	{
		cout << "can't open this file" << endl;
		return;
	}

	for (int i = 0; i < ignore_line; i++)
		infile.ignore(256, '\n');
	char buff[20], c;
	infile.get(c); infile.get(c);
	infile.get(buff, 10, c);  infile.get(c);
	infile >> s[0].DiameterScale; infile.get(c);
	//infile.get(buff, 20, c);  infile.get(c);
	infile.get(buff, 10, c);  infile.get(c);

	for (int j = 0; j < DIM; ++j)
		infile >> s[0].x[j];
	infile.ignore(256, '\n');

	for (int i = 1; i < N; i++)
	{
		infile.get(buff, 10, c);  infile.get(c);
		infile.get(buff, 10, c);  infile.get(c);
		infile >> s[i].DiameterScale; infile.get(c);
		//infile.get(buff, 20, c);  infile.get(c);
		infile.get(buff, 10, c);  infile.get(c);

		for (int k = 0; k < DIM; k++)
			infile >> s[i].x[k];
		infile.ignore(256, '\n');
	}

	infile.close();
	Ap.meanDiameter = Ap.mean_Diameter(s, N);
	for (int i = 0; i < N; ++i)
		s[i].DiameterScale = s[i].DiameterScale / Ap.meanDiameter;
}
void box::record_positions(const char* file)
{
	std::ofstream out;
	char filename[256];
	sprintf_s(filename, "E:\\soft_sphere\\jamming\\jamming\\jamming_code\\FIRE_test\\final_positions\\0.40\\final_positions_500N_%d.dat", 1);
	out.open(filename);
	out << "potential_energy is: " << std::setprecision(16) << MONITOR.Uenergy << std::endl;
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < DIM; ++j)
		{
			out << s[i].x[j] << "     ";
		}
		out << std::endl;
	}
	out.close();
}

//FIRE Algorithm need to used with an integration algorithm
void box::FIRErelax()
{
	int steps = 0;
	bool stop = false;
	double old_e = -INFINITY;

	//FORCE.HarmonicForce(s, epsilon, N, Length, invLength,MONITOR);
	//FORCE.HarmonicForcewithCelllist(s, FORCE.zeroforce, epsilon, sigma_mean, N, Length, nebrlist);
	FORCE.HarmonicForcewithNebrlist(s, epsilon, Ap, N, Length, invLength, nebrlist, MONITOR);
	//fprintf_s(stdout, "particle->force: %.17f  %.17f  %.17f\n", s[486].f[0], s[416].f[1], s[258].f[2]);
	clock_t start0 = clock();
	while(!stop)
	{
		
		//std::cout << steps << ' ' << std::setprecision(30) << "e :" << MONITOR.Uenergy << std::endl;
		//std::cout << steps << ' ' << std::setprecision(19) << "dt :" << dt << std::endl;
		//std::cout << steps << ' ' << std::setprecision(19) << "mean_sigma :" << Ap.meanDiameter << std::endl;
		//std::cout << steps << ' ' << std::setprecision(19) << "diameter :" << (s[0].r*2.0) << std::endl;
		old_e = MONITOR.Uenergy;
		//=========================================================================================
		//   Initial velocity of the system is 0, so in the first circle, the Fire optimimization
		//   step can be neglected.We modify Vn using Vn-1 and Fn-1,then we update Xn use Vn.
		//   When P<0, Vn will be reset as 0, Xn+1=Xn,Vn+1 have the same orientation as Fn+1,the 
		//   gradient;
		//=========================================================================================
		MOVE.Euler_integration(s, N, dt);
		
		adjustImg();
		//clock_t start = clock();
		//FORCE.HarmonicForce(s, epsilon, N, Length, invLength, MONITOR);
		//FORCE.HarmonicForcewithCelllist(s, FORCE.zeroforce, epsilon, sigma_mean, N, Length, nebrlist);
		FORCE.HarmonicForcewithNebrlist(s, epsilon, Ap, N, Length, invLength, nebrlist, MONITOR);
		//clock_t end = clock();
		//double time = (double)(end - start)/CLK_TCK;
		//cout << "GetTickCount:" << time << endl;

		MOVE.fire.update_para(dt, N, s);
		++steps;
		if (FORCE.aveforce <= 1e-14)
			stop = true;
	}
	clock_t end0 = clock();
	double time0 = (double)(end0 - start0)/CLK_TCK;
	cout << "GetTickCount:" << time0 << endl;
	//MONITOR.UenergySys(s, N, epsilon, Length, invLength);
	std::cout << steps << ' ' << std::setprecision(16) << "e :" << MONITOR.Uenergy << std::endl;
	MONITOR.computerCoordination(s, N, epsilon, Length, invLength);
	std::cout << "aveCoordNumWithRattler: " << MONITOR.aveCoordNumWithRattler << std::endl;
	std::cout << "aveCoordNumNoRattler: " << MONITOR.aveCoordNumNoRattler << std::endl;
}

void box::show()
{
	std::cout << "number of particles are : " << N << std::endl;
	std::cout << "box length is : " << Length << std::endl;
	std::cout << "epsilon is : " << epsilon << std::endl;
	std::cout << "mean Diameter is : " << Ap.meanDiameter << std::endl;
	std::cout << "time step is : " << dt << std::endl;
}


//========================================================================
//====   Understanding Molecular Simulation                 ==============
//====                             --- Daan Frenkel         ==============
//========================================================================
//-4.2.3 Integrating the Equations of Motion
//-Algorithm 6
//The Radial Distribution Function
void box::gr(int SWITCH)
{
	double rho = N / boxVolume;
	if (SWITCH == 0)
	{
		nhis = NHIS;
		ngr = 0;
		delg = Length[0] / ((double)(2 * nhis));
		for (int i = 0; i < nhis; ++i)
			g[i] = 0;
	}

	//calculate RDF by setting each particles as original point, 
	//then calculate average value of RDF for all particles.
	else if (SWITCH == 1)
	{
		++ngr;
		double similen = Length[0] / 2;
		for (int i = 0; i < N; ++i)
			for (int j = i + 1; j < N; ++j)
			{
				vec<DIM,double> vec_ij;
				double rij;
				vec_ij = distance_vector(s[i].x, s[j].x, Length, invLength);
				rij = scalarization(vec_ij);
				if (rij < similen)
				{
					int ig = int(rij / delg);
					g[ig] += 2;
				}
			}
	}
	else if (SWITCH == 2)
	{
		for (int i = 0; i < nhis; ++i)
		{
			r[i] = delg * (i + 0.5);
			double vb = (pow((i + 1), 3) - pow(i, 3))*pow(delg, 3);
			double nid = 4 * (1.0 / 3.0)*PI*vb*rho;
			g[i] = g[i] / (ngr*N*nid);
		}
	}
}
void box::integrate_i(int i, double& sumv2)
{
	vec<DIM, double> xx, H;
	vecMulti(xx, 2.0, s[i].x);
	Minus(xx, xx, s[i].xm);
	vecMulti(H, dt*dt, s[i].f);
	Add(xx, xx, H);
	s[i].v = distance_vector(xx, s[i].xm, Length, invLength);
	vecDevi(s[i].v, s[i].v, 2 * dt);
	Add(sumv, sumv, s[i].v);
	double a = 0;
	vecDot(a, s[i].v, s[i].v);
	sumv2 = sumv2 + a;
	s[i].xm = s[i].x;
	s[i].x = xx;
}
void box::integrate()
{
	for (int i = 0; i < DIM; ++i)
		sumv[i] = 0;
	double sumv2 = 0;
	for (int i = 0; i < N; ++i)
		integrate_i(i, sumv2);
	tempreture = sumv2 / (3 * N);
	etot = (MONITOR.Uenergy + 0.5*sumv2) / N;
}
void box::verletinit()
{
	epsilon = 1.0;
	sprintf_s(buff, "etot%d.dat", 1);
	ofstream endfile;
	endfile.open(buff);
	endfile.close();
	for (int i = 0; i < N; ++i)
		for (int j = 0; j < DIM; ++j)
			s[i].xm[j] = s[i].x[j];
}
void box::NVEdsys()
{
	verletinit();
	gr(0);
	int nsample = 1;
	enum { INIT = 0, SAMPLE = 1, RESULT = 2 };
	Diffusion diffusion(s, dt, N);
	diffusion.dif(INIT, nsample, s);
	for (int steps = 1; steps <= 50; ++steps)
	{
		MONITOR.UenergySys(s, N, epsilon, Length, invLength);
		MONITOR.KenergySys(s, N, epsilon, Length);

		FORCE.HarmonicForcewithCelllist(s, FORCE.zeroforce, epsilon, Ap, N, Length, invLength,nebrlist);
		if (steps % 10 == 0)
			gr(1);
		integrate();
		if (steps%nsample == 0)
			diffusion.dif(SAMPLE, nsample, s);
		//for (int i = 0; i < N; ++i)
			//PBC(i);
		cout << "etot: " << MONITOR.Kenergy << endl;
		ofstream endfile;
		endfile.open(buff, std::ios_base::app);
		//endfile << steps << ' ' << sumv[0] << ' ' << sumv[1] << ' ' << sumv[2] << endl;
		endfile.close();
	}
	gr(2);
	diffusion.dif(RESULT, nsample, s);
	diffusion.show();
	for (int i = 0; i < nhis; ++i)
		cout << r[i] << ' ' << g[i] << endl;
}
void box::NVTsys()
{
	gr(0);
	int nsample = 1;
	enum { INIT = 0, SAMPLE = 1, RESULT = 2 };
	Diffusion diffusion(s, dt, N);
	double T = 5.0;
	NVT NVTensem(T, N);
	diffusion.dif(INIT, nsample, s);
	for (int steps = 1; steps <= 500; ++steps)
	{
		MONITOR.UenergySys(s, N, epsilon, Length, invLength);
		//MONITOR.KenergySys(s, N, epsilon, Length);

		NVTensem.Chain(s, MONITOR.Kenergy, dt, N);
		FORCE.HarmonicForcewithCelllist(s, FORCE.zeroforce, epsilon, Ap, N, Length, invLength, nebrlist);
		NVTensem.Pos_Vel(s, MONITOR.Kenergy, dt, N);
		if (steps % 10 == 0)
			gr(1);
		if (steps%nsample == 0)
			diffusion.dif(SAMPLE, nsample, s);
		//for (int i = 0; i < N; ++i)
			//PBC(i);
		MONITOR.TemperaSys(s, N, epsilon, Length);
		cout << "Temperature: " << MONITOR.Temperature << endl;
		ofstream endfile;
		endfile.open(buff, std::ios_base::app);
		//endfile << steps << ' ' << sumv[0] << ' ' << sumv[1] << ' ' << sumv[2] << endl;
		endfile.close();
	}
	gr(2);
	diffusion.dif(RESULT, nsample, s);
	diffusion.show();
	for (int i = 0; i < nhis; ++i)
		cout << r[i] << ' ' << g[i] << endl;
}


void box::setBoxPara()
{
	double negHalf[3] = { 0.5,0.5,0.5 };
	vecHvoigtMulVec(boxLo, Length, negHalf);
	double tmp3[3] ={ Length[0],Length[1],Length[2] };
	vecAdd(boxHi, tmp3, boxLo);

	boxVolume = Length[0] * Length[1] * Length[2];

	invLength[0] = 1.0 / Length[0];
	invLength[1] = 1.0 / Length[1];
	invLength[2] = 1.0 / Length[2];

	invLength[3] =
		-Length[3] / (Length[1] * Length[2]);
	invLength[4] =
		Length[3] * Length[5] /
		(Length[0] * Length[1] * Length[2]) -
		Length[4] / (Length[0] * Length[2]);
	invLength[5] =
		-Length[5] / (Length[0] * Length[1]);
}
void box::adjustImg()
{
	for (int iatom = 0; iatom < N; iatom++) {
		particle *posXyz = &s[iatom];

		//Length is elements in L-tensor.If we want know which image box i'th atom is in, we can simple use formula L^-1*X ,X is 
		//positions of i'th atom.This formula is also adapted to distance between two particles.
		double lamda_x = invLength[0] * posXyz->x[0] +
			invLength[5] * posXyz->x[1] +
			invLength[4] * posXyz->x[2] + 0.5;
		double lamda_y = invLength[1] * posXyz->x[1] +
			invLength[3] * posXyz->x[2] + 0.5;
		double lamda_z = invLength[2] * posXyz->x[2] + 0.5;
		double shift_x = floor(lamda_x), shift_y = floor(lamda_y),
			shift_z = floor(lamda_z);

		double delta_x = Length[0] * shift_x + Length[5] * shift_y +
			Length[4] * shift_z;
		double delta_y = Length[1] * shift_y + Length[3] * shift_z;
		double delta_z = Length[2] * shift_z;
		posXyz->x[0] -= delta_x;
		posXyz->x[1] -= delta_y;
		posXyz->x[2] -= delta_z;
	}
}
void box::setBasicInfo() 
{
	MONITOR.updateSclInfo(s, N, boxVolume, Ap);
}
box::~box()
{
	delete[] s;
	delete nebrlist;
}

//#include "pch.h"
#include "FIRE.h"
#include<iostream>
#include<iomanip>
#include<omp.h> 
using namespace std;


FIRE::FIRE(const char* file,int ignore_line) : 
	box()
{
	initialize();

	CurrentStep = 0;
	NegativeStep = 0;

	P = 0;

	dtmax = 100 * dt;
	finc = 1.1;
	Nmin = 5;
	fdec = 0.5;
	alpha0 = 0.1;
	falpha = 0.99;
	alpha = alpha0;

	for (int i = 0; i < N*DIM; ++i)
	{
		F_norm.push_back(0);
		V.push_back(0);
	}
}

FIRE::FIRE(const char* file, int ignore_line, int N_i,double box_length) :
	box(N_i ,box_length)
{
	initialize();

	CurrentStep = 0;
	NegativeStep = 0;

	P = 0;

	dtmax = 100 * dt;
	finc = 1.1;
	Nmin = 5;
	fdec = 0.5;
	alpha0 = 0.1;
	falpha = 0.99;
	alpha = alpha0;

	for (int i = 0; i < N*DIM; ++i)
	{
		F_norm.push_back(0);
		V.push_back(0);
	}
}

FIRE::FIRE(const char* file, int ignore_line, int N_i, double box_length, double phi_i, int N_min) :
	box(N_i, box_length),
	phi(phi_i),
	Nmin(N_min)
{
	initialize();

	CurrentStep = 0;
	NegativeStep = 0;

	P = 0;

	dtmax = 100 * dt;
	finc = 1.1;
	fdec=0.5;
	alpha0 = 0.1;
	falpha = 0.99;
	alpha = alpha0;

	for (int i = 0; i < N*DIM; ++i)
	{
		F_norm.push_back(0);
		V.push_back(0);
	}
}

void FIRE::FIREinitial(char* file,int ignore_line)
{
	CurrentStep = 0;
	NegativeStep = 0;

	alpha = alpha0;

	for (int i = 0; i < N; ++i)
		for (int j = 0; j < DIM; ++j)
			s[i].v[j] = 0;
	read_positions(file, ignore_line);
}

double FIRE::P_value()
{
	double p0 = 0;
	p0 = scalar_product(V, F_norm);
	P = p0;
	return p0;
}

void FIRE::update_para()
{
	if (P > 0)
	{
		if (CurrentStep - NegativeStep > Nmin)
		{
			dt *= finc;
			dt = (dt < dtmax ? dt : dtmax);
			alpha *= falpha;
		}

		double VV = 0, FF = 0;
//#pragma omp parallel for private(j)
		VV = scalar_product(V, V);
		FF = scalar_product(F_norm, F_norm);
//#pragma omp barrier

		double VF = sqrt(VV / FF);
//#pragma omp parallel for
		for (int i = 0; i < DIM*N; ++i)
			V[i] = (1 - alpha)*V[i] + alpha * F_norm[i] * VF;
	}
	else
	{
		dt *= fdec;
		NegativeStep = CurrentStep;
		alpha = alpha0;
		for (int i = 0; i < DIM*N; ++i)
			V[i] = 0;
	}
	for (int i = 0; i < N; ++i)
		for (int j = 0; j < DIM; ++j)
			s[i].v[j] = V[i*DIM + j];
}

void FIRE::FIRErelax(double phi_i)
{
	phi = phi_i;
	inflate_volume(phi);
	dtmax = 100 * dt;

	//cell list
	//cellmolloc();
	//cellListBuild();
	//initcellnebr();
	//verlet list
	//for (int i = 0; i < N; ++i)
		//updateList(i);

	bool stop = false;
	double old_e = -INFINITY;

	zeroforceratio = 0;
//#pragma omp parallel for 
	for (int i = 0; i < N; ++i)
		s[i].f = force_i(i);

	zeroforceratio /= N;

	double start = omp_get_wtime();//获取起始时间 
	while (!stop) {

		//cout << CurrentStep << endl;
		system_property();
		cout << std::setprecision(16) << "e :" << potential_energy << endl;
		if (potential_energy / N < 1e-16 || fabs(potential_energy - old_e) / N < 1e-15)
			if (potential_energy != old_e)
				stop = true;
		old_e = potential_energy;
		//=========================================================================================
		//   Initial velocity of the system is 0, so in the first circle, the Fire optimimization
		//   step can be neglected.We modify Vn using Vn-1 and Fn-1,then we update Xn use Vn.
		//   When P<0, Vn will be reset as 0, Xn+1=Xn,Vn+1 have the same orientation as Fn+1,the 
		//   gradient;
		//=========================================================================================
		//the order of these three function need to be noticed
		for (int i = 0; i < N; ++i)
			Euler_integration(i);
			
		//cellListBuild();
			

		zeroforceratio = 0;
//#pragma omp parallel for 
		for (int i = 0; i < N; ++i)
			s[i].f = force_i(i);


		zeroforceratio /= N;

		for (int i = 0; i < N; ++i)
			for (int j = 0; j < DIM; ++j)
			{
				V[i*DIM + j] = s[i].v[j];
				F_norm[i*DIM + j] = s[i].f[j];
			}
		P_value();
		update_para();

		++CurrentStep;
	}

	//for (int i = 0; i < N; ++i)
		//PBC(i); 
	double volume = length * length*length;
	double fac = 0;
	for (int i = 0; i < N; ++i) {
		fac += pow(s[i].r, 3.0);
	}
	shearPhi = 4 * PI / 3 * fac / volume;
	std::cout << "shearPhi: " << shearPhi << std::endl;
	cout << CurrentStep << endl;
	std::cout << "potential_energy is: " << old_e << std::endl;
	double end = omp_get_wtime();
	cout << "计算耗时为：" << end - start;
	cout << endl;

	//cout << CurrentStep << endl;
	//system_property();
	cout << endl;
	celldelete();
}

void FIRE::compress(double dphi_i, double dphi_f, double phi_factor)
{
	double dphi = dphi_i;

	//system_property();
	double phi_i = phi;

	while (dphi < dphi_f)
	{
		inflate_volume(phi_i + dphi); //instantneously inflate spheres to the target density
		//Euler_integration();
		force();
		P_value();
		update_para();                                  //use FIRE to relax the energy
		//system_property();                         //update global variables
		dphi *= phi_factor;
		++CurrentStep;
	}
}

void FIRE::FIREcompress(double phi_i)
{
	bool start = false;
	double old_e = -INFINITY;
	phi = phi_i;
	inflate_volume(phi);
	force();
	while (zeroforceratio < 1.0) {

		cout << CurrentStep << endl;
		cout << endl;

		//=========================================================================================
		//   Initial velocity of the system is 0, so in the first circle, the Fire optimimization
		//   step can be neglected.We modify Vn using Vn-1 and Fn-1,then we update Xn use Vn.
		//   When P<0, Vn will be reset as 0, Xn+1=Xn,Vn+1 have the same orientation as Fn+1,the 
		//   gradient;
		//=========================================================================================
		//the order of these three function need to be noticed
		//Euler_integration();
		force();
		P_value();
		update_para();
		++CurrentStep;

		system_property();
		if (potential_energy / N < 1e-16 || (fabs(potential_energy - old_e) / N < 1e-15))
			start = true;
		old_e = potential_energy;
		while (start)
		{
			const double dphi_i = 1e-6;
			const double dphi_f = 1e-2;
			const double phi_factor = 1.1;
			std::cout << "Compressing ..." << std::endl;
			std::cout << std::endl;

			compress(dphi_i, dphi_f, phi_factor);
			start = false;

			inflate_volume(phi_i);
		}
	}
	cout << CurrentStep << endl;
	cout << "potential is " << potential_energy << endl;
	cout << "zeroforceratio is " << zeroforceratio << endl;
	cout << endl;
}

void FIRE::test()
{
	//cout << "====Initial Position of Particles====" << endl;
	//for (int i = 0; i < N; ++i)
	//{
	//	for (int j = 0; j < DIM; ++j)
	//		cout << std::setprecision(16) << s[i].x[j] << "   ";
	//	cout << endl;
	//}
	cout << "==========Parameters of FIRE===========" << endl;
	cout << "Currstep: " << CurrentStep << endl;
	cout << "NegativeStep: " << NegativeStep << endl;
	cout << "P: " << P << endl;
	cout << "dt: " << dt << endl;
	cout << "dtmax: " << dtmax << endl;
	cout << "finc: " << finc << endl;
	cout << "falpha: " << falpha << endl;
	cout << "fdec: " << fdec << endl;
	cout << "Nmin: " << Nmin << endl;
	cout << "alpha0: " << alpha0 << endl;
}

FIRE::~FIRE()
{
}

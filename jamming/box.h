#pragma once
#include"Device.h"
#include"NebrList.h"

#include"read_input.h"
#include<cmath>
#define FILE_BUFFER_LENGTH 30000
#define NHIS 20

class box
{
	friend read_input;
	friend class NebrList;
	friend class Diffusion;
	friend int main(int argc, char* argv[]);
public:
	box(int N_i, bool monitor_i, const char* monitor_file_i);
	void initialize();
	void setBasicInfo();
	void inflate_volume(double target_volFrac);
	void InitShear(double delg);
	void XzShearStrain();
	void read_positions(const char* file, int ignore_line);
	void record_positions(const char* file);

	Monitor MONITOR;  //monitor Energy,Uenergy,Kenergy.
	Move MOVE;        //move particle.
	Force FORCE;      //measure force between two particles.
	void FIRErelax();
	Shear SHEAR;
	Allparticle Ap;
	void NVTsys();

	void show();
	void setBoxPara();
	void adjustImg();
	//========================================================================
    //====   Understanding Molecular Simulation                 ==============
    //====                             --- Daan Frenkel         ==============
    //========================================================================
    //-4.2.3 Integrating the Equations of Motion
    //-Algorithm 6
	void verletinit();
	void integrate_i(int i, double& sumv2);
	void integrate();
	double tempreture;
	double etot;           //average total energy per particle
	vec<DIM,double> sumv;
	void NVEdsys();
	char buff[300];

	//-4.4 Computer Experiment
	//-Algorithm7
	//The Radial Distribution Function
	void gr(int SWITCH);
	int ngr;
	double delg;
	int nhis;
	double g[NHIS];
	double r[NHIS];

	virtual ~box();
	
protected:
	particle *s;
	int N;
	bool monitor;

	double dt;
	double Temperature;
	double r_verlet;

	NebrList *nebrlist;
private:
	//box parameter
	double Length[Hvoigt6];     //{[Length(0),Length(5),Length(4)],        {[Lx ,Lyx,Lzx],
								// [0,        Length(1),Length(3)],   =     [0  ,Ly ,Lzy],
								// [0,       ,0        ,Length(2)]}         [0  ,0  ,Lz ]}
	double boxVolume;
	double epsilon;

	double boxLo[DIM];
	double boxHi[DIM];
	double invLength[Hvoigt6];     //{[(1/Lx),(-Lyx/(Lx*Ly)),((Lyx*Lzy)/(Lx*Ly*Lz)-(Lzx/(Lx*Lz))],
								   // [0     ,(1/Ly)        ,(-Lzy/(Ly*Lz))                     ],
								   // [0     ,0             ,(1/Lz)                             ]}
};

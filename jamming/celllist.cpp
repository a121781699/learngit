#include"pch.h"
#include"box.h"
#include"Chain.h"
#include<iostream>
using namespace std;

Chain<particle>* cList = NULL;       //it present cell using a List recording information in the cell
particle* ptr = NULL;                //virtural particle in each cell,it present the head of cList
particle* s = new particle[50];      //real simulation particle

void cellListBuild()
{
	box b;
	int n = (int)(1.0 / 0.22);             //       n=(int)b.length/sigma
	cList = new Chain<particle>[n*n*n];
	ptr = new particle[n*n*n];

	double rn = 1.0 / (double)n;          //       rn=b.length/(double)n

	//set head of cell using virtrual particle
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			for (int k = 0; k < n; ++k)
			{
				double icellx = (double)k * rn;
				double icelly = (double)j * rn;
				double icellz = (double)i * rn;
				int inumber = k + j * n + i * n*n;

				ptr[inumber].x[1] = icellx;
				ptr[inumber].x[2] = icelly;
				ptr[inumber].x[3] = icellz;
				cList[inumber].Append(ptr[inumber]);
			}
		}
	}

	//estimate which cell particle i belong to
	for (int i = 0; i < 50; ++i)
	{
		int celli = (int)((s[i].x[0] + 0.5) / rn);
		int cellj = (int)((s[i].x[1] + 0.5) / rn);
		int cellk = (int)((s[i].x[2] + 0.5) / rn);
		int inumber = celli + cellj * n + cellk * n*n;
		cList[inumber].Append(s[i]);
	}


	delete[] cList;
	delete[] ptr;
}

int main()
{
	cellListBuild();
	delete[] s;
}
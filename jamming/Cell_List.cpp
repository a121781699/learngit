#include"pch.h"
#include"box.h"
#include"Chain.h"
#include<iostream>
#define N 500
#define neighboor 27
using namespace std;

box b;

int List[N];
particle* s = new particle[N];
int n = (int)(1.0 / 0.22);             //       n=(int)(b.length/sigma)     n is number of of cell of each direction(x or y or z)
int cnum = n * n*n;
double rn = 1.0 / (double)n;          //       rn=b.length/(double)n      rn is length of one cell
int *Head = new int[cnum];
Chain<int> *nebrlist = new Chain<int>[cnum];

void cellListBuild()
{
	//initialize Head
	for (int i = 0; i < N; ++i)
	{
		Head[i] = -1;
	}
	//set List(put index of particles into List)
	//set Head(reserve the index of last particle of List)
	for (int i = 0; i < N; ++i)
	{
		int celli = (int)((s[i].x[0] + 0.5) / rn);
		int cellj = (int)((s[i].x[1] + 0.5) / rn);
		int cellk = (int)((s[i].x[2] + 0.5) / rn);
		int inumber = celli + cellj * n + cellk * n*n;

		List[i] = Head[inumber];
		Head[inumber] = i;
	}
}

void initcellnebr(int icell, int ncell)
{
	box b;
	//postions of cells
	MyArray<double>* x = new MyArray<double>[cnum];
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			for (int k = 0; k < n; ++k)
			{
				int inum = i * n*n + j * n + k;
				x[inum][0] = k * rn;
				x[inum][1] = j * rn;
				x[inum][2] = i * rn;
			}
		}
	}
	MyArray<double> vec;
	double Rdia = sqrt(2 * rn*rn);
	for (int inum = 0; inum < n*n*n; ++inum)
	{
		for (int jnum = 0; jnum < n*n*n; ++jnum)
		{
			if (inum == jnum) continue;
			vec = b.distance_vector(x[inum], x[jnum]);
			double rij = vec.scalarization();
			if (rij <= Rdia)
			{
				nebrlist[inum].Append[jnum];
			}
		}
	}
	delete[] x;
}

int findneigh(int icell, int ncell)
{
	ChainNode<int> *ptemp = nebrlist[icell].first;
	for (int i = 0; i < ncell;++i)
		ptemp = ptemp->link;
	return *ptemp->data;
}

MyArray<double> forceij(int i, int j);
MyArray<double> force_i(int i);

int main()
{
	cellListBuild();
	delete[] Head;
	delete[] s;
}

MyArray<double> force_i(int i)
{
	int jcell = 0;
	MyArray<double> f;
	int celli = (int)((s[i].x[0] + 0.5) / rn);
	int cellj = (int)((s[i].x[1] + 0.5) / rn);
	int cellk = (int)((s[i].x[2] + 0.5) / rn);
	int icell = celli + cellj * n + cellk * n*n;
	for (int inebr = 0; inebr < neighboor; ++inebr)
	{
		jcell = findneigh(icell, inebr);
		int j = Head[jcell];
		while (j != -1)
		{
			if (i != j)
				f = f + forceij(i, j);
			j = List[j];
		}
	}
}
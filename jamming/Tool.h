#pragma once
#include"Allparticle.h"
#include"Shear.h"
#include<vector>

vec<DIM,double> distance_vector(vec<DIM,double>& x, vec<DIM,double>& y, double* Length,double* invLength);
void distance_vector(vec<DIM, double>& x, double* Length, double* invLength);
void periodicVec(vec<DIM, double>& Rij, double* Length);
double scalar_product(std::vector<double> const& a, std::vector<double> const& b);
vec<DIM,double> OrthoDistance(vec<DIM,double>& x, vec<DIM,double>& y, double* Length);

void vecHvoigtMulVec(double* vHmV, double* H, double* vecA);
void vecAdd(double* vSum, double* vecA, double* vecB);

template<class T>
void insert(T a[], int n, const T& x)
{
	int i;
	for (i = n - 1; i >= 0 && x < a[i]; --i)
		a[i + 1] = a[i];
	a[i + 1] = x;
}

template<class T>
void Insertsort(T a[], int n)
{
	for (int i = 1; i < n; ++i)
	{
		T t = a[i];
		insert(a, i, t);
	}
}


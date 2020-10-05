#include"pch.h"
#include"Tool.h"
#include <numeric>

double scalar_product(std::vector<double> const& a, std::vector<double> const& b)
{
	if (a.size() != b.size()) { throw std::runtime_error("different sizes"); }
	return std::inner_product(a.begin(), a.end(), b.begin(), 0.0);
} // scalar_product

vec<DIM,double> distance_vector(vec<DIM,double>& x, vec<DIM,double>& y, double* Length,double* invLength)
{
	vec<DIM,double> vec;
	Minus(vec, x, y);
	//(L^-1)*Vec get integer torwarding bottom,presented by [(L^-1)*Vec],[(L^-1)*Vec + 0.5*I]
	//which represent i'th cell(using 3D vector)where B sphere is(the origin is A sphere).Vec=B-A.
	//L*[(L^-1)*Vec] represent the postions of centre of i'th cell.
	double lamda_x = invLength[0] * vec[0] +
					 invLength[5] * vec[1] +
					 invLength[4] * vec[2] + 0.5;
	double lamda_y = invLength[1] * vec[1] +
					 invLength[3] * vec[2] + 0.5;
	double lamda_z = invLength[2] * vec[2] + 0.5;
	double shift_x = floor(lamda_x), shift_y = floor(lamda_y),
		shift_z = floor(lamda_z);

	double delta_x = Length[0] * shift_x + Length[5] * shift_y +
					 Length[4] * shift_z;
	double delta_y = Length[1] * shift_y + Length[3] * shift_z;
	double delta_z = Length[2] * shift_z;
	vec[0] -= delta_x;
	vec[1] -= delta_y;
	vec[2] -= delta_z;

	return vec;
}
void distance_vector(vec<DIM, double>& x, double* Length, double* invLength)
{
	//(L^-1)*Vec get integer torwarding bottom,presented by [(L^-1)*Vec],[(L^-1)*Vec + 0.5*I]
	//which represent i'th cell(using 3D vector)where B sphere is(the origin is A sphere).Vec=B-A.
	//L*[(L^-1)*Vec] represent the postions of centre of i'th cell.
	double lamda_x = invLength[0] * x[0] +
		invLength[5] * x[1] +
		invLength[4] * x[2] + 0.5;
	double lamda_y = invLength[1] * x[1] +
		invLength[3] * x[2] + 0.5;
	double lamda_z = invLength[2] * x[2] + 0.5;
	double shift_x = floor(lamda_x), shift_y = floor(lamda_y),
		shift_z = floor(lamda_z);

	double delta_x = Length[0] * shift_x + Length[5] * shift_y +
		Length[4] * shift_z;
	double delta_y = Length[1] * shift_y + Length[3] * shift_z;
	double delta_z = Length[2] * shift_z;
	x[0] -= delta_x;
	x[1] -= delta_y;
	x[2] -= delta_z;
}
void periodicVec(vec<DIM, double>& Rij, double* Length)
{
	if (Rij[2] > 0.5*Length[2])
	{
		Rij[2] -= Length[2];
		Rij[0] -= Length[4];
	}
	else if (Rij[2] < -0.5*Length[2])
	{
		Rij[2] += Length[2];
		Rij[0] += Length[4];
	}
	if (Rij[1] > 0.5*Length[1])
	{
		Rij[1] -= Length[1];
	}
	else if (Rij[1] < -0.5*Length[1])
	{
		Rij[1] += Length[1];
	}
	if (Rij[0] > 0.5*Length[0])
	{
		Rij[0] -= Length[0];
	}
	else if (Rij[0] < -0.5*Length[0])
	{
		Rij[0] += Length[0];
	}
}

vec<DIM,double> OrthoDistance(vec<DIM,double>& x, vec<DIM,double>& y, double* Length)
{
	int n;
	vec<DIM,double> vec;
	Minus(vec, x, y);
	for (int k = 0; k < DIM; ++k)
	{
		n = (int)(vec[k] * 1.0 / Length[k] + ((vec[k] >= 0.0) ? 0.5 : -0.5));   //a fast algorithm to compute the minimum image distance
		vec[k] -= n * Length[k];                                                   //-box_size/2.0 < vec_ij[k] < box_size/2.0
	}
	return vec;
}


void vecHvoigtMulVec(double* vHmV, double* H, double* vecA)
{
	vHmV[0] = H[0] * vecA[0] + H[5] * vecA[1] + H[4] * vecA[2]; 
	vHmV[1] = H[1] * vecA[1] + H[3] * vecA[2];                 
	vHmV[2] = H[2] * vecA[2];
}

void vecAdd(double* vSum, double* vecA, double* vecB)
{
	vSum[0] = vecA[0] + vecB[0];
	vSum[1] = vecA[1] + vecB[1];
	vSum[2] = vecA[2] + vecB[2];
}
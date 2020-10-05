#pragma once
//#include<cmath>
#define DIM 3
#define Hvoigt6 6
template<unsigned int D, class T>
class vec {
public:
	T Arr[D];
	
public:
	//construct
	vec() = default;

	//operator function
	T& operator[](const unsigned int index);

	//destructor
	~vec() = default;
};
template<unsigned int D, class T>
T& vec<D, T>::operator[](const unsigned int index) {
	return Arr[index];
}

template<unsigned int D, class T>
double scalarization(const vec<D,T> &arr)
{
	double r_square = 0;
	for (int i = 0; i < D; ++i)
		r_square += arr.Arr[i] * arr.Arr[i];
	double r = sqrt(r_square);
	return r;
}

template<unsigned int D, class T>
double square(const vec<D,T> &arr)
{
	double r_square = 0;
	for (int i = 0; i < DIM; ++i)
		r_square += arr.Arr[i] * arr.Arr[i];
	return r_square;
}

template<unsigned int D, class T>
void Add(vec<D,T> &arr0, const vec<D,T> &arr1, const vec<D,T> &arr2) {
	for (int i = 0; i < DIM; ++i) {
		arr0.Arr[i] = arr1.Arr[i] + arr2.Arr[i];
	}
}

template<unsigned int D, class T>
void Minus(vec<D,T> &arr0, const vec<D,T> &arr1, const vec<D,T> &arr2) {
	for (int i = 0; i < DIM; ++i) {
		arr0.Arr[i] = arr1.Arr[i] - arr2.Arr[i];
	}
}

template<unsigned int D, class T>
void vecDot(double &A, const vec<D,T> &arr1, const vec<D,T> &arr2) {
	A = arr1.Arr[0] * arr2.Arr[0] + arr1.Arr[1] * arr2.Arr[1] + arr1.Arr[2] * arr2.Arr[2];
}

template<unsigned int D, class T,class B>
void vecMulti(vec<D,T> &arr0, const B &A, const vec<D,T> &arr1) {
	for (int i = 0; i < DIM; ++i) {
		arr0.Arr[i] = arr1.Arr[i] * A;
	}
}

template<unsigned int D, class T>
void vecDevi(vec<D,T> &arr0, const vec<D,T> &arr1, const T &A){
	for (int i = 0; i < DIM; ++i) {
		arr0.Arr[i] = arr1.Arr[i] / A;
	}
}

template<unsigned int D, class T>
void vecScaleAdd(vec<D,T> &arr0, const vec<D,T> &arr1, const T &A, const vec<D,T> &arr2)
{
	for (unsigned int i = 0; i < DIM; ++i)
		arr0.Arr[i] = arr1.Arr[i] + A * arr2.Arr[i];
}
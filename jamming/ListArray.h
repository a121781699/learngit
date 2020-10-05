#pragma once
#include<iostream>
using namespace std;
#define ELEMENT 270000
#define LISTNUM 10000
template<class T>
class ListArray
{
	friend class box;
	friend class Force;
	friend class NebrList;
private:
	T Array[ELEMENT + 2];
	T* First[LISTNUM + 2];
	T* Last[LISTNUM + 2];
	//int List_Index;
public:
	ListArray();
	//int GetIndex() const{ return List_Index; }
	//void SetIndex(int index) { List_Index = index; }
	ListArray<T>& Move_Left(int i);  //i'th List move left
	ListArray<T>& Move_Right(int i);  //i'th List move right
	int IsEmpty_Left(int i);   //if have enough space between i'th and (i+1)'th list
	int IsEmpty_Right(int i);
	ListArray<T>& Insert(int i, int k, const T x);      //insert element x behind k'th in i'th list 
	T& operator[](int index);
	void output(ostream& out) const;
	~ListArray();
};
template<class T>
ListArray<T>::ListArray()
{
	First[0] = Last[0] = &Array[0];
	First[LISTNUM + 1] = Last[LISTNUM + 1] = &Array[ELEMENT + 1];
	for (int i = 1; i <= LISTNUM; ++i) {
		First[i] = &Array[i];
		Last[i] = &Array[i];
	}
}
template<class T>
ListArray<T>::~ListArray()
{}
template<class T>
ListArray<T>& ListArray<T>::Move_Left(int i)
{
	if (i == 0)
		return *this;
	for (int* k = First[i]; k < Last[i]; ++k)
		*(k - 1) = *k;
	--First[i];
	--Last[i];
	return *this;
}

template<class T>
ListArray<T>& ListArray<T>::Move_Right(int i)
{
	if (i == LISTNUM + 1)
		return *this;
	for (int* k = Last[i] - 1; k >= First[i]; --k)
		*(k + 1) = *k;
	++First[i];
	++Last[i];
	return *this;
}

template<class T>
T& ListArray<T>::operator[](int index)
{
	if (index<1 || index>ELEMENT) throw - 1;
	return Array[index];
}

template<class T>
int ListArray<T>::IsEmpty_Left(int i)
{
	int n = 1;
	for (i; i > 0; --i)
	{
		if (First[i] > Last[i - 1])
			return n;
		++n;
	}
	return 0;
}

template<class T>
int ListArray<T>::IsEmpty_Right(int i)
{
	int n = 1;
	for (i; i < LISTNUM + 1; ++i)
	{
		if (Last[i] < First[i + 1])
			return n;
		++n;
	}
	return 0;
}

template<class T>
ListArray<T>& ListArray<T>::Insert(int i, int k, const T x)
{
	if (i<1 || i>LISTNUM) throw - 1;
	int listsize = Last[i] - First[i];
	if (k<0 || k>listsize) throw - 1;
	int inum = IsEmpty_Left(i);
	int jnum = IsEmpty_Right(i);
	if (inum == 0 || (inum >= jnum && jnum != 0))
	{
		if (jnum == 0)
			throw - 1;
		for (jnum; jnum > 1; --jnum)
			Move_Right(i + jnum - 1);
		for (int l = listsize; l > k; --l)
			*(First[i] + l) = *(First[i] + l - 1);
		*(First[i] + k) = x;
		++Last[i];
		return *this;
	}
	else
	{
		for (inum; inum > 1; --inum)
			Move_Left(i - inum + 1);
		for (int l = 0; l < k; ++l)
			*(First[i] + l - 1) = *(First[i] + l);
		*(First[i] + k - 1) = x;
		--First[i];
		return *this;
	}
	throw - 1;
}

template<class T>
void ListArray<T>::output(ostream& out) const
{
	for (int i = 1; i <= LISTNUM; ++i)
	{
		if (First[i] != Last[i])
		{
			out << i << ": " << "   ";
			for (int* j = First[i]; j < Last[i]; ++j)
				out << *j << "   ";
			out << "\n";
		}
	}
}

template<class T>
ostream& operator<<(ostream& out, const ListArray<T> x)
{
	x.output(out); return out;
}
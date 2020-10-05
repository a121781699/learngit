#pragma once
#pragma once
#include<iostream>
#include"ChainNode.h"
using namespace std;

template<class T>
class Chain
{
	template<class T> friend class ChainIterator;
	friend class box;
public:
	Chain() { first = 0; }
	~Chain();
	bool IsEmpty() const { return first == 0; }
	int Length() const;
	bool Find(int k, T& x) const;
	int Search(const T& x) const;
	Chain<T>& Delete(int k, T& x);
	Chain<T>& Insert(int k, const T x);
	Chain<T>& Append(T& x);
	void Output(ostream& out) const;
	void Erase();
private:
	ChainNode<T> *first; // ָ���һ���ڵ��ָ��
	ChainNode<T> *last;
};

template<class T>
Chain<T>::~Chain()
{//�������������������ɾ�������е����нڵ�
	Erase();
}


template<class T>
int Chain<T>::Length() const
{// ���������е�Ԫ������
	ChainNode<T> *current = first;
	int len = 0;
	while (current) {
		len++;
		current = current->link;
	}
	return len;
}

template<class T>
bool Chain<T>::Find(int k, T& x) const
{// Ѱ�������еĵ�k��Ԫ�أ������䴫����x
	//��������ڵ�k��Ԫ�أ��򷵻�f a l s e�����򷵻�t r u e
	if (k < 1) return false;
	ChainNode<T> *current = first;
	int index = 1; // current������
	while (index < k && current) {
		current = current->link;
		index++;
	}
	if (current) { x = current->data; return true; }
	return false; // �����ڵ�k��Ԫ��
}

template<class T>
int Chain<T>::Search(const T& x) const
{// Ѱ��x���������x���򷵻�x�ĵ�ַ
 //���x���������У��򷵻�0
	ChainNode<T> *current = first;
	int index = 1; // current������
	while (current && current->data != x) {
		current = current->link;
		index++;
	}
	if (current) return index;
	return 0;
}

inline int OutOfBounds()
{
	cerr << "size of k out of bounds!" << endl;
	return -1;
}

template<class T>
Chain<T>& Chain<T>::Delete(int k, T& x)
{// �ѵ�k��Ԫ��ȡ��x��Ȼ���������ɾ����k��Ԫ��
	//��������ڵ�k��Ԫ�أ��������쳣 OutOfBounds
	if (k < 1 || !first)
		throw OutOfBounds(); // �����ڵ�k��Ԫ��
			// p���ս�ָ���k���ڵ�
	ChainNode<T> *p = first;
	// ��p�ƶ�����k��Ԫ�أ�����������ɾ����Ԫ��
	if (k == 1) // p�Ѿ�ָ���k��Ԫ��
		first = first->link; // ɾ��֮
	else {// ��qָ���k - 1��Ԫ��
		ChainNode<T> *q = first;
		for (int index = 1; index < k - 1 && q; index++)
			q = q->link;
		if (!q || !q->link)
			throw OutOfBounds(); //�����ڵ�k��Ԫ��
		p = q->link; // ���ڵ�k��Ԫ��
		if (p == last) last = q;
		q->link = p->link;// ��������ɾ����Ԫ��
	}		//�����k��Ԫ�ز��ͷŽڵ�p
	x = p->data;
	delete p;
	return *this;
}

template<class T>
Chain<T>& Chain<T>::Insert(int k, const T x)
{// �ڵ�k��Ԫ��֮�����x
	//��������ڵ�k��Ԫ�أ��������쳣O u t O f B o u n d s
	// ���û���㹻�Ŀռ䣬�򴫵�N o M e m�쳣
	if (k < 0) throw OutOfBounds();
	// p���ս�ָ���k���ڵ�
	ChainNode<T> *p = first;
	//��p�ƶ�����k��Ԫ��
	for (int index = 1; index < k && p; index++)
		p = p->link;
	if (k > 0 && !p) throw OutOfBounds(); //�����ڵ�k��Ԫ��
	// ����
	ChainNode<T> *y = new ChainNode<T>;
	y->data = x;
	if (k > 1) {// ��p֮�����
		y->link = p->link;
		p->link = y;
	}
	else {// ��Ϊ��һ��Ԫ�ز���
		y->link = first->link;
		first = y;
	}
	if (!y->link) last = y;
	return *this;
}

template<class T>
void Chain<T>::Output(ostream& out) const
{// ������Ԫ�����������
	ChainNode<T> *current;
	for (current = first; current; current = current->link)
		out << current->data << ' ';
}
//���� <<
template <class T>
ostream& operator<<(ostream& out, const Chain<T>& x)
{
	x.Output(out); return out;
}

template<class T>
void Chain<T>::Erase()
{
	//ɾ�������е����нڵ�
	ChainNode<T> *next;
	while (first) {
		next = first->link;
		delete first;
		first = next;
	}
	first = NULL;
	last = NULL;
}

template < class T >
Chain<T>& Chain<T>::Append(T& x)
{
	//������β�����x
	ChainNode<T> *y;
	y = new ChainNode<T>;
	y->data = &x; y->link = 0;
	if (first) {//����ǿ�
		last->link = y;
		last = y;
	}
	else // ����Ϊ��
		first = last = y;
	return *this;
}

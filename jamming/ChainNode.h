#pragma once
#pragma once
template<class T>
class ChainNode
{
	template<class T> friend class Chain;
	template<class T> friend class ChainIterator;
	template<class T> friend class CircleChain;
	friend class box;
private:
	T* data;
	ChainNode<T> *link;
};



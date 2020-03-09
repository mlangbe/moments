/*
	operations on moment tensors, basis transform code gen + impl.

	This file is part of the source code used for the calculation of the moment invariants
	as described in the dissertation of Max Langbein
	https://nbn-resolving.org/urn:nbn:de:hbz:386-kluedo-38558

	Copyright (C) 2020 TU Kaiserslautern, Prof.Dr. Hans Hagen (AG Computergraphik und HCI) (hagen@cs.uni-kl.de)

	Author: Max Langbein (m_langbe@cs.uni-kl.de)

	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef POINTCLOUD_H
#define POINTCLOUD_H
#ifndef GNUCC
#pragma once
#endif

#include<string>
#include<ostream>
#include "tensor.h"

typedef long double mytype;


/** an iterator over the components of a
* tensor as computed in add3doords
*/
class tsetiter
{

	


	int 
		N, //max.order
		i, // N - xorder
		j, // i - yorder
		k; // j - zorder
public:
	tsetiter(int NN)
	{
		i=j=k=N=NN;
	}

	void operator ++ ()
	{
		if(k >0 )--k;
		else {
			if(j>0)
			{ 
				--j;
				k=j;
			}
			else{
				--i;
				k=j=i;
			}		
		}
	}


	//get the number of indices that are 1, 2, and 3 in o
	//and return the sum of them (=the total order)
	int getOrders(int o[3])
	{
		o[0]= N-i;
		o[1]= i-j;
		o[2]= j-k;
		return N-k;
	}


};





/* return the sum of the number of independent components that 
* the symmetric tensors of order 0 .. N have
*/
inline int getAsiz(int N)
{
	return (((N+6)*N+11)*N)/6 +1;
}


/* calculate the tensor component index
* in a tensorset
* from the number of 1,2s,3s
* tsetiter should be usable for the inversion
* (only for 3d total symmetric tensors)
*\param N: maximum order in the tensorset
*\param I: number of indices that are 0
*\param J: number of indices that are 1
*\param K: number of indices that are 2
*/
inline int calcIndex(int N,int I, int J, int K)
{
	return ( ( I -3*N-6 )*I  +3*(N+2)*(N+2)-1 )*I / 6
					+ ( 2*(N-I) -J+3 )* J / 2
					+ K;
}


/** calculate the index of the component of a 3d symmetric tensor
*with order = I+J+K
*\param I: number of indices that are 0
*\param J: number of indices that are 1
*\param K: number of indices that are 2
*/
inline int calcIndex2(int I,int J,int K)
{
	return I*(I+J+K+1)- (I-1)*I/2+ J;
}

/** calculate the number of independent components
*of a total symmetric tensor of order n and dimension 3
*/
inline int getTsiz(int n)
{
	return (n+1)*(n+2)/2;
}


/**
* evaluate the folding of  the Tensor of order ord
* in the tensorset with max order N
* with the vector x
* e.g order 3: sum_ijk x_i x_j x_k A(3)_ijk
*/
template<class T,int ord,int N>
T evalTotalFold(const T*A,const T*x);

/* tcreate a partial specialization of evalTotalFold*/
void createTotalFold(std::ostream &out, int N, int ord, const std::string& T);

mytype* add3dcoords(int N,mytype *A,const mytype *x,mytype weight=1);

template<int N>
struct translate
{
	translate(mytype*newA, const mytype*A, const mytype *t);
};


/** optimzed version for N>=1 */
template<int N>
struct add3dcoords2
{
	add3dcoords2(mytype *A,const mytype *x);
};

/**
calculate (n over k)  =( n!/( k! (n-k)! ))
*/
inline double nuk(int n,int k)
{
  double ret = 1;
  int i;
  if(k>n/2)k=n-k;
  for(i=0;i<k;i++){ret*=n-i;ret/=i+1;}
  return ret;
}


/*
 calculate the number of multiindices 
 that have i times 0, j times 1 and k times 2
*/
inline double calcNumMultiIndices(int i,int j,int k)
{
	return nuk(i+j+k,i)*nuk(j+k,j);
}


/**
*a tensor that represents one of the tensors in a tensorSet3d
*/
template<class T>
class tensFromTset:public tens<T>
{
	#ifndef OLD_VC
	using odtens::o;
	using odtens::ord;
	#endif

	const T * tset;
	int tsetord;

public:

	T operator[] (const int*index) const
	{
		int numinds[3];
		for(int i=0;i<3;++i)
			numinds[i]=0;

		for(int i=0;i<o;++i)
			++numinds[index[i]];
		
		return tset[calcIndex(tsetord,numinds[0],numinds[1],numinds[2])];	
	}

	tensFromTset(int ord,int tsord,const T* ts)
		:tens<T>(ord,3),tsetord(tsord),tset(ts)
	{		
	}

};


template<class T>
class tensFromTsetVector:public tens<T>
{
	#ifndef OLD_VC
	using odtens::ord;
	using odtens::dim;
	#endif

	const T * tsetv;
	int tsetord;

public:

	T operator[] (const int*index) const
	{

		tensFromTset<T> tt(
			ord()-1, 
			tsetord, 
			//the indices in vector are inverted
			tsetv + (dim() - index[ord()-1] -1 )  *  getAsiz(tsetord) );
		
		return tt[index];	
	}

	tensFromTsetVector(int ord,int tsord,const T* tsv)
		:tens<T>(ord,3),tsetord(tsord),tsetv(tsv)
	{		
	}

};






/*
\param A: the basis-transformed tensor set
\param B: the original tensor set
\param M: the basis-transformation matrix
\post: A is the basis-transformed version of B
*/
template<class T>
void basisTransform(int N,T*A,const T*B,T M[3][3])
{
//x,y,z are the number of indices in A 
//that are x, y or z respectively.

	int *indbuf=new int[N];

	for(int x=0;x<=N;++x)
	{
		for(int y=0;y<=N-x;++y)
		{
			for(int z=0;z<=N-x-y;++z,++A)
			{
				int ord=x+y+z;

				tensFromTset<T> t(ord,N,B);
				basisTransformed<T> bt(t,&M[0][0]);

				//re-create indices
				int *ib=indbuf;
				for(int i=0;i<x;++i,++ib)
					*ib=0;
				for(int i=0;i<y;++i,++ib)
					*ib=1;
				for(int i=0;i<z;++i,++ib)
					*ib=2;

				*A = bt[indbuf];

				//not yet implemented!
				assert(false);
			}		
		}			
	}

	delete[] indbuf;

}

/*
* a version for a vector of symmetric tensor sets
* -> the last index isn't symmetric
\param N: the max order in the symmetric tensors
\param A: the basis-transformed tensor set
\param B: the original tensor set
\param M: the basis-transformation matrix
\post: A is the basis-transformed version of B
*/
template<class T>
void basisTransformV(int N,T*A,const T*B,T M[3][3])
{
//x,y,z are the number of indices in A 
//that are x, y or z respectively.

	int *indbuf=new int[N+1];

	for(int vidx=0;vidx<3;++vidx)
	for(int x=0;x<=N;++x)
	{
		for(int y=0;y<=N-x;++y)
		{
			for(int z=0;z<=N-x-y;++z,++A)
			{
				int ord=x+y+z+1;

				tensFromTsetVector<T> t(ord,N,B);
				basisTransformed<T> bt(t,&M[0][0]);

				//re-create indices
				int *ib=indbuf;
				for(int i=0;i<x;++i,++ib)
					*ib=0;
				for(int i=0;i<y;++i,++ib)
					*ib=1;
				for(int i=0;i<z;++i,++ib)
					*ib=2;
				
				//inverted indices
				*ib=2-vidx;

				*A = bt[indbuf];
			}		
		}			
	}

	delete[] indbuf;

}



void createBasisTransform(std::ostream&o,int N);

void testadd3dc();

#endif

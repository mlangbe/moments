/*
	representations of 3d symmetric tensors sedts (order 0..N) and their calulations from tensor products,
	plus translations.

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
/** a set of 3D symmetric tensors with max.order N */

#ifndef GNUCC
#pragma once
#endif
#ifndef TENSORSET_H
#define TENSORSET_H
#include<memory.h>
#include "pointCloud.h"

template<int N=4,class T=long double>
struct tensorSet3D{
	
	typedef T tval;
	static const int maxord = N;



	static const int ncomps=(((N+6)*N+11)*N)/6 +1;
	
	tval A[ ncomps ];

	/* calculate the tensor component index
	* from the number of 1,2s,3s
	*/
	template<int I,int J,int K>
	struct calcIndex
	{

		static const int 
			idx= ( ( I -3*N-6 )*I  +3*(N+2)*(N+2)-1 )*I / 6
						+ ( 2*(N-I) -J+3 )* J / 2
						+ K;
	};

	//calcindex(N,1,0,0)
	


	void translated(const tensorSet3D<N,T> &ts,
									const tval*t)
	{
	  typename  ::translate<N>::translate(A,ts.A,t);
	}

	void translate(const tval*t)
	{
		tensorSet3D buf=*this;
		translated(buf,t);
	}
	
	static const int _1A0= calcIndex<1,0,0>::idx;
	static const int _1A1= calcIndex<0,1,0>::idx;
	static const int _1A2= calcIndex<0,0,1>::idx;

	void getCenter(mytype cent[3])
	{
		cent[0]=A[_1A0]/A[0];
		cent[1]=A[_1A1]/A[0];
		cent[2]=A[_1A2]/A[0];	
	}

	void transinv()
	{
		tval t[3];
		getCenter(t);
		t[0]=-t[0];t[1]=-t[1];t[2]=-t[2];
		translate(t);								
	}


	void add3dcoords(const tval*x)
	{
		::add3dcoords2<N>(A,x);
	}


	tensorSet3D& reset()
	{
		memset(A,0,sizeof(A));
		return *this;
	}

	tensorSet3D& operator += (const tensorSet3D&t)
	{
		for(int i=0;i<ncomps;++i)
			A[i] += t.A[i];
		return *this;
	}

};


#endif

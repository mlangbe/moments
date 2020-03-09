/*
	implementations of invariant computing methods from moment tensors.

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

#ifndef GNUCC
#pragma once
#endif

#ifndef MOMENTSET_H
#define MOMENTSET_H

/**
base class for momentset,
which is common header for all momentsets

for every id, 
the metainformation has to be given.

\param T: number type(could be polynomials)
\param nmoments: number of moments
\param momentsetid: id of the moment set for this number of moments
*/
template<class T,int nmoments,int momentsetid=0>
class momentSetBase{

public:


	typedef  T valuetype;
	static const int n=nmoments;
	
	valuetype M[n];

	struct metaentry
	{
		//value order
		int vord;
	
		//domain order
		int dord;

		//optional: description of this moment
		const char * des;
	};

	struct meta{
		
		//number of entries in the input array
		int numinputs;
		
		//dimension,order of the source tensors
		int dim,ord;

		static const int n=momentSetBase::n;		
		//one entry for each moment
		metaentry me[n];	

		//optional: descriptions
		const char*des;
	};

	/*pointer to meta-information: hast to be provided
		*by every specialization!*/
	static const meta metainfo;


	
	void scaleValue(const valuetype& sc)
	{
		for(int i=0;i<n;++i)
			{
				const metaentry & e =metainfo.me[i];
				M[i]*=pow(sc,e.vord);
			}	
	}

	void scaleDomain(const valuetype& sc)
	{
		for(int i=0;i<n;++i)
			{
				const metaentry & e =metainfo.me[i];
				M[i]*=pow(sc,e.dord);
			}	
	}


	valuetype maxNormDomain()
	{
		valuetype start=0;
		for(int i=0;i<n;++i)
			{
				const metaentry & e =metainfo.me[i];
				valuetype t = pow(fabs(M[i]),1.0/e.dord);
				if(t>start)start=t;
			}				
		return start;
	}

	valuetype maxNormValue()
	{
		valuetype start=0;
		for(int i=0;i<n;++i)
			{
				const metaentry & e =metainfo.me[i];
				valuetype t = pow(fabs(M[i]),1.0/e.vord);
				if(t>start)start=t;
			}				
		return start;
	}


	valuetype valnormOf(int i)
	{
		const metaentry & e =metainfo.me[i];
		return pow(fabs(M[i]),1.0/e.vord);
	}

	valuetype domnormOf(int i)
	{
		const metaentry & e =metainfo.me[i];
		return pow(fabs(M[i]),1.0/e.dord);
	}

	

};

/*
\param implementationid: 
if multiple implementations
for the compute function exist,
the id of the implementation
*/
template
<
	class T,
	int nmoments,
	int momentsetid=0,
	int implementationid=0
>
class momentSet : public momentSetBase<T,nmoments,momentsetid> 
{
public:
	typedef momentSetBase<T,nmoments,momentsetid> base;
#ifndef OLD_VC
	using typename momentSetBase<T,nmoments,momentsetid>::valuetype;
	using momentSetBase<T,nmoments,momentsetid>::M;
#endif
	/** this function has to be implemented for every	 specialization!*/
	static void compute(valuetype*M,const valuetype*A);

	//a comment about the type:
	static const char * const comment;

	void compute(const valuetype*A)
	{	compute(M,A);	}

};


typedef momentSet<long double,28> faceMomSet;

#endif

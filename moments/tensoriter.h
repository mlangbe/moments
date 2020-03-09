/*
	iterators over the independent components of tensors with various symmetry types.

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
#ifndef TENSORITER_H
#define TENSORITER_H



/* the type of symmetry
		ASSYMM: no symmetry assumed
		SYMM: symmetric
		ANTISYMM: totally antisymmetric

		for the interchanging of the multi indices

		if children.size()!=0, the type is
		used for the interchange of index groups,
		(one group consists of all indices of all subtrees)
	*/
enum _symmtype_ { ASSYMM, SYMM, ANTISYMM , TREE };



/** calculate the max.value for the lower-level
* iterators 
*/

template<int d,unsigned bpi,int o>
struct indsmax { 
	static const unsigned imax 
		= ( indsmax<d,bpi,o-1>::imax << bpi )+d-1;
};

template<int d,unsigned bpi>
struct indsmax<d,bpi,1>
{
	static const unsigned imax = d-1;
};

/*max index for antisymm(strictly decreasing)*/

template<int d,unsigned bpi,int o>
struct indsmaxanti { 
	static const unsigned imax 
		= ( indsmaxanti<d,bpi,o-1>::imax << bpi )+d-o;

	static const unsigned imin 
		= indsmaxanti<d,bpi,o-1>::imin | ( (o-1)<<(bpi*(o-1)) );
};

template<int d,unsigned bpi>
struct indsmaxanti<d,bpi,1>
{
	static const unsigned imax = d-1;
	static const unsigned imin = 0;
};



/** an iterator that
* iterates through all independant components
* in a symmetric tensor
* used primarily for building index map;
* in difference to symminditer the 
* indices are sorted downwards now
* (x[0]>=x[1]>=x[2]),
* still the last is the fastest moving index
 *(if possible)
 *\param dim dimension
 *\param ord order
 *\param nbits number of bits for every index
*/
template<int dimen,int order,unsigned bitsperindex=2>
struct symmii
{
  public:

	static const unsigned 
		ord=order,dim=dimen,nbits = bitsperindex,
		ordm1=ord-1, dimm1=dim-1, indbits=(1<<nbits)-1;

	//--------------------

	///the indices of the tensor, sorted
	///(last one in least-significant bits)
	unsigned inds;

	unsigned idx;

	//--------------------

	static int symmetrytype() {return SYMM;}

	symmii()
	{
		//dim-1 should fit in indbits bits
		assert( (dimm1 & indbits)==dimm1);

		//there should be enough bits for the indices to fit
		assert( sizeof( unsigned ) *8 >= nbits*ord ) ;

		inds=0;idx=0;
	}	
		
	//is the iterator at its maximum ?
	bool ismax() const
	{ return ( inds&indbits )== dimm1  ;  }

	void reset(){inds=0;idx=0;}

	bool operator != (const symmii &si) const
	{ return inds!= si.inds  ;  }



	unsigned packinds(const unsigned*inds) const
	{
		unsigned ret=0;
		const unsigned*indsend=inds+ordm1;
		for(;inds!=indsend;++inds)
			ret<<=nbits,ret|=*inds;
		return ret;
	}

	void unpackinds(unsigned indexset,unsigned*inds) const
	{
		unsigned *i= inds+ordm1+1;
		do{
			--i;
			*i = indexset&indbits;
			indexset>>=nbits;
		}while(i!=inds);
	}

	void getindices(unsigned*indices) const{
		unpackinds(inds,indices);	
	}

	const symmii& operator -- ()
	{ 

		assert(inds!=0);
		--idx;

		unsigned nb = 0;

		//look how many indices are 0
		while( ((inds>>nb)& indbits ) == 0 )
			nb += nbits;

		//for sake of efficiency:
		if(nb==0)
		{ --inds; return *this; }

		//decrease the first one not 0
		inds -= 1<<nb;

		//now set all lower indices to 0
		inds &= ~ (( 1<<nb) -1 );

		//set all lower indices to the one just decreased
		nb = inds & (indbits << nb);
		do{
			nb>>=nbits;
			inds|=nb;
		}while(nb);

		return *this;
	}

	const symmii& operator ++()
	{ 
		assert(!ismax());

		++idx;

		//look how many indices are at max.value
		//= are at the value of the next upper index

		unsigned nb, buf;
		if(inds==0)
			nb= ordm1*nbits;
		else{
			nb=0, buf=(inds>>nbits)^inds ;
			while( (buf & indbits ) ==0)
				nb+=nbits, buf>>=nbits;
		}

		//for sake of efficiency:
		if(nb==0)
		{ ++inds; return *this; }

		//decrease the first one not at min value
		inds += 1<<nb;

		//now set all lower indices to 0
		inds &= ~((1<<nb)-1);

		return *this;
	}

};


/** an iterator that
* iterates through all independant components
* in a symmetric tensor
* used primarily for building index map;
* in difference to symminditer the 
* indices are sorted downwards now
* (x[0]>=x[1]>=x[2]),
* still the last is the fastest moving index
 *(if possible)
 *\param dim dimension
 *\param ord order
 *\param nbits number of bits for every index
*/
template<int dimen,int order,unsigned bitsperindex=2>
struct assymmii :public symmii<dimen,order,bitsperindex>
{

	#ifndef OLD_VC
	using base=symmii<dimen,order,bitsperindex>;
	using base::inds;
	using base::indbits;
	using base::nbits;
	using base::idx;
	using base::dimm1;
	using base::ordm1;
	#endif

	static int symmetrytype() {return ASSYMM;}


	static const unsigned
		imax = indsmax<dimen,bitsperindex,order>::imax;
		
	//is the iterator at its maximum ?
	bool ismax() const
	{ return inds == imax;  }



	const assymmii& operator -- ()
	{ 

		assert(inds!=0);
		--idx;

		unsigned nb = 0;

		//look how many indices are 0
		while( ((inds>>nb)& indbits ) == 0 )
			nb += nbits;

		//for sake of efficiency:
		if(nb==0)
		{ --inds; return *this; }

		//decrease the first one not 0
		inds -= 1<<nb;

		//set all lower indices to dim-1
		inds &= ~ (( 1<<nb) -1 );
		inds |= (( 1<<nb) -1 ) & imax ;

		return *this;
	}

	const assymmii& operator ++()
	{ 
		assert(!ismax());

		++idx;

		//look how many indices are at max.value

		unsigned nb=0;
		while(( (inds>>nb) & indbits ) ==dimm1)
			nb+=nbits;
				
		//increase the first one not at min value
		inds += 1<<nb;

		//now set all lower indices to 0
		inds &= ~((1<<nb)-1);

		return *this;
	}

};

/** an iterator that
* iterates through all independant components
* in a antisymmetric tensor
* used primarily for building index map;
* in difference to symminditer the 
* indices are sorted downwards now
* (x[0]>=x[1]>=x[2]),
* still the last is the fastest moving index
 *(if possible)
 *\param dim dimension
 *\param ord order
 *\param nbits number of bits for every index

 **ATTENTION: not yet tested!!,not yet working**
*/
template<int dimen,int order,unsigned bitsperindex=2>
struct antisymmii :public symmii<dimen,order,bitsperindex>
{

	#ifndef OLD_VC
	using base=symmii<dimen,order,bitsperindex>;
	using base::inds;
	using base::indbits;
	using base::nbits;
	using base::idx;
	using base::dimm1;
	using base::ordm1;
	using base::dim;
	using base::ord;
	#endif

	static int symmetrytype() {return ANTISYMM;}


	static const unsigned
		imax = indsmaxanti<dimen,bitsperindex,order>::imax,
		imin = indsmaxanti<dimen,bitsperindex,order>::imin;
		
	//is the iterator at its maximum ?
	bool ismax() const
	{ return inds == imax;  }

	void reset()
	{ inds=imin;idx=0; }

	antisymmii()
	{ 
		assert(dim>=ord);
		reset(); 	
	}



	const antisymmii& operator -- ()
	{ 
		assert(inds!=imin);
		--idx;

		if((inds&indbits)>0)
		{--inds;return *this;}

		//first index is 0:
		//look how many indices are at their min.value
		//= are at the value of the next lower index+1
		unsigned nb = nbits;
		unsigned buf = ( (inds>>nbits) - inds ) & ((1<<nbits*ordm1)-1);
		while( (buf & indbits ) == 1 ) 
			nb+=nbits, buf>>=nbits;

		//decrease the first one not at min value
		inds -= 1<<nb;

		//now set all right indices to the maximum decreasing order possible
		buf = (inds>>nb)&indbits;
		//set the indices right of the decreased one to zero
		inds &= ~ (( 1<<nb) -1 );
		//fill with decreasing buf
		do{
			nb-=nbits;--buf;
			inds |=  buf<<nb ;
		}while(nb!=0);

		return *this;
	}

	const antisymmii& operator ++()
	{ 
		assert(!ismax());
		++idx;	
		

		//look how many indices are at their max.value
		//= are at the value of the next upper index-1
		unsigned nb = 0,buf = ( (inds>>nbits) - inds ) & ((1<<nbits*ordm1)-1);
		while( (buf & indbits ) == 1 ) 
			nb+=nbits, buf>>=nbits;

		//increase the first one not at min value
		inds += 1<<nb;

		//now set all right indices to the minimum decreasing order
		inds &= ~ (( 1<<nb) -1 );
		inds |=  ((1<<nb)-1) & imin ;

		return *this;
	}


};


template<class base,int n>
struct symmiicomposed
{
	base its[n];
	int idx;

	static const int 
		dim = base::dim, 
		ord = base::ord*n
		;
	
	bool ismax() const
	{
		return its[n-1].ismax();	
	}

	void reset()
	{
		for(int i=0;i<n;++i)
			its[i].reset();
		idx=0;
	}

	bool operator != (const symmiicomposed &x) const
	{
		for(int i=0;i<n;++i)
			if(its[i]!=x.its[i])
				return true;

		return false;
	}


	void getindices(unsigned*indices) const{		
		for(int i=0;i<n;++i,indices+=base::ord)
			its[i].getindices(indices);	
	}

	symmiicomposed & operator++()
	{
		assert( !ismax() );

		++idx;

		int i;
		for(i=n-1;i>0;--i)
			if( its[i] != its[i-1] )
				break;

		if(i==n-1)
		{	++its[i]; return *this; }

		++its[i];
		for(--i;i>=0;--i)
			its[i].reset();			
		return *this;
	}

};

template<class base,int n>
struct antisymmiicomposed:public symmiicomposed<base,n>
{
#ifndef OLD_VC
	using symmiicomposed<base,n>::its;
	using symmiicomposed<base,n>::idx;
	using symmiicomposed<base,n>::ismax;
#endif

	void reset()
	{
		its[n-1].reset();
		for(int i=n-2;i>=0;--i)
			its[i]=its[i+1], ++its[i];
		idx=0;
	}

	antisymmiicomposed()
	{ reset(); }

	bool ismax()
	{ 
		if(!its[0].ismax())
			return false;

		base x= its[n-1];
		for(int i=n-2;i>=0;--i)
		{
			++x;
			if(x!=its[i])
				return false;
		}

		return true;

	}

	antisymmiicomposed & operator++()
	{
		assert( !ismax() );

		++idx;

		int i;
		for(i=n-1;i>0;--i)
		{
			++its[i];
			if(its[i]!=its[i-1])
				break;
		}
		if(i==0)
			++its[0];

		if(i==n-1)
		{ return *this; }

		//set all lower indices to the minimum monotonic decreasing order
		its[n-1].reset();
		for(int k=n-2;k>i;--k)
		{
			its[k]=its[k+1];
			++its[k];
		}			
		return *this;
	}

};

template<class base1,class base2>
struct assymmiicomposed
{
	static const int
		ord=base1::ord+base2::ord,
		dim=base1::dim;

	base1 a;
	base2 b;

	unsigned idx;

	assymmiicomposed()
	{
		assert(base1::dim==base2::dim);
	}

	bool ismax() const
	{
		return a.ismax() && b.ismax();	
	}

	void reset()
	{
		a.reset();b.reset();
		idx =0;
	}

	bool operator != (const assymmiicomposed &x) const
	{
		return a!=x.a||b!=x.b;
	}

	void getindices(unsigned*indices) const{		
		a.getindices(indices);
		b.getindices(indices + a.ord);
	}

	assymmiicomposed & operator++()
	{
		assert(!ismax());

		++idx;

		if(a.ismax())
			a.reset(),++b;
		else
			++a;

		return *this;
	}

};

#endif

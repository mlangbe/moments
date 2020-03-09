/*
	some simple bitwise operations needed.

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


#pragma once
#ifndef BITOPS_H
#define BITOPS_H


int myfrexp(float&x)
{
	int & u=*((int*)&x);

	int ret=((u>>23)&0xff)-0x7f;
	
	u = u & ~(0xff<<23) | (0x7f<<23);

	return ret;
}



inline int bigEndian()
{
	double x0=-1;
	return (  *((char*)&x0)<0);
}


int myfrexp(double&x)
{

	int & u=*(((int *)&x)+  !bigEndian() );

	int ret=((u>>20)&0x7ff) -0x3ff;
	
	u = u & ~(0x7ff<<20) | (0x3ff<<20);

	return ret;
}


template<class T,int bits=sizeof(T)*8,int block=sizeof(T)*4> 
struct blockmask
{
private:
	static const T f = blockmask<T,bits,block*2>::e; 
	static const T ff = ( f >> block) &  f ;	
public:	
	// e has now the following structure block bits 1, then block times  0 then again block times 1 and so on.
	static const T e =   (ff << block*2) | ff; 	
};

template<class T,int bits>
struct blockmask<T,bits,bits>
{
	static const T e=~(T)0;
};







template<class T,int bits=sizeof(T)*8,int block=sizeof(T)*4> 
struct bitcounter
{
	
	// e has now structure block bits 1, then block bits 0 then again block bits 1.
	static const T e =   blockmask<T,bits,block>::e; 	
	static T countbits1(T i)
	{
		i=bitcounter<T,bits,block/2>::countbits1(i);
		i = (i&e) + ((i>>block)&e);
		return i;
	}

};

template<class T,int bits>
struct bitcounter<T,bits,0>
{
	static T countbits1(T i)
	{
		return i;
	}
};


template<class T,int bits>
struct bitcounter<T,bits,bits>
{
	static const T e=~(T)0;
};

template<class T>
T countbits1(T i)
{
	return bitcounter<T,sizeof(T)*8,sizeof(T)*4>::countbits1(i);	
}


template<class T,int bits=sizeof(T)*8,int block=sizeof(T)*4>
struct bitfiller
{
	static T fillbitsright(T i)
	{
		i=bitfiller<T,bits,block/2>::fillbitsright(i);
		i|=i>>block;
		return i;
	}
};

template<class T,int bits>
struct bitfiller<T,bits,0>
{
	static T fillbitsright(T i)
	{
		return i;
	}
};


template<class T>
T log2i(T i)
{
	return countbits1(bitfiller<T>::fillbitsright(i));
}


	/**inverse gray code */
	template<class uint> uint invgray(uint y)
	{
		uint x=y;
		x = x^(x>>1);
		x = x^(x>>2);
		x = x^(x>>4);
		x = x^(x>>8);
		x = x^(x>>16);		
		if( (uint(1)>>32)==0 )
			x = x^(x>>32);
		if( (uint(1)>>64)==0 )
			x = x^(x>>64);
		return x;
	}



#endif

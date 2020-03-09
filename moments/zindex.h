/*
	z-index computation methods used e.g. in smooth fractal dimension calculation.

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

#ifndef ZINDEX_H
#define ZINDEX_H
#pragma once

template<int i>
struct pow3
{
	static const int c=pow3<i-1>::c*3;
};

template<>
struct pow3<0>
{
	static const int c=1;
};

template<int i>
struct log3
{
	static const int e=log3<i/3>::e+1;
};

template<>
struct log3<0>
{
	static const int e=0;	
};



template<class T,int blocksiz,int nblocks > 
struct blockmask3
{
	// e has now the following structure block bits 1, then block times  0 then again block times 1 and so on.
	static const T e =  (blockmask3<T,blocksiz,nblocks-1>::e << 3*blocksiz) | (( ((T)1) <<blocksiz )-1);	
};

template<class T,int blocksiz>
struct blockmask3<T,blocksiz,0>
{
	static const T e=0;
};



template<class T,int bits,int level>
struct step{

	static void st( T&x,T&y,T&z)
	{
		static const int n = pow3<level>::c; 

		if(n<bits){	
			T xx=x,yy=y,zz=z;
			static const T a0 =  blockmask3<T,n, (bits +n*3-1)/(n*3)>::e;
			x = (xx&a0) | ((yy&a0)<<n);
			static const T a1=a0<<n;
			y = ((xx&a1)>>n) | (yy&a1) | ((zz&a1)<<n);
			if(2*n<bits){
				static const T a2=a0<<2*n;	
				x|= ((zz&a0)<< 2*n);
				z = ((xx&a2)>>2*n) | ((yy&a2)>>n) | (zz&a2);
				if(3*n<bits){
					step<T,bits,level+1>::st(x,y,z);
				}
			}
		}

	}

};


template<class T,int bits>
struct step<T,bits,8>{
	static void st( T&x,T&y,T&z)
	{
	}
};

template<class T>
T zindex(T x,T y,T z)
{
	step<T,sizeof(T)*8,0>::st(x,y,z);
	return x;
}



template<class T,int bits,int level>
struct invstep{

	static void st( T&x,T&y,T&z)
	{

		static const int n = pow3<level>::c; 

		if(n<bits)
		{	
			T xx=x,yy=y,zz=z;
			static const T c0 =  blockmask3<T,n, (bits +n*3-1)/(n*3)>::e;
			x = (xx&c0)|((yy&c0)<<n);
			static const T c1=c0<<n;
			y = ((xx&c1)>>n) | (yy&c1) | ((zz&c1)<<n);
			if(2*n<bits)
			{
				x|= ((zz&c0)<<2*n);
				static const T c2=c0<<2*n;
				z = ((xx&c2)>>2*n) | ((yy&c2)>>n) | (zz&c2);
			}
			else
			{
				z=0;
			}	 
		}

		invstep<T,bits,level-1>::st(x,y,z);
	}

};


template<class T,int bits>
struct invstep<T,bits,-1>{
	static void st( T&x,T&y,T&z)
	{
	}
};

template<class T>
void invzindex(T& x,T& y,T& z ,T in)
{
	x=in;y=0;z=0;
	invstep<T,sizeof(T)*8,log3<sizeof(T)*8>::e+1>::st(x,y,z);
}



    template<class uint>
	static inline uint zord(uint x, uint y, uint z)
	{
		uint xx,yy,zz;
		
		uint xyz = x|y|z;
		//binary ...01001001001001001001001001001001
		static const uint A0=0111111111u; 
		static const uint a0=(((A0<<27)|A0)<<27)|A0;
		xx = (x&a0) | ((y&a0)<<1) | ((z&a0)<<2);
		if( (xyz>>1)==0 ) 
			return xx;
		static const uint a1=a0<<1;
		yy = ((x&a1)>>1) | (y&a1) | ((z&a1)<<1);
		if( (xyz>>2)==0 ) 
			return xx|(yy<<3);
		static const uint a2=a0<<2;
		zz = ((x&a2)>>2) | ((y&a2)>>1) | (z&a2);
		
		//binary ...00111000000111000000111000000111
		static const uint B0=07007007;
		static const uint b0=(((B0<<27)|B0)<<27)|B0; 
		x =  (xx&b0) | ((yy&b0)<<3) | ((zz&b0)<<6);
		if( (xyz>>3)==0 ) 
			return x;
		static const uint b3=b0<<3;
		y =  (xx&b3)>>3 | ((yy&b3)) | ((zz&b3)<<3);
		if( (xyz>>6)==0 ) 
			return x|(y<<9);
		static const uint b6=b0<<6;
		z =  (xx&b6)>>6 | ((yy&b6)>>3) | ((zz&b6));		

		//binary ...11000000000000000000000111111111
		static const uint C0=0777;
		static const uint c0=(((C0<<27)|C0)<<27)|C0; 
		xx =  (x&c0) | ((y&c0)<<9) | ((z&c0)<<18);
		if( (xyz>>9)==0 ) 
			return xx;
		static const uint c9=c0<<9;
		yy =  (x&c9)>>9 | ((y&c9)) | ((z&c9)<<9);
		if( (xyz>>18)==0 ) 
			return xx | (yy<<27);
		static const uint c18=c0<<18;
		zz =  ((x&c18)>>18) | ((y&c18)>>9) | ((z&c18));		

		//binary ...00111111111111111111111111111111
		static const uint d0=0777777777; 
		if((xyz>>27)==0 )		
			return xx | (yy<<27) | (zz<<54);	
		else
			return ~0Lu;
	}



	template<class uint>
	static inline uint zord2(uint x, uint y, uint z)
	{
		uint xx,yy,zz;
		
		//binary ...01001001001001001001001001001001
		static const uint A0=0111111111u; 
		static const uint a0=(((A0<<27)|A0)<<27)|A0;
		xx = (x&a0) | ((y&a0)<<1) | ((z&a0)<<2);
		static const uint a1=a0<<1;
		yy = ((x&a1)>>1) | (y&a1) | ((z&a1)<<1);
		static const uint a2=a0<<2;
		zz = ((x&a2)>>2) | ((y&a2)>>1) | (z&a2);
		
		//binary ...00111000000111000000111000000111
		static const uint B0=07007007;
		static const uint b0=(((B0<<27)|B0)<<27)|B0; 
		x =  (xx&b0) | ((yy&b0)<<3) | ((zz&b0)<<6);
		static const uint b3=b0<<3;
		y =  (xx&b3)>>3 | ((yy&b3)) | ((zz&b3)<<3);
		static const uint b6=b0<<6;
		z =  (xx&b6)>>6 | ((yy&b6)>>3) | ((zz&b6));		

		//binary ...11000000000000000000000111111111
		static const uint C0=0777;
		static const uint c0=(((C0<<27)|C0)<<27)|C0; 
		xx =  (x&c0) | ((y&c0)<<9) | ((z&c0)<<18);
		static const uint c9=c0<<9;
		yy =  (x&c9)>>9 | ((y&c9)) | ((z&c9)<<9);
		static const uint c18=c0<<18;
		zz =  ((x&c18)>>18) | ((y&c18)>>9) | ((z&c18));		
		return xx | (yy<<27) | (zz<<54);	
	}


	

	/** optimized inv zord*/
	template<class uint> void invzord(uint&x,uint&y,uint&z,uint idx)
	{
		uint xx,yy,zz;
		
		static const uint D0=0777777777;
		xx=idx&D0;
		yy=(idx>>27)&D0;
		zz=(idx>>54)&D0;
	
		//binary ...11000000000000000000000111111111
		static const uint C0=0777;		
		static const uint c0=(((C0<<27)|C0)<<27)|C0,c9=c0<<9,c18=c0<<18; 

		x = (xx&c0) | ((yy&c0)<<9) | ((zz&c0)<<18);
		y = ((xx&c9)>>9) | (yy&c9) | ((zz&c9)<<9);
		z = ((xx&c18)>>18) | ((yy&c18)>>9) | (zz&c18); 

		//binary ...00111000000111000000111000000111
		static const uint B0=07007007;
		static const uint b0=(((B0<<27)|B0)<<27)|B0,b3=b0<<3,b6=b0<<6; 

		xx = (x&b0) | ((y&b0)<<3) | ((z&b0)<<6);
		yy = ((x&b3)>>3) | ((y&b3)) | ((z&b3)<<3);
		zz = ((x&b6)>>6) | ((y&b6)>>3) | ((z&b6));

		//binary ...01001001001001001001001001001001
		static const uint A0=0111111111u; 
		static const uint a0=(((A0<<27)|A0)<<27)|A0,a1=a0<<1,a2=a0<<2;
	
		x=(xx&a0) | ((yy&a0)<<1) | ((zz&a0)<<2);
		y=((xx&a1)>>1)| (yy&a1) | ((zz&a1)<<1);
		z=((xx&a2)>>2)| ((yy&a2)>>1) | ((zz&a2));

	}



#endif

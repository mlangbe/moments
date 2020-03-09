/*
	Some nice approximations to mathematical functions.

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


#ifndef APPROXMATH_H
#define APPROXMATH_H
#pragma once

#include<cmath>
#include<stdio.h>
#include "bitops.h"


template<int n,class num>
inline num pow2n(num x)
{
	for(int i=0;i<n;++i)
		x*=x;
	return x;
}







/**
approximate log_2(i) with precision 2^-8
*/
struct logtab
{

protected:
	static const int pw=8;
	static const int tabsiz=1<<pw;
	float table[tabsiz]; 
public:	
	logtab()
	{
		printf("\ninit logtab\n");
		for(int i=1;i<tabsiz;++i)
			table[i]=log((float)i);
		table[0]=-1e+30;
	}
	




	inline float logi(unsigned i) const
	{
		unsigned cnt=0;

		unsigned dp=pw/2;
		while(i>=tabsiz)
		{
			i>>=dp;
			cnt+=dp;
		}
		return cnt*table[2]+table[i];
	}

	inline float logf(float x) const
	{
		static const int mul=tabsiz>>1;

		int i;
		i=myfrexp(x);
		float ret=i*table[2];
		x*=mul;
		ret -= table[mul];
		int ii=(int)x;
		ret+=table[ii];
		x/=ii;	

		x-=1;


		return ret + (-.5*x+1)*x ;
	}

};


template<typename T, int l2>
struct glogtab
{

protected:
	static const int pw=l2;
	static const int tabsiz=1<<pw;
	T logtab[tabsiz]; 
	T divtab[tabsiz];
public:	
	glogtab()
	{
		printf("\ninit logtab\n");
		for(int i=1;i<tabsiz;++i){
			logtab[i]=(T)::log(i);
			divtab[i]=((T)1)/((T)i);
		}
		logtab[0]=(T)-HUGE_VAL;
		divtab[0]=(T)HUGE_VAL;
	}
	




	inline T logi(unsigned i) const
	{
		unsigned cnt=0;

		unsigned dp=pw/2;
		while(i>=tabsiz)
		{
			i>>=dp;
			cnt+=dp;
		}
		return cnt*logtab[2]+logtab[i];
	}

	inline T log(T x) const
	{
		static const int mul=tabsiz>>1;

		//invariant: x * exp(ret) = const
		int i;
		//x = frexp(x,&i);
		i=myfrexp(x);
		T ret=i*logtab[2];

		x*=mul;
		ret -= logtab[mul];

		int ii=(int)x;

		ret+=logtab[ii];
		x*=divtab[ii];	

		x-=1;		
		return (((-3 *x + 4)*x -6)*x + 12)*x  * divtab[12]  + ret;
	}

};





inline float mylogi(int i)
{
	static const logtab l;

	return l.logi(i);
}

inline float mylogf(float i)
{
	static const logtab l;

	return l.logf(i);
}






template<int n,typename flp>
struct sqtab{
	/** contains 2^( 2^(1-i)  )*/
	flp table[n];	
	//1/entries of table
	flp itable[n];	
	sqtab()
	{
		flp otab=2;
		for(int i=0;i<n;++i)
		{
			flp nt = sqrt(otab);
			table[i]= nt;				
			itable[i]=1/nt;
			otab=nt;
		}
	}	
};


//try to have exact logf. just for fun, no extreme performance expected, as too many 
//non-inline basic functions are called
inline float logf(float x1)
{
	static const float log_2=(float)log(2);
	static const float _1_12 =(float)(1.0/12);

	//rule:  nsq *(order+1) = bits of mantissa - log2(order+1)
	static const int nsq=4;
	static const sqtab<nsq,float> sq;
	
	
	float log2_x =  myfrexp(x1);
	float x=x1;
	float _2_i1 =0.5;//2^-(i+1)	
	//loop condition: x  *  2 ** (log2_x) stays constant
	for(int i=0;i<nsq; ++i,_2_i1 *=.5)
	{
		if(x > sq.table[i]) //x >  2 ** ( + 2 ** (-(i+1)) )
		{			
			x *= sq.itable[i]; //x *=   2 ** ( - 2 ** (-(i+1)) )
			log2_x += _2_i1 ; //log2x += 2** (-i+1)
		}
	}	
	
	x-=1;		
	return   log_2 * log2_x  +  _1_12 * (((-3 *x + 4)*x -6)*x + 12)*x ;
}

//try to have exact logf. just for fun, no extreme performance expected, as too many 
//non-inline basic functions are called
inline double logd(double x)
{
	static const double log_2=(double)logl(2);
	static const double _1_60 =(double)(1.0/(long double)60);

	//bits to be cropped off by bisection with the table values before applying taylor series.
	//rule:  nsq *(order+1) = bits of mantissa - log2(order+1)
	static const int nsq=7;
	static const sqtab<nsq,double> sq;
	
	int ee;
	x=frexp(x,&ee);	
	double log2_x =  ee;//myfrexp(x);
	double _2_i1 =0.5;//2^-(i+1)	
	//loop condition: x  *  2 ** (log2_x) stays constant
	for(int i=0;i<nsq; ++i,_2_i1 *=.5)
	{
		if(x < sq.itable[i]) //x <  2 ** ( - 2 ** (-(i+1)) )
		{			
			x *= sq.table[i]; //x *=   2 ** ( + 2 ** (-(i+1)) )
			log2_x -= _2_i1 ; //log2x -= 2** (-i+1)
		}
	}	
	
	x-=1;		
	return   log_2 * log2_x  +  _1_60 *
	(((((
	 -10.*x //order 6: -1/6
	 +12.)*x
	 -15.)*x
	 +20.)*x
	 -30.)*x
	 +60 )*x;
	;
}







/**
approximate 2^x with precision of 2^-9 assuming 0<=x<31
*/
inline float pow2 (float x) 
{
		static const float c1=log(2.f)/2;
		static const float c2=c1*c1/2;
		int ix=(int)(x+.5);
		
		float r = (x-ix);

		//error is approx  x^3/6

		float exr =( c2*r +c1)* r + 1.f;

		//float rr=ret*ret;

		//*((unsigned*)&rr)+=ix<<23;return rr;

		return exr*exr*(1<<ix);
}




template<int pw,class num>
inline num bellcurve(const num& x)
{
    //approximate exp(x) by (1 + x/n)^n,
	//consequently bell(x)=exp(-.5*x*x) by (1  - x*x/(2n))^n.
	//x^(2^n) = squared n times
	//static const int pw = 8;
	static const int n = 1<<pw;
	static const num factor=0.5/n;
	//by factor=1.0 / 64, this means that it gets zero at x=16.
	num y= 1.0 - x*x*factor;

	

	if(y<0) 
		return 0;


	return pow2n<pw>(y);
}

template<class num>
inline num bellcurve(const num& x)
{
	return bellcurve<7>(x);
}

#endif

/*
    long double implemntations, once used for infnum on 128-bit long doubles on a big-endian architecture.


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


#ifndef MYLMATH_HH
#define MYLMATH_HH

#ifndef ldexpl

#ifndef __sparc
#error "defs in mylmath only for sparc architecture"
#endif



#include<float.h>
#include<sys/int_types.h>
#include<stdio.h>

inline long double ldexpl(long double x,int e)
{
 int16_t &ex = *((int16_t*)&x);

 e+= (ex & 0x7fff)-0x3fff;

 if(e>=LDBL_MAX_EXP | e<LDBL_MIN_EXP | x==0 ){
   if(x==0){ 
     return 0;
   }
   if(e<LDBL_MIN_EXP){
     printf("denormalized: e:%d x:%Lf\n",e,x);
     return 0;
   }
   return LDBL_MAX;
 }


 ex = ex &0x8000 | (e+0x3fff)&0x7fff;

 return x;

}

inline long double frexpl(long double x,int*e)
{
 int16_t &ex = *((int16_t*)&x);
 *e = (ex&0x7fff) - 0x3ffe; 
 ex = (ex&0x8000) | 0x3ffe;
 return x;
}


inline long double floorl(long double x)
{

  long double ret = x;
  char*p = (char*)&ret;

  int zbits =  ( *((int16_t*)p)  & 0x7fff ) - int16_t(0x3ffe)  ;

  if(zbits<=0) 
    return -(x<0);

  zbits = int(LDBL_MANT_DIG)-zbits;

  if(zbits<=0) //if no bits of mantissa have to be set 0
    return x;
 
  int zbytes= zbits>>3;
  zbits &= 7;

  //set mantissa bits 0
  memset(p + (sizeof(ret)- zbytes) , 0 , zbytes );
  p[sizeof(ret)-zbytes-1] &= 0xff << zbits;       

  //-(ret>x) means to decrease it by one if rounded up 
  return ret-(ret>x);
}

#endif

#endif

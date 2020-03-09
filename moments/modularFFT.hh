/*
    modular fft which may be used to optimize the multiplication of very large numbers.

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

#ifndef MODULARFFT_HH
#define MODULARFFT_HH

#include "primtest.hh"
/**
 *class that implememts a restclass(Z/mZ)
 *\pre: myint is a signed integer,
 * m fits into myfloat without rounding,
 * m is positive
 */
template<class myint,class myfloat,myint m>
struct restclass
{
  myint a;

  inline restclass(){}
  inline restclass(const myint&c){ a=c; }

  static const int signbit=sizeof(myint)*8-1;

  inline void operator *=(const myint&b){

    myint x = myint(myfloat(a)*b/m)+1;
    a= (a*b - x*m);
    //    if(a<0)a+=m;
    a+=m&(a>>signbit);
    a+=m&(a>>signbit);
  }

  /** compute the modular inverse 
   *\pre greatest common divisor of (*this,m)==1
   *\post return * *this ==1
   */
  restclass inverse() const;

  inline void operator /=(const myint& b)
  { *this *= b.inverse(); }
  
  inline operator myint() const
  {return a;}

  inline void operator +=(const myint&b)
  {a+=b-m;a+=m&(a>>signbit);/*if(a<0)a+=m;*/}

  inline void operator -=(const myint&b)
  {a-=b;a+=m&(a>>signbit);/*if(a<0)a+=m;*/}


};

#define RESTCLASS_BINARY_FROM_UNARY(op)\
template<class myint,class myfloat,myint m> \
inline \
restclass<myint,myfloat,m> \
operator op(restclass<myint,myfloat,m> a,\
	    const restclass<myint,myfloat,m>& b) \
{ a op##= b; return a; } 

RESTCLASS_BINARY_FROM_UNARY(+)
RESTCLASS_BINARY_FROM_UNARY(-)
RESTCLASS_BINARY_FROM_UNARY(*)
RESTCLASS_BINARY_FROM_UNARY(/)

#undef RESTCLASS_BINARY_FROM_UNARY


template<class myint,class myfloat,myint m>
restclass<myint,myfloat,m>
restclass<myint,myfloat,m>::inverse() const
{



  if( a==1 | a==m-1 ) return a;

  myint x0,x1, r0,r1, rem,quot;

  if(euklid(a,m,r0,r1)==1)
    return r0;
  else
    throw;

  quot = m/a;
  x1 = a;
  r1 = 1;

  //x0 = m % xx, r0= -m/xx;
  x0 = m - quot*a;
  r0 = -quot;

  //if |x1| < |x0|, exchange 1 and 0;
  if( (x0>x1)^(x1<0) )
    {rem=x0;x0=x1;x1=rem;
    rem=r0;r0=r1;r1=rem;}

  while( x0!=1 & x0!=-1 & x0!=0 ){

    quot=x1/x0;
    rem = x1 - x0 * quot;
    x1=x0;x0=rem;
    rem = r1 - r0 * quot;
    r1=r0;r0=rem;
  } 

  //     cout<<x1<<"="<<((long long)r1*xx%m)<<"="<<r1<<"*"<<xx<<endl;
  //     cout<<x0<<"="<<((long long)r0*xx%m)<<"="<<r0<<"*"<<xx<<endl<<endl;
  r0*=x0;

  
  if(r0<0)r0+=m;
  return r0;
}



struct modularpolynFFTBase
{
  typedef long long myint;

  //number of 2's 
  //in faktorization of myprime-1
  static const int  maxpw=48;

  //number of binary digits of myprime
  static const int  pbits=63;

  //prime number to make modulus;
  static const myint
  myprime = (myint(0x3fdc)<<maxpw)+1;
  
  //myprime should fit into myfloat without rounding
  typedef long double myfloat;

  typedef restclass<myint,myfloat,myprime> myrc;

  //element with ord(x[i])=i
  myrc eiwpw[maxpw+1];

  //elements with ord(x[i])=i
  //and x[i]=eiwpw[i]^-1 
  myrc eiwpwinv[maxpw+1];

  //2^-pw (mod myprime)
  myrc minv[maxpw+1];

  modularpolynFFTBase();

private:

  static inline int spiegbits(int x,int n);

  template<class ta,class tb>
  static inline ta pow(ta a,tb b);

  void xchvals(myrc * x,int pw,bool invert) const;

public:
  //different versions of FFT
  void FFT_v1(myrc * x, 
	      int pw,
	      bool invert=false) const;

  void FFT_v2(myrc * x, 
	      int pw,
	      bool invert=false) const;

  void FFT_v3(myrc * x, 
	      int pw,
	      bool invert=false) const;

  void FFT(myrc * x, 
	   int pw,
	   bool invert=false) const;
}; 


extern const 
modularpolynFFTBase modularPolynFFT;

int testrestclass();
int testmodularFFT(int pw);

#endif

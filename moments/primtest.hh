/*
	test some prime-number operations

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

#ifndef PRIMTEST_HH
#define PRIMTEST_HH

#include "infnum.hh"

typedef infnum::t1elem ulint;


struct gencandidate
{

  typedef unsigned long ulint;  
  ulint n,phin;    
  ulint *a;

  gencandidate(ulint level);

  ~gencandidate()
  { delete[]a; }
  
}; 

template<class Tint>
Tint gen(const gencandidate&g, Tint i );

extern gencandidate gc5; 

template<class Tint>
Tint gen(const gencandidate&g, Tint i )
{
  return Tint(i/g.phin)*g.n + Tint(g.a[i%g.phin]);
}

//gen(gc5,i) gives now all numbers 
//not dividable by the first 5 prime numbers



template<int nn=10000>
struct primedb
{
  static const ulint n=nn;
  
  ulint p[n];
    
  inline ulint * end(){ return p+n; }
  inline ulint * begin(){return p; }
  
  
  primedb();


  //tests if n primzahl
  //1 = ja, 0 = nein, -1 = weiss nicht
  template<class Tint>
  int test( const Tint& n );

};

extern primedb<> primzahlen;
extern gencandidate gc;

class mythreadvars
{
  infnum * next; //next number to test, if 0,there's none
  int nextfakt; //next prime faktor to test for ord

  bool notready; //not ready ?

  infnum primz;//first prime number found, if 0, there's none
  infnum a;//number to test for order

  pthread_mutex_t mymutex;
  pthread_cond_t nextvalid,nextempty;

  int primfakts[primedb<>::n];
  int nfakts;
  

public:
  mythreadvars();
  ~mythreadvars();

private:

  //The ..._v functions are necessary to call 
  //pthread_create

  inline void testord();
  static void* testord_v(void*x);

  inline void geninfnum();

  inline void testprime();
  static void* testprime_v(void*x);

public:
  

  /*
   *initiates multiple threads to search for a randomly selected prime
   *\return a reference on primz, 
   *which is then the newly created prime number 
  */
  infnum& getprime();

};


template<class Tint>
void genfactorizable(Tint&ret,double lnret)
{
  double lnvect[primzahlen.n];
  double powvect[primzahlen.n];
  
  double sumpow=0;
  double integ,fract;
  for(int i=0;i<primzahlen.n;i++){
    powvect[i] = drand48();    
    sumpow += powvect[i];
  }

  ret=1;
  for(int i=0;i<primzahlen.n;i++){   

    powvect[i] *= lnret / sumpow / log(primzahlen.p[i]) ;
    fract = modf(powvect[i],&integ);
    integ+= ( fract > drand48() );

    for(int j=0;j<integ;++j)
      ret*=primzahlen.p[i];

  }
  
}


int genprime(infnum&n);


template<class Tint>
Tint euklid(Tint a,Tint b,Tint&ai,Tint&bi);

template<class Tint>
Tint ggt(Tint a,Tint b);

template<class Tint>
Tint jakobi(Tint a,Tint q);

template<class Tint>
Tint powmod(Tint a,Tint exponent,Tint n);

template<class Tint>
Tint makeEqualDist(const Tint&n);

template<class Tint>
int W(const Tint& a, const Tint& n);

#include "primtest.hpp"


#endif

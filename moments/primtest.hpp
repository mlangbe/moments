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

#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<iostream>
using namespace std;

#include "infnum.hh"

template<class Tint,class myint>
/**
 * division with rest:
 * divide numer by denom.
 *\post: numer contains quotient;denom contains rest
 *\param numer: enumerator;will be replaced by quotient;
 *\param denom: denominator;will be replaced by rest;
 */
void mydiv(Tint&numer,myint&denom)
{
  myint ret = numer%denom;
  numer /= denom;
  denom=ret;
}

/**specializations*/

void mydiv(infnum&numer,infnum::t1elem &denom);

void mydiv(infnum&numer,infnum&denom);

void mydiv(long int &numer,long int&denom);


//generates nn prime numbers
template<int nn>
primedb<nn>::primedb()
{
  p[0]=2;

  int next=1;
  
  for(int i=3;next<n;i+=2){
    
    int j,bound = int(ceil(sqrt(i)));
    
    for( j=0 ; p[j] <= bound  & j<next ; ++j )
      if ( i % p[j]==0 ) { j=-1; break; }
    
    if ( j>=0 ){
      p[next++] = i;
    }
    
  }
}

//tests if n primzahl
//1 = ja, 0 = nein, -1 = weiss nicht

template<int nn>
template<class Tint>
int primedb<nn>::test( const Tint& n )
{
  int bound; 
  double fbound = ceil(sqrt(float(n)));
  bound = int(fbound);
  if(bound < fbound )
    bound = int(~0u>>1);

  ulint *i;
  for (i = begin();i!=end() & *i < bound ;++i){
    if( n % *i == 0 )      
      return 0;
  }

  if( *i >= bound )
    return 1; 
  
  return -1;
}



template<class Tint>
Tint euklid(Tint a,Tint b,Tint&ai,Tint&bi)
{
  Tint f,c[2][3]={{a,1,0},{b,0,1}} ;
  Tint*bigc,*lowc,*xch;

  if(a>b)
    {bigc=c[0];lowc=c[1];}
  else
    {bigc=c[1];lowc=c[0];}

  while( *lowc!=0 ){

    f = bigc[0] / lowc[0];
    bigc[0] -= f * lowc[0];
    bigc[1] -= f * lowc[1];
    bigc[2] -= f * lowc[2];
    assert(bigc[0]== bigc[1]*a + bigc[2]*b);
        

    xch=lowc;lowc=bigc;bigc=xch;    
  }

  ai=bigc[1];bi=bigc[2];

  return bigc[0]; // =ggt
}

template<class Tint>
Tint ggt(Tint a,Tint b)
{
  if(a==0) return b;

  if(a<b)b%=a;

  for(;;){

    if(b==0) return a;

    a%=b;

    if(a==0) return b;

    b%=a;

  }

}


template<class Tint>
Tint jakobi(Tint a,Tint q)
{
 
  //  assert( ggt(a,q) == 1 )

  int ltz=0;

  Tint *pa=&a,*pq=&q,*xch;

  unsigned long ia,iq,m;

  static const unsigned long 
    sh    = sizeof(ia)*8-1,
    sh1   = 1<<sh;

  Tint q0=q;	

  while(*pa!=1)
    {
      //..................................
      //modulo rechnung

      if(*pa>*pq) *pa %= *pq;

      //..................................

      //zweierpotenzen von a abspalten
      //zur effizienzsteigerung in unsigned int rechnen

      //. . . . . . . . . . . . . . . .. 
      // anzahl zweierpotenzen bestimmen und abspalten:

      ia= *pa % sh1; 
      while(!ia)
	{ ia=sh1; mydiv(*pa,ia);  }

      ia|= sh1;//zur terminierung
      iq = 1;	
      m=0;//m=anzahl zweier ungerade?
      while(!(ia&iq))
	{ m^=1 ; iq<<=1;  }

      *pa/=iq;
      //. . . . . . . . . . . . . . . .. 

      //beitrag zu ltz berechnen
      
      iq = *pq % 8; 
      
      ltz ^= ( (iq*iq-1)>>3 ) & m ;

      //................................

      if (*pa==1) break;
      //................................

      //reziprozitaet nutzen      

      ia = *pa % 8;
      ltz ^= ((ia-1)>>1) & ((iq-1)>>1) & 1 ;

      xch=pa;pa=pq;pq=xch;

    }
  
 
  if( ltz )return  q0-Tint(1);
  return 1;

}

template<class Tint>
Tint powmod(Tint a,Tint exponent,Tint n)
{

  Tint ret=1;


  if(exponent<0){
    fprintf(stderr,"powmod:I don't want negative exponents\n");
    throw;
  }

//   if(n>(1LLU<<(sizeof(n)*8/2))){
//     fprintf(stderr,"powmod: n too big:"
// 	    " product of numbers will overflow before mod");
//     throw;
//   }

  if(exponent!=0)for(;;){

    if(a>=n)a%=n;

    unsigned long odd =2;

    mydiv(exponent,odd);

    if (odd){   
      ret *= a;
      if (ret>=n) ret %= n;
    }

    if(exponent==0) break;
   
    a*=a;


  }
    
  return ret;
}

template<class Tint>
Tint makeEqualDist(const Tint&n)
{
  return (Tint(rand())*n/((long unsigned)RAND_MAX) + Tint(rand())) % n ;
}

template<class Tint>
int W(const Tint& a, const Tint& n)
{
  static const Tint my1=1;
  if(ggt<Tint>(a,n)!=1)return 0;  
  return ( powmod<Tint>(a, (n-my1)/2LU, n) == jakobi<Tint>(a,n) ); 
}


template<class Tint>
int pt(const Tint& n,int ntries=2)
{
  static const Tint m1=1,m2=2;
  for(;ntries;--ntries)
    if( ! W(makeEqualDist(n-m2)+m1,n) )
      return 0;

  return 1;
}

//gives the different prime faktors of a and removes them from a
//it is not guarenteed that theprime faktors are complete,
//(only the first ntries numbers are tested)

template<class Tint>
int primfaktors(Tint&a,int*f,int nf,int ntries=1000)
{
  
  int nfakts=0;

  for(ulint*i=primzahlen.begin() ; i!=primzahlen.end() ;i++){
    
    if ( a % *i == 0 ){
      f[nfakts++] = *i;
      do{ a /= *i; } while(a % *i == 0 );
    }

    --ntries;
    if( ntries<=0 | nfakts==nf | a==1 )
      return nfakts;    
  }
  
  int k=primzahlen.n;
  while(gen(gc,k)<=*(primzahlen.end()-1))k++;
 
  for(; a!=1 & k < ntries & nfakts < nf ; k++ ){

    int c = gen(gc,k);
    
    if (  a % c == 0 && primzahlen.test(c)==1 ){
      f[nfakts++] = c;
      do{ a /= (long unsigned)c; } while(a % c == 0 );
    }

  }

  return nfakts;

}

template<class Tint>
Tint phi(Tint a)
{
  Tint phia=1;
  
  for(ulint*i=primzahlen.begin() ; 
      i!=primzahlen.end() & a>Tint(*i)*Tint(*i) ;
      i++){
    
    if ( a % *i == 0 ){
      phia *= *i-1;a/= *i;
      while(a % *i == 0 ){ a /= *i; phia *= *i; } 
    }

  }
  
  Tint k=primzahlen.n;
  while(gen(gc,k)<=Tint(*(primzahlen.end()-1)))++k;
 
  for(Tint c = gen(gc,k) ; a > c*c ; k++, c=gen(gc,k) ){
    
    if ( a % c == 0 ){
      phia *= c-Tint(1); a/= c;
      while(a % c == 0 ){ a /= c; phia *= c; } 
    }

  }

  if(a>1) phia*=a-Tint(1);

  return phia;
}


template<class Tint>
int ord(const Tint&n,const Tint&a,const Tint&n1,int*f,int nf)
{
  static const Tint m1=1;
  int ret=1;
  for(int i=0;i<nf;i++){
    //    cout<<"testing: "<<i<<": a^((n-1)/"<<f[i]<<")==1 (mod n)"<<endl;
    if(powmod(a,(n-m1)/((unsigned long)f[i]),n)==m1)
      {ret=0;break;}
  }
  
  if(ret & n1!=m1)
    ret = (powmod(a,(n-m1)/n1,n)!=m1) ;
  
  return ret;
}



//tests if n prim; 1 = ja, 0 = nein, -1 = weiss nicht
template<class Tint>
int prim(Tint n)
{
  static const Tint my1=1,my2=2;

  //cout<<n<<flush;

  //test by looking for divisors in first 1000 prime numbers
  int d=primzahlen.test(n);

  if(d==1){
    // cout<<"d:no;"<<flush;
    return 1;
  }
  if(d==0){
    //cout<<"[d: yes]"<<flush;
    return 0;
  }

  if(!pt(n))
    {
      //printf("p:no;\n");
      return 0;
    }
  //  else printf("p:yes;");

  //.....................................................................

  //primfaktorzerlegung von n-1 soweit wie moeglich

  int primfakts[primzahlen.n];

  Tint a,n1 = n-my1 ;

  //try the first 10000 numbers 
  int nfakts=primfaktors(n1,primfakts,primzahlen.n,10000);

  if( n1!=my1 && prim(n1)!=1 )
      return -1;        

  //...............................................................
  //test if there is an a with ord_n(a)=n-1

  for(int j=0;j<log(n)*10;j++){

    a = makeEqualDist(n-my2)+my1;

    if(powmod(a,n-my1,n)!=my1)
      { //printf("e:no,j:%d;\n",j);
	return 0;}

    if(ord(n,a,n1,primfakts,nfakts))
      { //printf("e:yes,j:%d;\n",j);
	return 1;}

  }

  printf("dont know\n");

  return -1;

}

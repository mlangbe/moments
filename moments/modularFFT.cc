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

#include<iostream>
#include<assert.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

#include "modularFFT.hh"

#include "infnum.hh"
 
template<class ta,class tb>
ta modularpolynFFTBase::pow(ta a,tb b)
{
  ta ret = b&1? a : ta(1) ;
  for(b>>=1;b;b>>=1){
    a*=a;
    if(b&1)
      ret*=a;
  }
  return ret;
}

const modularpolynFFTBase::myint
  modularpolynFFTBase::myprime;


modularpolynFFTBase::modularpolynFFTBase()
{
  //buffer vars, needed for enough precision
  myrc a,b,c;

  int j=0;

  for(int i=2; i<10000;i++)
    {	

      // a= b^(phi(p)/2^maxpw),
      // so that ord(a) is a power of 2
      b= pow(myrc(i),myprime>>maxpw);
      
      for(a=b,j=0; j<maxpw & a!=myprime-1 & a!=1 ; ++j)
	a*=a;

      if( a==myprime-1 && j==maxpw-1 ) break;

    }

  assert(j==maxpw-1);

  //save b in a;
  a=b;
  for(int i=maxpw;i>=0;--i)
    eiwpw[i]=b,  b*=b;

  assert(eiwpw[1]==myprime-1  && eiwpw[0]==1);
      
  //get inverse of a (mod p) 
  //and write it to b
  b=a.inverse();
  assert( (a*=b , a==1) );
    
  for(int i=maxpw;i>=0;--i)
    eiwpwinv[i]=b, b*=b;

  assert(eiwpwinv[1]==myprime-1 && eiwpwinv[0]==1);
    
  //get inverse of 2 (mod p) 
  //and write it to a
  a=myrc(2).inverse();
  assert( (b=a, b*=2 , b==1) );

  b=a;minv[0]=1;
  for(int i=1;i<=maxpw;++i)
    minv[i]=b, b*=a;

            
}

inline int modularpolynFFTBase::spiegbits(int x,int n)
{
  int ret=0;
  for(;n!=0 & x!=0;n--,x>>=1){ret<<=1;ret|=(x&1);}
  return ret<<n;
}


void modularpolynFFTBase::xchvals(myrc * x, 
				  int pw,
				  bool invert) const
{
  if (invert){
    //c=2^-pw)
    myrc c=minv[pw];
    for(int k=0;!(k>>pw);k++){
      int i=spiegbits(k,pw);
      if(i>=k){
	myrc a=x[i], b=x[k];
	a*=c;b*=c;
	x[i]=b; x[k]=a;}
    }          
  }   
  else
    for(int k=0;!(k>>pw);++k){
      int i=spiegbits(k,pw);
      if(i>=k){myrc a=x[i];x[i]=x[k];x[k]=a;}
    }
}




//here, help functions for FFT_v3 are stored
namespace modular_fft_v3
{

  typedef modularpolynFFTBase::myrc myrc;

  
  //log2( sizeof(L1-Cache)/sizeof(myrc) )
  static const int pcache = 10;


#define INNERLOOP_EIWN(j,d,eiwn)\
myrc *jd = j+d;\
myrc a= *j; a -= *jd; a *= eiwn;\
*j += *jd;\
*jd = a;\
//

//if eiwn==1 (mod m)
#define INNERLOOP_1(j,d)\
myrc *jd = j+d;\
myrc a = *j;a -= *jd; \
*j += *jd;\
*jd = a;\
//


  //makes all the rest of the partitions
  void manyparts(myrc * x,int n,int d,const myrc*peiw)
  {
    myrc *j, *jend ;

    for(;d>1;d>>=1,--peiw){

      int d2 = d<<1;
      //eiwn=1
      for(j=x,jend=j+n ; j!=jend ; j+=d2)
	{INNERLOOP_1(j,d)}
      
      myrc eiw=*peiw,eiwn=eiw;
      for(int k=1 ; k!=d-1 ; ++k, eiwn*=eiw)
	//eiwn= eiw^k
	for(j=x+k,jend=j+n ; j!=jend ; j+=d2)
	  {INNERLOOP_EIWN(j,d,eiwn)}
	  
      //eiwn= eiw^(d-1)
      for(j=x+d-1,jend=j+n ; j!=jend ; j+=d2)
	{INNERLOOP_EIWN(j,d,eiwn)}
    }

    //eiwn=1,d=1
    for(j=x,jend=j+n ; j!=jend ; j+=2)
      {INNERLOOP_1(j,1)}
  }

  //partitions two arrays of length 2*d in two parts of length d
  void twoparts(myrc * x, 
		int d,
		const myrc * peiw
		)
  {
    int d2=d<<1;

    myrc *j;

    //eiwn=1
    j=x;
    {INNERLOOP_1(j,d)}
    j+=d2;
    {INNERLOOP_1(j,d)}

    if(d>1){
      myrc eiw=*peiw,eiwn=eiw;
      myrc *jd2,*jend=x+d-1;
      for(j=x+1,jd2=j+d2 ; j!=jend ; ++j,++jd2,eiwn*=eiw){
	//eiwn= eiw^k
	{INNERLOOP_EIWN(j,d,eiwn)}
	{INNERLOOP_EIWN(jd2,d,eiwn)}
      }

      //eiwn= eiw^(d-1)
      {INNERLOOP_EIWN(j,d,eiwn)}
      {INNERLOOP_EIWN(jd2,d,eiwn)}

      --peiw;d>>=1;

      if( d2 >= (1<<pcache) ){
	twoparts(x+d2,d,peiw);
	twoparts(x,d,peiw);    
      }
      else{
	manyparts(x,d2<<1,d,peiw);
      }
    }
  }


  //partitions array of length 2*d in two parts of length d
  //and invokes subprocedures for further partition
  inline void firstpart(myrc * x, 
			int d,
			const myrc * peiw
			)
  {

    {INNERLOOP_1(x,d)}
    
    if(d>1){
      myrc eiw=*peiw,eiwn=eiw;
      myrc*j=x+1,*jend=x+d-1;
      for(; j!=jend ; ++j, eiwn*=eiw)
	{INNERLOOP_EIWN(j,d,eiwn)}
      
      {INNERLOOP_EIWN(j,d,eiwn)}
      
      twoparts(x,d>>1,peiw-1);
    }
    
  }

#undef INNERLOOP_1
#undef INNERLOOP_EIWN

};//end namespace modular_fft_v3






void modularpolynFFTBase::FFT_v3(myrc * x, 
			      int pw,
 			      bool invert) const
{
  //so that 2^pw fits into int
  assert(pw<sizeof(int)*8-1);
  assert(pw<=maxpw);

  
  modular_fft_v3::
     firstpart(x, 1<<(pw-1), 
	       invert? &eiwpwinv[pw] : &eiwpw[pw]); 

  xchvals(x,pw,invert);
    
}


void modularpolynFFTBase::FFT_v1(myrc * x, 
			      int pw,
 			      bool invert) const
{

  //minimal value for log_2(nx); should be 
  // smaller than size of
  //  1st-level-cache divided by  sizeof(myrc)
  // and should be a power of 2
  static const int pminnx=10;

  //log_2 of number of times an eiw is used.
  //this number  should
  // correspond to half of the associativity
  // of the cache
  static const int pd2mul=1;

  //so that 2^pw fits into int
  assert(pw<sizeof(int)*8-1);
  int n=1<<pw;

  assert(pw<=maxpw);
  //pointer on element with order pw
  const myrc *peiw = invert? &eiwpwinv[pw] : &eiwpw[pw];

  for( int d=n>>1; d; d>>=1, --peiw ){
    int d2=d<<1;
    int nx;

    if(d2 > (n>>pd2mul))
      nx=n;
    else{
      //so that one eiw is used at least 4 times
      //(except first loop)
      nx = d2<<pd2mul;

      //if accesses can still be cache-local,
      //re-use eiwn's more often
      if(nx < (1<<pminnx))
	nx = 1<<pminnx;

      if(nx>n)
	nx=n;
    }
      
    for( myrc*jx=x ; jx != x+n ; jx+=nx ){

      myrc *j, *jend ;

      //eiwn=1
      for(j=jx,jend=j+nx ; j!=jend ; j+=d2){		
	myrc *jd = j+d;      
	myrc a = *j;
	*j += *jd;
	a -= *jd; *jd = a;
      }

      myrc eiw=*peiw,eiwn=eiw;
      for(int k=1 ; k<d-1 ; k++, eiwn*=eiw){
	//eiwn= eiw^k
	for(j=jx+k,jend=j+nx ; j!=jend ; j+=d2){		
	  myrc *jd = j+d;
	  myrc a= *j;
	  *j += *jd;
	  a -= *jd; a *= eiwn;
	  *jd = a;
	}
      }

      //eiwn= eiw^(d-1)
      if(d>1)
	for(j=jx+d-1,jend=j+nx ; j!=jend ; j+=d2){		
	  myrc *jd = j+d;
	  myrc a= *j;
	  *j += *jd;
	  a -= *jd; a *= eiwn;
	  *jd = a;
	} 
    }
  }

  xchvals(x,pw,invert);
    
}




void modularpolynFFTBase::FFT_v2(myrc * x, 
				 int pw,
				 bool invert) const
{
  assert(pw<=maxpw);

  int n=1<<pw;

  //pointer on element with order pw
  const myrc *peiw = invert? &eiwpwinv[pw] : &eiwpw[pw];

  for( int d=n>>1; d; d>>=1, --peiw ){
    int d2=d<<1;
    //i=0
    for(myrc*j=x;j<x+n;j+=d2){		
      myrc *jd = j+d;      
      myrc a = *j;
      *j += *jd;
      a -= *jd; *jd = a;
    }
    myrc eiw=*peiw,eiwn=eiw;
    for(int k=1 ; k<d ; k++){
      for(myrc* j=x+k;j<x+n;j+=d2){		
	myrc *jd = j+d;
	myrc a= *j;
	*j += *jd;
	a -= *jd; a *= eiwn;
	*jd = a;
      }
      if(k<d-1)eiwn*=eiw;
    }
  }

  xchvals(x,pw,invert);    
}


void modularpolynFFTBase::FFT(myrc * x, 
			      int pw,
 			      bool invert) const
{
  FFT_v3(x,pw,invert);
}


const modularpolynFFTBase modularPolynFFT;

int testrestclass()
{
  typedef modularpolynFFTBase::myrc myrc;
  typedef modularpolynFFTBase::myint myint;
  typedef modularpolynFFTBase::myfloat myfloat;
  static const myint myprime=modularpolynFFTBase::myprime;

  myrc A,B,C;
  infnum a,b,c,p;p.fromfloat(myprime);
  

  for(int i=0;i<1000000;++i){

    if(i%1000==0)
      cout<<i<<'\r'<<flush;
    A=myrc(((myint(random())<<30)+random())%myprime);
    B=myrc(((myint(random())<<30)+random())%myprime);
    
    C=A;C*=B;

    a=myfloat(A);b=myfloat(B);

    c = a*b % p;

    //    cout<<"a,b"<<a<<' '<<b<<endl;

    myint ua=A,ub=B,uc=myint(myfloat(c)),aa,bb,diff;
    aa = myint(myfloat(a*b/p));
    bb = myint(myfloat(A)*B/myprime);

    
 //    for(i=-2;i<=2;++i){
//       diff = ua*ub - (bb+i)*myprime - uc ;
//       if(diff==0)
// 	{
// 	  cout<<"yeah:"<<i<<endl;
// 	  break;
// 	}
//     }
//     if(i>2)
//       {cout<<"noyeah"<<endl;return -1;}


//     cout<<"aa,bb"<<aa<<' '<<bb<<endl;

    myint f=myint(C)-myint(myfloat(c));

    if(f!=0){
      cout<<"wrong  !!"<<C<<" - "<<c<<"="<<hex<<f<<dec<<endl;
    }

  }
  



}



int testmodularFFT(int pw)
{
  typedef modularpolynFFTBase::myrc myrc;
  typedef modularpolynFFTBase::myint myint;
  int n=1<<pw;
  myrc m = modularpolynFFTBase::myprime; 
  myrc*a=new myrc[n],*oa=new myrc[n],*i,*j;



  i=a+n,j=oa+n;
  while(i!=a)
    *(--j)=*(--i)=myrc(drand48()*m);


  cout<<"measuring times and checking errors:"<<endl;

  for(int p=pw-4;p<=pw;++p){

    int nn=1<<p;
    int nrep=0;

    if(p<=17)
      nrep = (1<<(22-p))/p;

    float t,tmul=1./CLOCKS_PER_SEC/(2*nrep+1);

    cout<<"log2(n): "<<p<<endl;

    clock_t t1,t2;

    t1=clock();

    for(int ii=0;ii<nrep;++ii){
      modularPolynFFT.FFT_v2(a,p);
      modularPolynFFT.FFT_v2(a,p,true);      
    }
    modularPolynFFT.FFT_v2(a,p);

    t2=clock();
    t=(t2-t1)*tmul;
    
    cout<<(t2-t1)<<"fft_v2 used "<<t
	<<" secs, div by n*log2(n):"
	<<t*1e9/ldexp(p,p)<<" ns"<<endl; 
    
    t1=clock();
    for(int ii=0;ii<nrep;++ii){
      modularPolynFFT.FFT_v1(a,p,true);
      modularPolynFFT.FFT_v1(a,p);      
    }
    modularPolynFFT.FFT_v1(a,p,true);

    t2=clock();
    t=(t2-t1)*tmul;

    cout<<(t2-t1)<<"fft_v1 used "<<t
	<<" secs, div by n*log2(n):"
	<<t*1e9/ldexp(p,p)<<" ns"<<endl; 

    i=a+nn;j=oa+nn;
    while(i!=a){
      --i,--j;
      if(myint(*i)!=myint(*j)){
	cout<<"error FFT2 !"<<endl;throw;
      }
    }

    t1=clock();
    modularPolynFFT.FFT_v3(a,p);
    for(int ii=0;ii<nrep;++ii){
      modularPolynFFT.FFT_v3(a,p,true);
      modularPolynFFT.FFT_v3(a,p);      
    }
    t2=clock();
    t=(t2-t1)*tmul;
    
    cout<<(t2-t1)<<"fft_v3 used "<<t
	<<" secs, div by n*log2(n):"
	<<t*1e9/ldexp(p,p)<<" ns"<<endl; 

    modularPolynFFT.FFT(a,p,true);

    i=a+nn;j=oa+nn;
    while(i!=a){
      --i,--j;
      if(myint(*i)!=myint(*j)){
	cout<<"error FFT3 !"<<endl;throw;
      }
    }

  }



  delete[]a;
  delete[]oa;
}

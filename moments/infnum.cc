/*
	implementations for numbers with infinite number of bits "infnum"

	This file is part of the source code used for the calculation of the moment invariants
	as described in the dissertation of Max Langbein
	https://nbn-resolving.org/urn:nbn:de:hbz:386-kluedo-38558

	Copyright (C) 2020 TU Kaiserslautern, Prof.Dr. Hans Hagen (AG Computergraphik und HCI) (hagen@cs.uni-kl.de), Max Langbein

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
#include<string.h>
#include<stdio.h>
#include<stdlib.h>
#include<ctype.h>
#include<cassert>
#include<iostream>
#include<sstream>
#include<stack>

//needed so that our division algorithm works
//even if the vc++ compiler has option /fp:fast set
#pragma float_control(precise,on)

using namespace std;

#include"infnum.hh"
//#include"modularFFT.hh"

//provides ldexpl,floorrl,frexpl
//#include "mylmath.hh"
#ifdef __sparc
#define ldexpl(a,b) ldexp(a,b)
#define floorl(a) floor(a)
#define frexpl(a,b) frexp(a,b)
#endif

//#ifdef _WIN32_WINNT
#define copy_n(a,n,b) std::copy(a,a+n,b)
//#endif 



int infnum::ndig=10;
long cost=0;


const infnum infnum::zero=0;
const infnum infnum::one=1;
const infnum infnum::two=2;

void infnum::operator<<=(tindex i)
{
  if(!N())return;

  tindex ii=i/bits1elem;
  int im=i-ii*bits1elem;
  sh+=ii;

  if(!im)return;

  if(im<0)
    im+=bits1elem,--sh;

  t1elem top = p.back()>>(bits1elem-im);
  for(iter j=p.end()-1;j!=p.begin();--j)
    *j= ( *j << im ) | ( *(j-1) >> (bits1elem-im) );
  p.front() <<= bits1elem-im;
  if(top)p.push_back(top);

}



void infnum::operator*=(infnum::t1elem i)
{
  if(!N())return;
  t2elem l=i,m=0;
  for(int b=0;b<N();b++){
    m+=l*p[b];
    p[b]=(t1elem)m;
    m>>=bits1elem;
  }
  if(m){resize(N()+1);p[N()-1]=(t1elem)m;}
}



infnum& infnum::operator+=(const infnum&a)
{
  if(this==&a){infnum temp=a;*this+=temp;return*this;}
  if(!a.N())return *this;
  if(!N()){*this=a;return*this;}
  if(sign!=a.sign){
    sign=-sign;(*this)-=a;sign=-sign;
    return *this;
  }
  tindex ud=sh-a.sh,od=a.N()-N()-ud;

  // ud0=ud>0?ud:0, od0=od>0?od:0
  tindex ud0 = ud & -tindex(ud>0),od0= od & -tindex(od>0);

  tindex on=N();
  resize(N()+od0+ud0);

  iter myp;
  citer ap,apend=a.p.end();

  if(ud0){//if lowest t1elem of a is below lowest of this
    //-> ud>0
    //resize this at bottom, move old values ud0 items to top
    copy_backward(p.begin(),p.begin()+on, p.end()-od0);
    //memmove(&p[ud0],&p[0],on*sizeof(p[0]));
    sh = a.sh;

    if(ud0>=a.N()){//if this and a do not overlap:
      copy(a.p.begin(),a.p.end(),p.begin());
      //memcpy(&p[0],&a.p[0],a.N()*sizeof(p[0]));

      fill(p.begin()+a.N(),p.begin()+ud0,0);
      //memset(&p[a.N()],0,(ud0-a.N())*sizeof(p[0]));
      return *this;
    }
    copy_n(a.p.begin(),ud0,p.begin());
    //      memcpy(&p[0],&a.p[0],ud0*sizeof(p[0]));
    myp = p.begin()   + ud0;
    ap  = a.p.begin() + ud0;
  }
  else{//ud<=0:
    myp = p.begin() - ud;
    ap  = a.p.begin();
  }


  if(od0){//if highest element of a is above highest of this    

    if(od0>=a.N()){//if this and a do not overlap
      fill(p.begin()+on,p.end()-a.N(),0);
      copy(a.p.begin(),a.p.end(),p.end()-a.N());
      return *this;
    }
    //fill high elements of this with high elements of a
    apend-=od0;
    copy(apend,a.p.end(),p.end()-od0);
  }


  //a addieren
  t2elem m=0;
  while(ap!=apend){
    m+=*ap;
    m+=*myp;
    *myp=t1elem(m);
    m>>=bits1elem;++myp;++ap;
  }

  //uebertrag durchlaufen lassen
  while(myp!=p.end() & m!=0){
    m+=*myp;
    *myp=t1elem(m);
    m>>=bits1elem;++myp;
  }

  if(m)
    p.push_back((t1elem)m);

  return *this;
}



void infnum::basicSubtract(const infnum&a)
{
  //assumes that |this| > |a|
  //and that signs are equal

  tindex ud=sh-a.sh,b;
  if(ud>0){resizebottom(ud);ud=0;}

  //subtraktion durch addition von ~a +1
  t2elem m=1;
  for(b=-ud;b<a.N()-ud;b++){
    m+=(t1elem)~(a.p[b+ud]);
    m+=p[b];
    p[b]=t1elem(m);
    m>>=bits1elem;
  }
  //uebertrag durchlaufen lassen
  for(;b<N()&&!m;b++){
    m+=(t1elem)(-1);
    m+=p[b];
    p[b]=t1elem(m);
    m>>=bits1elem;
  }

}


infnum& infnum::operator-=(const infnum&a)
{
  //put the arguments into "normal forms"
  corsiz();a.corsiz();
  if(this==&a){resize(0);return *this;}
  if(!a.N())return *this;
  if(!N()){*this=a;sign=-sign;return *this;}
  if(sign!=a.sign){
    //    cout<<" acase ";
    sign=-sign;(*this)+=a;sign=-sign;
    return *this;
  }

  int v=absvgl(*this,a);
  if(v==0){
    resize(0);
    return *this; 
  }
  if(v<0){
    //    cout<<" bcase ";
    infnum temp=a;temp.basicSubtract(*this);
    *this=temp;sign=-sign;
    return *this;
  }

  basicSubtract(a);
  return *this;

}

//--------------------------------------------------------
//basic memory functions:
void infnum::resizebottom(tindex sht)
  //fuegt am Ende sht Leerstellen an
  //(Wenn sht<0 werden -sht Leerstellen weggenommen)
{
  if(-sht>=N())
    {p.resize(0);sh=0;sign=1;return;}
  sh-=sht;
  if(sht<0)
    copy(p.begin()-sht,p.end(),p.begin());
  p.resize(N()+sht);
  if(sht>0){
    copy_backward(p.begin(),p.end()-sht,p.end());
    fill_n(p.begin(),sht,0);
  }
}

//"normalize" number
void infnum::corsiz() const
{
  infnum*th=const_cast<infnum*>(this);

  if(!N()){th->sign=1;th->sh=0;return;}

  citer a = p.end()-1, b=p.begin();

  //get new borders
  while( a!=p.begin() & !*a )--a;
  while( b!=a & !*b )++b;


  if( b == a & !*a ){
    th->p.resize(0);th->sh=0;th->sign=1;
    return;
  }

  //a zeigt jetzt hinter ende 
  //des neuen intervalls
  ++a;

  if(b!=p.begin()){
    copy(b,a,th->p.begin());
    th->sh+=b-p.begin();
  }

  th->p.resize(a-b);

}



void infnum::resize(tindex nsiz)
{
  tindex on=p.size();
  p.resize(nsiz);
  if(nsiz>on)
    fill(p.begin()+on,p.end(),0);
}

void infnum::operator=(const infnum&b)
{
  p=b.p;
  sign=b.sign;
  sh=b.sh;
}


infnum::t1elem&infnum::operator[](tindex i)
{
  i-=sh;
  if(i<0){resizebottom(-i);i=0;}
  if(i>=N())resize(i+1);
  return p[i];
}

//--------------------------------------------------------------------

char*infnum::strxd(char*s,tindex i,int d) const
{
  if(d<=1)return s;
  s[i]=0;
  t1elem d2=1,b2,ul;
  while((ul=d2*d)>d2)d2=ul;
  infnum e=*this;
  tindex b,j,k=0;
  for(j=i-1;j&&e!=infnum();j--){
    b2=e.t1elemdiv(d2);
    e.corsiz();
    for(;b2&&j;j--,b2/=d){
      b=b2%d;
      if(b<10)s[j]=b+'0';
      else s[j]=b-10+'A';
    }
  }
  while(j){s[j]='0';j--;}
  if(i<2)return s;
  while(e!=infnum()){k++;
  b2=e.t1elemdiv(d2);
  e.corsiz();
  copy_n(s+2,i-2,s+1);
  for(;b2&&j;j--,b2/=d){
    b=b2%d;
    if(b<10)s[1]=b+'0';
    else s[1]=b-10+'A';
  }
  }
  if(k>0) {
    k+=i-2;
    copy_n(s+3,i-3,s+2);
    s[2]=',';
    sprintf(s+(i-20),"*%d^%d",d,k);
  }
  if(sign<0)s[0]='-';else s[0]='+';
  return s;
}



void infnum::fromfloat(infnum::tfloat l)
{
  assert(bits2elem==2*bits1elem);
  p.resize(0);sh=0;sign=1;
  if(!l)return;
  if(l<0){sign=-1;l=-l;}

  //if l= infty
  if( 2*l==l & !(l==0) ){
     cerr<<"infty!!"<<l<<endl;
     throw;
   }

  int e;
  frexpl(l,&e);
  resize(sizeof(l)/sizeof(t1elem)+1);
  if(e>0)
    sh = (e-1)/bits1elem -N() +1;
  else
    sh = e/bits1elem -N(); 

  //  printf("sh:%d n:%u l:%Lf",sh,n,l);
  l=ldexpl(l,-(sh+N()-1)*bits1elem);

  if(floorl(l)>(infnum::tfloat)maxt1elem | floorl(l)==0)
  {
   cerr<<"sth wrong:"<<l<<endl;
   throw;
  }
  static const infnum::tfloat lmul=ldexpl(1,bits1elem);
  
  iter j=p.end();
  do{
    --j;
    *j=(t1elem)l;
    //printf("\np[%d]:%u, l:%Lg\n",j-p,*j,l);
    l-=*j;
    l*=lmul;
  }while(j!=p.begin() & l!=0);

  if(j!=p.begin()){
    sh+=j-p.begin();
    p.erase(p.begin(),j);
    //    copy(j,p.end(),p.begin());
    //resize(p.end()-j);
  }

  if(l!=0){
    cerr<<"sth wrong:didn't get all bits of l:"<<l<<endl;
    assert(l==0);
  }
}


//division des vor-komma-teils durch t
//rueckgabewert:rest;
infnum::t1elem infnum::t1elemdiv(t1elem t)
{
  infnum::t2elem m=0;
  iter 
    pa = p.end(),
    paend = (sh < 0 ? p.begin()-sh : p.begin() ) ;

  while(pa>paend ){
    --pa;
    m<<=bits1elem;
    m|=*pa;
    *pa=t1elem(m/t);
    m%=t;
  }

  if( p.size() && *(p.end()-1)==0 )
    p.resize(p.size()-1);

  if( m!=0 & sh>0 ){

    tindex osh=sh; resizebottom(sh);
    pa = p.begin() + osh; paend=p.begin();
    while(pa!=paend){
      --pa;
      m<<=bits1elem;
      *pa=t1elem(m/t);
      m%=t;
    }
  }

  return t1elem(m);
}

infnum::t1elem infnum::operator%(t1elem t) const
{

  int e;

  //if t is a power of 2:
  if( frexp((long double)t,&e) == 0.5 )
    return (*this)[0] & ((1<<(e-1))-1);

  t2elem m=0;
  tindex a,amin=max(-sh,0);

  for(a=N()-1;a>=amin;a--){
    m<<=bits1elem;
    m|=p[a];
    m%=t;
  };

  for(a=sh-1;a>=0&&m;a--){
    m<<=bits1elem;
    m%=t;
  }

  return t1elem(m);
}

void infnum::sqr()
{
  if(!N())return;
  sh+=sh;
  sign=1;
  resize(2*N());
  for(tindex j=N()-1;j>=0;j--){
    tindex k=max(j+1-N()/2,0),jmk=j-k;
    t1elem pj=0;
    while(k<=jmk){
      t2elem sumh=0,suml=0;
      for(t1elem c=infnum::maxt1elem/4-1;k<jmk&&c;k++,jmk--,c--){
	t2elem m=p[jmk];
	m*=p[k];
	suml+=t1elem(m);
	sumh+=m>>bits1elem;
      }
      suml<<=1;sumh<<=1;//gemischte terme kommen immer doppelt vor
      if(k==jmk){
	t2elem m=p[k];
	m*=m;
	suml+=t1elem(m);
	sumh+=m>>bits1elem;
	k++,jmk--;
      }
      suml+=pj;
      pj=t1elem(suml);
      sumh+=suml>>bits1elem;
      for(tindex l=j+1;sumh;l++){ //uebertrag nach oben tragen
	sumh+=p[l];
	p[l]=t1elem(sumh);
	sumh>>=bits1elem;
      }
    }
    p[j]=pj;
  }
}


//in situ version of operator *=
void infnum::muleq(const infnum&i)
{


  if(&i==this){
    infnum tmp=i; *this *= tmp;
    return;
  }
  corsiz();i.corsiz();

  if(!i.N()||!N())return;
  sh+=i.sh;
  sign*=i.sign;
  resize(N()+i.N());
  for(tindex j=N()-1;j>=0;j--){
    tindex k=max(j-N()+i.N()+1,0),//damit k >=0 und j-k <= N()-i.N()-1 
      kend=min(i.N(),j+1);//damit k < i.N()  und k < j+1 und j-k > -1 und 
    t2elem suml=0;
    while(k<kend){//terme auf position j aufsummieren
      t2elem sumh=0;
      // maximal maxt1elem-1 terme auf einmal summieren
      // (damit summe in t2elem hineinpasst)
      for(t1elem c=maxt1elem-1; k<kend & c!=0 ; k++,c-- ){
	t2elem m=((t2elem)p[j-k])*i.p[k];
	suml+=m&maxt1elem;
	sumh+=m>>bits1elem;
      }
      sumh+=suml>>bits1elem;
      suml&=maxt1elem;

      for(tindex l=j+1;sumh;l++){ //uebertrag nach oben tragen
	sumh+=p[l];
	p[l]=(t1elem)sumh;
	sumh>>=bits1elem;
      }
    }
    p[j]=t1elem(suml);
  }
}

void infnum::simplemul(const infnum&a,const infnum&b)
{

  a.corsiz();b.corsiz();

  if(!(a.N()&&b.N()))return;

  t2elem suml=0,sumh=0,sumsumh,m;

  sh=a.sh+b.sh;
  sign=a.sign*b.sign;
  resize(a.N()+b.N());

  for(tindex j=0;j<N();j++){
    tindex 
      kend = min(b.N(),j+1),  //damit k < b.N() und j-k > -1
      k    = max(j-a.N()+1,0);//damit k >= 0    und j-k <= a.N()-1

    sumsumh=0;

    while(k<kend){
      for(t1elem c=maxt1elem-1;k != kend & c != 0;++k,--c){
	m  = a.p[j-k];
	m *= b.p[k];
	suml+=m&maxt1elem;
	sumh+=m>>bits1elem;
      }
      sumh+=suml>>bits1elem;
      suml&=maxt1elem;
      sumsumh+=sumh>>bits1elem;
      sumh&=maxt1elem;
    }
    p[j]=t1elem(suml);
    suml=sumh;
    sumh=sumsumh;
 }
}


//help fcn for fftmul

//fill in every element of a nbits bits from b
template<class outiter,class initer>
inline void extractbits(outiter a,outiter aend,
		   initer b,initer bend,int nbits)
{
  static const int maxbits = sizeof(*b)*8;

  int bitpos=0;
  long belembits=(1l<<nbits)-1;

  for( ; bitpos<maxbits ; ++a,bitpos+=nbits )
    *a = (*b >> bitpos) & belembits;
  bitpos-=maxbits;

  for(++b ; b!=bend ; ++b)
    {
      //fill lowest bits of b into high-order bits of old a
      //because if bitpos>0, a didnt get enough bits from b
      *(a-1) |= (*b & ((1<<bitpos)-1)) << (nbits-bitpos);
      for( ; bitpos<maxbits ; ++a,bitpos+=nbits )
	*a = (*b >> bitpos) & belembits;
      bitpos-=maxbits;
    }

  //fill rest of a with zeros
  fill(a,aend,0);

}

//(a = sum_i { b[i] 2^(nbits*i) }
template<class outiter,class initer>
inline void distributetobits(outiter a,outiter aend,
			initer b,initer bend,int nbits)
{
  static const int maxbits = sizeof(*a)*8;

  unsigned long long sum = 0;
  unsigned long belembits=(1<<nbits)-1;
  int bitpos=0;

  for(; a!=aend & b!=bend ; ++a ){
    *a = (sum & belembits)>>(nbits-bitpos);

    for(; b!=bend & bitpos<maxbits ;++b,bitpos+=nbits){
      sum>>=nbits;      
      sum+= *b;
      *a |= (sum & belembits)<<bitpos;
    }    
    bitpos-=maxbits;
  }

  if(bitpos<0)bitpos+=maxbits;
  sum >>= nbits-bitpos;

  for(; a!=aend  ; ++a,sum>>=maxbits )
    *a = sum;

  assert(sum==0);

}


#ifdef MODULARFFT_HH
void infnum::fftmul(const infnum&x,const infnum&y)
{

  typedef modularpolynFFTBase::myint myint;
  typedef modularpolynFFTBase::myrc myrc;

  static const myint 
    m = modularpolynFFTBase::myprime;

  int belem = bits1elem;//bits per 1 element

  double minxy = x.N()<y.N()?x.N():y.N();


  /* while( (number of elements in sum)
     x (maximum product value)  > m  )  */
  while( ldexp( minxy*bits1elem/belem+1,2*belem ) > m ) 
    --belem;

  assert(belem>0);

  cout<<"using "<<belem<<"bits per entry"<<endl;
  int belembits = (1<<belem)-1;


  //number of elements in output vector
  tindex nelemout = (x.N()+y.N()+belem-1)/belem * bits1elem;

  int pw;
  if(frexp(nelemout,&pw)==0.5)
    --pw;

  cout<<"using fft with 2^"<<pw<<"elements"<<endl;

  tindex n=1<<pw;

  myint *a=new myint[n], *b=new myint[n];
  myrc *rca=(myrc*)a,*rcb=(myrc*)b;

  //extract each belem bits of x to an entry of a

  extractbits(a,a+n,x.p.begin(),x.p.end(),belem);
  modularPolynFFT.FFT(rca,pw);

  extractbits(b,b+n,y.p.begin(),y.p.end(),belem);
  modularPolynFFT.FFT(rcb,pw);

  for(int i=0;i<n;++i)
    rca[i]*=rcb[i];

  delete[]b;

  modularPolynFFT.FFT(rca,pw,true);

  

  p.resize(x.N()+y.N()+2);

  distributetobits(p.begin(),p.end(),a,a+n,belem);

  sh= x.sh+y.sh;
  sign = x.sign*y.sign;
       
}

#endif


void infnum::recursivemul(const infnum&x,const infnum&y)
{  
  x.corsiz();y.corsiz();
 
  assert(x.N()>1 && y.N()>1);

  resize(0);
  p.reserve(x.N()+y.N());
 
  //a should be smaller:
  bool xgt = x.N()>y.N();
  const infnum 
    &a = xgt? y : x, 
    &b = xgt? x : y;
      
  infnum a0,a1,a0pa1,b0,b1,a0b0,a1b1;

  tindex am = a.N()/2,am2 = am<<1;      

  b0.p.reserve(2*am+1);
  a0pa1.p.reserve(a.N()-am+1);
  a0b0.p.reserve(2*am+1);
  a1b1.p.reserve(a.N()+1);

  a0.p.resize(am);
  copy_n(a.p.begin(),am,a0.p.begin());

  a1.p.resize(a.N()-am);
  copy(a.p.begin()+am,a.p.end(),a1.p.begin());

  a0pa1=a0;
  a0pa1+=a1;

  citer bstart=b.p.begin(), bend = bstart + (b.N()/am2)*am2 ;

  for( ; bstart!=bend ; bstart+=am2, this->sh-=am2 ){

    b0.p.resize(am);
    b0.sh=0;
    copy_n(bstart,am,b0.p.begin());

    b1.p.resize(am);
    b1.sh=0;
    copy(bstart+am,bstart+am2,b1.p.begin());

    if(am>=minrecurse){
      a0b0.recursivemul(a0,b0); a1b1.recursivemul(a1,b1);
      b1+=b0;  b0.recursivemul(a0pa1,b1); 
    }
    else{
      a0b0.simplemul(a0,b0); a1b1.simplemul(a1,b1);
      b1+=b0;  b0.simplemul(a0pa1,b1); 
    }

    b0-=a0b0; b0-=a1b1;

    *this+=a0b0;
    b0.addNsh(am); *this+=b0;
    a1b1.addNsh(2*am); *this+=a1b1;    

  }

  if(bstart!=b.p.end()){
    b0.p.resize(b.p.end()-bstart);
    b0.sh=0;
    copy(bstart,b.p.end(),b0.p.begin());    

    if( b0.N() >= minrecurse )
      a0b0.recursivemul(a,b0);
    else
      a0b0.simplemul(a,b0);

    a0b0.addNsh(-a.sh);
    *this+=a0b0;
  }
  
  addNsh( a.sh+b.sh + (bstart-b.p.begin()) );
  sign=a.sign*b.sign;

}




infnum sqr(const infnum&a)
{
	if(!a.N())return infnum::zero;
  infnum ret;
  infnum::t2elem sumh=0,suml=0;
  ret.sh=a.sh*2;
  ret.sign=1;
  ret.resize(2*a.N());
  for(infnum::tindex j=0;j<ret.N();j++){
    infnum::tindex k=max(j+1-a.N(),0),jmk=j-k;
    infnum::t2elem sumsumh=0;
    while(k<=jmk){
      infnum::t2elem sumlx=0,sumhx=0;
      for(infnum::t1elem c=infnum::maxt1elem/4-1;k<jmk&&c;k++,jmk--,c--){
	infnum::t2elem m=a.p[jmk];
	m*=a.p[k];
	sumlx+=infnum::t1elem(m);
	sumhx+=m>>infnum::bits1elem;
      }
      suml+=sumlx<<1;sumh+=sumhx<<1;//gemischte terme kommen immer doppelt vor
      if(k==jmk){
	infnum::t2elem m=a.p[k];
	m*=m;
	suml+=infnum::t1elem(m);
	sumh+=m>>infnum::bits1elem;
	k++,jmk--;
      }
      sumh+=suml>>infnum::bits1elem;
      suml=infnum::t1elem(suml);
      sumsumh+=sumh>>infnum::bits1elem;
      sumh=infnum::t1elem(sumh);
    }
    ret.p[j]=infnum::t1elem(suml);
    suml=sumh;
    sumh=sumsumh;
  }
  return ret;
}

void infnum::cutnround(tindex i)
{
  if(i>=N())return;
  infnum iadd;
  if(p[N()-i-1]>>(bits1elem-1)){
    iadd.resize(1);
    iadd.sh=N()-i+sh;
    iadd.p[0]=1;
  }
  resizebottom(i-N());
  (*this)+=iadd;
}




int absvgl(const infnum&a,const infnum&b) 
{
  a.corsiz();b.corsiz();

  if(!a.N()|!b.N()){
    if(a.N())return 1;
    if(b.N())return -1;
    return 0;
  }

  infnum::tindex i=a.N()+a.sh,j=b.N()+b.sh;
  if(i>j)return 1;
  if(i<j)return -1;

  infnum::citer 
    pa = a.p.end()-1,
    pb = b.p.end()-1,
    pa0= a.p.begin()+(a.N()>b.N())*(a.N()-b.N());
  while(*pa==*pb & pa!=pa0){--pa;--pb;}
  if(*pa>*pb)return 1;
  if(*pa<*pb)return -1;

  if(a.sh>b.sh)return -1;//wenn a weniger
  //"nachkomma"-stellen hat
  if(a.sh<b.sh)return 1;
  return 0;
}


void infnum::npow(t1elem a,t1elem e)
{
  
  static const infnum::tfloat
    maxp=ldexp(1.l,sizeof(infnum::tfloat));

  infnum::tfloat p=a, ret;
  if(e&1)ret=p;
  else ret=1;
  for(e>>=1;e!=0 & p<maxp ;e>>=1){
    p*=p;
    if(e&1){ret*=p;}
  }
//   cout<<"ret:"<<infnum(ret)<<" p="<<infnum(p)<<endl;

  *this=(infnum)ret;
  for(infnum pi=(infnum)p;e;e>>=1){
    //    pi.sqr();
    pi = pi*pi;

    if(e&1)(*this)*=pi;

  }
  corsiz();
}


void infnum::write(ostream&a)
{
  corsiz(); 

  if(sign<0){a.put('-');sign=1;}

  //used to invert sequence of chars 
  stack<char> st;
  char c;
  do{
    c = t1elemdiv(ndig); 
    c+= c>9 ? 'A'-10 : '0' ;
    st.push(c) ;
  }while( Nsh()>0 );

  while(!st.empty())
    {a.put(st.top());st.pop();}
      
  if( vglz() !=0 ){
    a.put('.');
    for(tindex i=a.precision();i!=0 & vglz()!=0;--i){
      *this *= ndig;
      c=(*this)[0];(*this)[0]=0;
      a.put(c + ( c>9 ? 'A'-10 : '0' ) );
    }
  }
}




inline void digout(ostream&a,int dig,int&firstn)
{
  if(firstn|dig){
    firstn=1;
    if(dig>=10)a<<char(dig-10+'A');
    else a<<char(dig+'0');
  }
}

void infnum::writebin(ostream&a,int pd)
{
  if(!N()){a.put('0');return;}
  if(sign<0)a<<'-';
  int shx=int((long(N()+sh)*bits1elem-1)
	      /pd*pd%bits1elem),
    xand=int(long(1)<<pd)-1,
    shrev,if0=0;
  tindex i=N()-1;
  for(;shx>=0;shx-=pd)
    digout(a,(p[i]>>shx)&xand,if0);
  for(i--;i>=0;i--){
    shx+=bits1elem;
    shrev=bits1elem-shx;
    if(shrev<pd){
      digout(a,((p[i]>>shx)|(p[i+1]<<shrev))&xand,if0);
      shx-=pd;
    }
    for(;shx>=0;shx-=pd)
      digout(a,(p[i]>>shx)&xand,if0);
  }
  if(-shx<pd)digout(a,(p[0]<<(-shx))&xand,if0);
  if(sh){
    a<<" * 10^";
    (infnum( ((infnum::tfloat)sh)*bits1elem/pd) ) .writebin(a,pd);
  }
}

ostream&operator<<(ostream&a,infnum x) {
//   if(a.flags()&ios::basefield){
//     if(a.flags()&ios::hex)x.writebin(a,4);
//     else if(a.flags()&ios::oct)x.writebin(a,3);
//     else if(a.flags()&ios::dec){
//       int ondig=infnum::ndig;
//       infnum::ndig=10;
//       x.write(a);
//       infnum::ndig=ondig;}
//     return a;
//   }
//   int pn=1;
//   while(!((infnum::ndig>>pn)&1))pn++;
//   if((1<<pn)==infnum::ndig)
//     x.writebin(a,pn);
//   else {
//     x.write(a);
//   }
  x.write(a);
  return a;
}


istream& operator>>(istream &a,infnum&x) {

  infnum::t1elem xdig;
  char c;
  int exp=0;
  x.resize(0);x.setNsh(0);x.setsign(1);

  if(a.flags()&ios::hex)xdig=16;
  else if(a.flags()&ios::oct)xdig=8;
  else if(a.flags()&ios::dec)xdig=10;
  else xdig=infnum::ndig;

  do{c=a.get();} while(isspace(c));

  if(c=='-') { x.setsign(-1);c=a.get(); }
  if(c=='+') { c=a.get();}  

  while(isalnum(c)){    
    x*=xdig;
    infnum::t1elem dig = (isalpha(c)?toupper(c)-'A'+10:c-'0');    
    if(dig<xdig)
      x+=infnum(dig);
    else
      break;
    c=a.get();
  }

  if(c=='.'){
    c=a.get();
    while(isalnum(c)){    
      --exp;
      x*=xdig;
      infnum::t1elem dig = (isalpha(c)?toupper(c)-'A'+10:c-'0');    

      if(dig<xdig)
				x+=(infnum)dig;
      else
	break;

      c=a.get();
    }    
  }

  if(c=='e'|c=='E'){
    int exp2;
    a>>exp2;
    exp+=exp2;
  }
  else
    a.putback(c);
  
  if(exp){
    infnum e;
    e.npow(xdig,abs(exp));

    if(exp>0){
      x*=e;    
    }
    else{
      infnum::tindex 
	ostnkomm=infnum::stnkomm;
      infnum::stnkomm=-(e.Nsh()+3);
      cout<<"setting stnkomm to"<<x.stnkomm<<endl;
      infnum rest=x;
      rest.quot(e,x);    
      infnum::stnkomm=ostnkomm;
    }
  }

  x.corsiz();

  return a;
}

void infnum::valxd(const char*s,int d)
{  
  infnum::t1elem xdig;
  int exp=0;
  resize(0);setNsh(0);sign=1;

  while(isspace(*s))
    ++s;

  if(*s=='-') { sign=-1;++s; }
  if(*s=='+') { ++s;}  

  while(isalnum(*s)){
    if(*s=='E'|*s=='e' & *(s+1)=='+' | *(s+1)=='-')
      break;
    (*this)*=ndig;
    (*this)+=infnum
      (infnum::t1elem
       (isalpha(*s)?toupper(*s)-'A'+10:*s-'0'));    
    ++s;
  }

  if(*s=='.'){
    ++s;
    while(isalnum(*s)){    
      --exp;
      (*this)*=ndig;
      (*this)+=infnum((t1elem)(isalpha(*s)?toupper(*s)-'A'+10:*s-'0'));    
      ++s;
    }    
  }

  if(*s=='e'|*s=='E'){
    ++s;
    int exp2=0;
    int esign=1;
    if(*s=='-') { esign=-1;++s; }
    if(*s=='+') { ++s;}  
    while(isalnum(*s)){    
      exp2 *= ndig;
      exp2 += isalpha(*s)?toupper(*s)-'A'+10:*s-'0';    
      ++s;
    }        
    exp+=exp2;
  }

  if(exp){
    infnum e;
    e.npow(ndig,abs(exp));

    if(exp>0){
      (*this)*=e;    
    }
    else{
      infnum::tindex 
	ostnkomm=infnum::stnkomm;
      infnum::stnkomm=-(e.Nsh()+3);
      cout<<"setting stnkomm to"<<stnkomm<<endl;
      infnum rest=*this;
      rest.quot(e,*this);    
      infnum::stnkomm=ostnkomm;
    }
  }

  corsiz();

//   string ss(s);
//   istringstream i(ss);
//   int oxd = ndig;
//   ndig = d;
//   i >> *this;
//   ndig = oxd;
}


infnum sqrt(const infnum&a)
{
  static const infnum::t1elem m2=2;
	const infnum& m1=infnum::one;

  infnum ret,d;
  if(a.sgn()<0){cout<<"sqrt of minus";abort();}
  ret=a;
  ret.setNsh(a.Nsh()/2);
  do{
    d=(a/ret-ret)/m2;
    d.corsiz();
    ret+=d;
  }while(d.N()!=0);
  infnum retx2m1=ret*m2-m1;d=ret*ret-a;
  while(d.sgn()>0){d-=retx2m1;retx2m1-=infnum(2);--ret;}
  return ret;
}

#include<math.h>
void sqrt(const infnum&a,infnum&b,infnum&c)//a=b*b*c,b maximal
{
  static const infnum::t1elem m2=2;
	static const infnum m1=infnum((infnum::t1elem)1u);
  infnum sqrta=sqrt(a),g;
  b=m1;c=a;
  infnum::t1elem uf,maxuf=0xff;
  if(infnum((infnum::tfloat)(maxuf*maxuf))>sqrta)maxuf=infnum::t1elem(sqrt((long double)sqrta[0]));
  for(uf=1;uf<maxuf;uf++)
    if((g=c,g.t1elemdiv(uf*uf))==0){c=g;b*=uf;}
  infnum f(uf),fm2p1=f*m2+m1,sqrf=infnum(uf)*uf,m;
  while(sqrf<sqrta){
    m=c;
    m.quot(sqrf,g);
    if(m.vglz()==0){c=g;b*=f;}
    sqrf+=fm2p1;fm2p1+=infnum(2);++f;
  }
}



infnum nuk(infnum::t1elem nn,infnum::t1elem k)
{
  infnum ret=infnum(1);
  infnum::t1elem i;
  if(k>nn/2)k=nn-k;
  for(i=0;i<k;i++){ret*=nn-i;ret/=i+1;}
  return ret;
}

infnum ggt(infnum a,infnum b)
{
  a.corsiz();
  for(;;){
    b.corsiz();
    if(!b.N())return a;
    a%=b;
    a.corsiz();
    if(!a.N())return b;
    b%=a;
  }
}

infnum npow(const infnum&x,infnum::t1elem e)
{  
  if(!e)return infnum(1);
  infnum ret,p=x;
	if(e&1)ret=x;else ret=infnum::one;
  for(e>>=1;e;e>>=1){
    p.sqr();
    p.corsiz();
    if(e&1)ret*=p;
    ret.corsiz();
  }
  return ret;
}

void binsum(infnum&zsum,infnum&tsum, //z�hler,nenner des ergebnisses
	    infnum::t1elem zp,infnum::t1elem tp, //z�hler,nenner der
	    //einzelwahrscheinlichkeit										lichkeit
	    infnum::t1elem nn, //anzahl versuche
	    infnum::t1elem x1,infnum::t1elem x2 //grenzen f�r summe
	    )
{
  if((tp=0)|(zp>tp)|(x2>nn)|(x1>x2))return;
  infnum::t1elem i,zq=tp-zp;
  tsum=npow(infnum(tp),nn);
  infnum pr=npow((infnum)zp,x1);
	pr*=npow((infnum)zq,nn-x1);
  for(i=0;i<x1;i++)
    {pr*=nn-i;pr/=i+1;}
  for(i=x1;i<=x2;i++) {
    zsum+=pr;
    if(0xffff/zp>nn-i)pr*=(nn-i)*zp;
    else {pr*=zp;pr*=nn-i;}
    if(0xffff/zq-1>=i)pr/=(i+1)*zq;
    else{pr/=i+1;pr/=zq;}
  }
  pr=ggt(zsum,tsum);
  zsum/=pr;tsum/=pr;
}

float gettop(infnum&x)
{//gibt oberste stelle an,
 //weitere stellen werden in genauigkeit verbaut
  return x.top(1)+ldexpl(x.top(2),-infnum::bits1elem)+ldexpl(x.top(3),-infnum::bits2elem);
}



 infnum::tindex infnum::getMaxFloatNsh()
{
	tindex ret=1;
	tfloat f= ((t1elem)-1);
	
	while( ldexp(ldexp(f,ret),-ret)==f )
		++ret;

	return ret-2;
}

const infnum::tindex infnum::maxfloatNsh = infnum::getMaxFloatNsh();

infnum::tfloat infnum::tofloat(tindex ex) const
{
  //number of elems needed for representation of tfloat
  static const tindex floatelems 
    =(sizeof(infnum::tfloat)+sizeof(t1elem)-1)
    / sizeof(t1elem)  + 1;

  corsiz();
  tindex e = N() - floatelems;
  e &= -t1elem(e>=0); //equivalent zu: if (e<0) e=0;
  infnum::tfloat ret=0;
  static const infnum::tfloat lmul=ldexpl(1,bits1elem);
  
  citer i=p.end(),iend = p.begin()+e;

  while(i!=iend){
    --i;
    ret = ret*lmul + *i;
  }

  tfloat oret=ret;
  if(ret*2==ret & !(ret==0) )
    {cerr<<"produced infty"<<endl;
    throw;}
  
  ret= ldexpl(ret,(sh+e+ex)*bits1elem)*sign;

  if(ret*2==ret && !(ret==0) )
    {cerr<<"produced infty2 "<<maxfloatNsh<<endl;
    throw;}
  
//   if(e>0 && p[e-1]!=0)
//     {cerr<<"inexact due to limited number of bits in long double"<<e<<endl;}


  return ret;
}

//Anzahl Stellen Ergebnis quot
//wenn negativ, anzahl Nachkommastellen
infnum::tindex infnum::stnkomm=0;


infnum oquot(infnum&z,infnum&t)
{
  infnum ret;
  infnum::tindex tsign=t.sgn(),zsign=z.sgn();
  ret.setsign(zsign/tsign);
  t.setsign(1);z.setsign(1);
  t.corsiz();z.corsiz();
  infnum::tindex shdig;
  if(infnum::stnkomm>0)shdig=infnum::stnkomm-z.Nsh()+t.Nsh();
  else shdig=-infnum::stnkomm;
  t.setNsh(t.Nsh()-shdig);//wenn stnkomm>0, Anzahl stellen ret=stnkomm
  1/t.N(); //abbruch bei division durch null
  if(!z.N())return ret;
  ret.resize(z.Nsh()-t.Nsh()+1);
  float ttop=gettop(t);
  while(z>=t){
    infnum::tindex ridx=z.Nsh()-t.Nsh();
    float xmul=gettop(z)/ttop;
    if(xmul<1){ridx--;xmul=ldexpl(xmul,infnum::bits1elem);}
    if(ridx<0)break;
    infnum::t1elem rnd=(infnum::t1elem)xmul;
    t.addNsh(ridx);
    z-=t*rnd;
    while(z.vglz()<0)
      {z+=t;rnd--;}
    while(z>t)
      {z-=t;rnd++;}
    ret[ridx]=rnd;
    t.addNsh(-ridx);
  }
  t.setsign(tsign);z.setsign(zsign);
  t.addNsh(shdig);ret.addNsh(-shdig);
  return ret;
}



/** 
 * help function for quot
 *\pre : divisor,*this are normalized (by corsiz),
 *
 *\post : (*this) + t*(ret-ret_old) = (*this)_old, (ret-ret_old) is integer ,
 * 0<=*this<divisor
 * XXX_old means here the value of XXX before the method was executed
 *
 *\param *this: contains "Zaehler", after execution contains rest
 *\param t: divisor
 *\retval ret: solution
 */
void infnum::intQuot(const infnum& t,infnum& ret)
{
  tindex maxnsh=2; //maxfloatNsh/4;
  
  infnum::tfloat dl,tl;  
  infnum xdl,xdlt;
  tindex tnsh=t.Nsh(),mynsh=Nsh();  
  tindex mysign=sign;

  tl=t.tofloat(-tnsh);

  //while exponent difference is too big to calculate float directly,
  //use "normalized" floats
  while( mynsh - tnsh > maxnsh  ){    

    dl=this->tofloat(-mynsh);

		tfloat dltl=dl/tl;
		
		if(dltl*2 == dltl)
		{
			cerr<<"infnum::intQuot: produced infty with :"<<dl<<" / "<<tl<<endl;
			throw;
		}

    xdl.fromfloat(dltl);
    xdl.addNsh( mynsh - tnsh );  

    ret += xdl;
    xdlt.mul(xdl,t);  
    (*this)-= xdlt;  

    corsiz();

    mynsh = Nsh();
  }


  while( absvgl(*this,t) >=0 ){

    long double dl2=this->tofloat(-tnsh); 
		long double dl3=dl2/tl;
		dl=floorl(dl3);
		if(dl*2 == dl && dl!=0 )
		{
			cerr<<"infnum::intQuot: produced infty "<<dl<<" with "<<dl2<<"/"<<tl<<" floorl of "<<dl3<<endl;
			throw;
		}


    xdl.fromfloat(dl);
    ret += xdl;
    
    xdlt.mul(xdl,t);
  
    (*this)-= xdlt;
        
  }

  //if the rest hasnt got the same sign as the counter
  if(this->vglz() * mysign < 0 )
    if(t.sign==mysign)
      {*this += t; --ret;}
    else
      {*this -= t; ++ret;}
    
}


/**
 *\pre: divisor!=0
 *\post : *this + divisor*ret == oldthis &&
 * ( stnkomm > 0 && ret.Nsh() == stnkomm ||
 *   stnkomm<= 0 && ret.sh    == stnkomm  )  &&
 *   0 <= |*this| < |divisor| 
 *     *2 ^ (bitst1elem * (stnkomm<=0 ? 
 *                        stnkomm : 
 *                        oldthis.Nsh()-divisor.Nsh()-stnkomm  
 *                       )
 *\param *this: contains "Zaehler", after execution contains Rest
 *\param divisor: divisor
 *\retval ret: solution of division
 */
void infnum::quot(const infnum& divisor,infnum& ret)
{

  infnum&t = *const_cast<infnum*>(&divisor);

  t.corsiz();

  if(!t.N()){
    cerr<<"division by zero"<<endl;
    throw;
  }

  corsiz();

  if(!N()||t.Nsh()>Nsh()) // wenn *this==0 or exponent of t > exponent of *this
    return; 

  //if stnkom<=0, ret has tnkomm digits(=elements) after komma,
  //else stnkomm gives the (approximated+-1) total number of digits for ret;
  //(=the division output)
  tindex sttot = stnkomm>0 ? stnkomm : Nsh() - t.Nsh() - stnkomm;

  tindex tsh=t.sh;

  t.setNsh(Nsh()-sttot);
  
	ret=infnum::zero;
  intQuot(t,ret);

  ret.sh += t.sh-tsh;

  t.sh=tsh;
}





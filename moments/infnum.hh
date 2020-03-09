/*
	numbers with infinite number of bits "infnum"

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

#ifndef __INFNUM_HH

#define __INFNUM_HH

#include<iostream>
#include<vector>
#include<cmath>
#include<cassert>
//#include "modularFFT.hh"

#ifdef showcost
extern long cost;
#endif




class infnum{

public:

  // basic integer type in which the number is stored
  typedef unsigned long t1elem;
  // integer type with double as many bits;
  //used for multiplication of two t1elems 
  typedef unsigned long long t2elem;
  //integer type for indexing,shifting ...
  typedef int tindex;


  static const int bits1elem = sizeof(t1elem)*8;
  static const int bits2elem = sizeof(t2elem)*8;
  static const t1elem maxt1elem = t1elem(-1);

#ifdef __sparc
  //float type to which infnum should be convertible;
  //should be the biggest float type the machine supports
  typedef double tfloat;

  //maximum Nsh so infnum is convertible to tfloat;
  //used in quot
  //static const int maxfloatNsh = (1<<10)/bits1elem;
#else
  typedef long double tfloat;
  //maximum Nsh so infnum is convertible to tfloat
  //used in quot
#endif

  /*maximum Nsh so infnum is convertible to tfloat
	*(= max.exponent /  bitst1elem
  *used in quot*/
	static int getMaxFloatNsh();

	static const int maxfloatNsh;


  //base for ascii conversion (in val,str)
  static int ndig;

  //this var steers the division results:
  //if > 0 total number of elements 
  //else -(number of elements after decimal point)
  static tindex stnkomm;

	//constants for one , zero, and two
	static const infnum zero;
	static const infnum one;
	static const infnum two;

private:
  ///vector holding basic elements
  typedef std::vector<t1elem> container;
  typedef container::iterator iter;
  typedef container::const_iterator citer;

  container p;  

  ///precision ( -(number elements right of .) )
  tindex sh;

  tindex sign;

public:

  inline void setsign(int s);
  inline int sgn() const {return sign;}

  /** empty constructor */
  inline infnum();

  /** constructor from number */
  inline infnum(tfloat);

	/** constructor from (unsigned number */
  inline void fromt1elem(t1elem i)
	{	sh=0;sign=1;if(i){p.resize(1);p[0]=i;} }

	/** constructor from (unsigned number */
  inline void fromint(int i)
	{	
		assert(sizeof(t1elem)==sizeof(int));
		sh=0;
		if(i>0)
		{sign=1;p.resize(1);p[0]=i;} 
		else if(i<0)
		{sign=-1;p.resize(1);p[0]=-i;} 	
	}

	/** constructor from string */
  inline explicit infnum(const char*);

  /**copy constructor*/
  infnum(const infnum&);

  /** destructor */
  inline ~infnum();


  infnum& operator+=(const infnum&);

private:
  //help function for operator -=
  void basicSubtract(const infnum&);
public:
  infnum& operator-=(const infnum&);

  //....................................
  // von +=,-= abgeleitete funktionen

  inline infnum operator-() const
  { infnum ret=*this;ret.sign=-sign; return ret; }

  /** prefix decrement*/
  inline infnum& operator--()
  { return((*this)-=one); }
  
  /** prefix increment*/
  inline infnum& operator++()
  { return((*this)+=one); }

  /** postfix decrement*/
  inline infnum operator -- (int)
  { infnum x=*this;(*this)-=one;return x;}
  
  /** postfix increment*/
  inline infnum operator ++ (int)
  { infnum x=*this;(*this)+=one;return x;}


  //....................................

  /** division;
   *\param t: divisor
   *\param this: ergebnis der division 
   *\return rest
   */
  t1elem t1elemdiv(t1elem t);

  /** modulobildung.
   *\param divisor: teiler
   *\return: abgerundeter rest der division
   */
  t1elem operator %(t1elem divisor) const;

  /** division durch t
   *\pre: *this ist zaehler
   *\post: *this is rest der division, 
   *  ret ergebnis der division,
   *  ret ist integer
   *\param divisor: teiler
   *\retval ret: ergebnis der division
   */
  void intQuot(const infnum&divisor,infnum&ret);

  /** division durch t
   *\pre: *this ist zaehler
   *\post: *this is rest der division, 
   *  ret ergebnis der division,
   * stellenanzahl haengt von stnkomm ab
   *\param divisor: teiler
   *\retval ret: ergebnis der division
   */
  void quot(const infnum&divisor,infnum&ret);


  inline void operator/=(t1elem t)
  {t1elemdiv(t);}

  inline void operator/=(tfloat t)
  {infnum tmp;quot(infnum(t),tmp);*this=tmp;}

  inline void operator/=(const infnum&b)
  {infnum tmp;quot(b,tmp);*this=tmp;}

  inline void operator%=(const infnum&t)
  {infnum tmp;quot(t,tmp);}

  inline void operator%=(t1elem t)
  {infnum tmp(*this % t); *this=tmp;}


  //....................................


  //bitshifts
  void operator>>=(tindex);
  void operator<<=(tindex);

  //private:
  //....................................

  /*help functions for mul*/

  void simplemul(const infnum&,const infnum&);

  //minimum number of elemnets so that recursive mul is applied
  static const int minrecurse=64;
  void recursivemul(const infnum&,const infnum&);

#ifdef MODULARFFT_HH
  /**dont use, inefficient*/
  void fftmul(const infnum&,const infnum&);
#endif

  /** in-situ-multiplication (*this*=b */
  void muleq(const infnum &);


  /** multiplication */
  inline void mul(const infnum&,const infnum&);

  inline infnum operator *(const infnum &b) const
  { infnum ret; ret.simplemul(*this,b); return ret; }

  inline void operator*=(const infnum&x)
  { infnum a;a.simplemul(*this,x);*this = a; }

  void operator*=(t1elem);

  void lmul(t2elem);

  //...............................................
  //basic access ,copy,and resize functions

  void operator=(const infnum&);

  ///exchange two infnums
  inline void swap(infnum&xch);

  inline tindex N() const
  {return p.size();}

  /** number of elements left of dec.point */
  inline tindex Nsh() const
  {return p.size()+sh;}

  inline void setNsh(tindex e) 
  {sh=e-p.size();}

  /** multiply by 2^(e*bits1elem) */
  inline void addNsh(tindex e) 
  {sh+=e;}

  //...................................................

  //access base elements of number
  t1elem&operator[](tindex);
  inline t1elem operator[](tindex) const;

  inline t1elem operator & (t1elem t) const
  { return (*this)[0] & t; }

  //access most significant element of number
  inline t1elem&top(tindex i)
  {return p[N()-1-i];}

  inline t1elem top(tindex i) const
  {return p[N()-1-i];}


  //...................................................

  //convert infnum to float, optionally exponentiate by 
  //2^(bits1elem*ex);
  tfloat tofloat(tindex ex=0) const;

  void fromfloat(tfloat);

  inline operator tfloat() const
  {return tofloat();}


  //...................................................

  void npow(t1elem,t1elem);
  void sqr();
  void cutnround(tindex);

  //...................................................
  //to and from ascii conversion

  /** gibt zahl zur basis ndig aus,
   * vernichtet dabei zahl
   */
  void write(std::ostream&);

  /** gibt zahl zur basis 2^pd aus */
  void writebin(std::ostream&,int pd);

  /**gibt zahl zur basis d mit max.
   *i stellen aus
   *\param s:string
   *\param i:min. anzahl elemente in s
   *\param d:basis der zahldarstellung (dezimal=basis 10)
   *\return: s
   */
  char*strxd(char*s,tindex i,int d) const;

  /**umwandlung in ascii-text
   *\param s: char-array
   *\param i: anzahl elemente in s
   */
  char*str(char*s,int i)
  {return strxd(s,i,ndig);}


  /** liest zahl zur basis d ein
   *\param s:string,der zahl enthaelt
   *\param d: basis der zahldarstellung
   */
  void valxd(const char*s,int d);

  /**liest zahl zur basis 2^pd ein
   *\param s:string,der zahl enthaelt
   *\param d: basis der zahldarstellung
   */
  void valbin(const char*s,int pd);

  /**
   *liest zahl zur basis ndig ein (default:10)
   *\param s: string, der Zahl enthaelt
   */
  inline void val(const char*s)
  { valxd(s,ndig); }
  
  //.........................................................


  /** removes leading null-elements
   * and so makes representation more efficient 
   * this function is declared const,
   * because only the representation changes,not the logic
   */
  void corsiz() const;

  /**aendert groesse*/
  void resize(tindex);

  /*andert anzahl kommastellen(auch ins negative*/
  void resizebottom(tindex);

  //.........................................................

  /** compares two numbers, ignoring the sign
   *\param a, \param b : number to compare
   *\return: <0 if a<b,>0 if a>b, =0 if a==b
   */
  friend int absvgl(const infnum&a,const infnum&b);


  friend inline tfloat log(const infnum&a);

  /** compare this number to zero*/
  int vglz() const;

  /** quadrat */
  friend infnum sqr(const infnum&);

};

inline void infnum::operator>>=(infnum::tindex i)
{ *this <<= -i; }

#ifndef max
inline infnum::tindex max(infnum::tindex a,infnum::tindex b)
{return a>b?a:b;}
inline infnum::tindex min(infnum::tindex a,infnum::tindex b)
{return a<b?a:b;}
#endif

inline void infnum::mul(const infnum&x,const infnum&y)
{

  infnum * mythis = this;

  if( this==&x | this==&y )
    mythis = new infnum;

  if( x.N()<minrecurse | y.N()<minrecurse )
    mythis->simplemul(x,y); 

  else{

#ifdef MODULARFFT_HH

    tindex 
      sumn = x.N()+y.N(),
      maxn=max(x.N(),y.N()),
      minn=sumn-maxn;       
    int lsumn; 
    frexp((long double)sumn,&lsumn);

    if ( maxn*sqrt((long double)minn) < 2.5*sumn*lsumn )
      mythis->recursivemul(x,y);
    else
      mythis->fftmul(x,y);    
#else
    mythis->recursivemul(x,y);
#endif
  }
  
  if(mythis!=this){
    this->swap(*mythis);
    delete mythis;
  }
   
}











inline infnum::tfloat log(const infnum&a)
{
  static const infnum::tfloat lnb 
    = log(2.l)*infnum::bits1elem;
  return log(a.tofloat(-a.Nsh())) + lnb*a.Nsh();
}

inline void infnum::swap(infnum&xch)
{
  p.swap(xch.p);
  tindex x;
  x=xch.sign;xch.sign=sign;sign=x;
  x=xch.sh;xch.sh=sh;sh=x;
}


inline infnum::infnum(tfloat t)
{ fromfloat(t);}

inline infnum::infnum(const infnum&i)
{ *this=i; }


infnum sqrt(const infnum&a);
  
infnum sqr(const infnum&a);
infnum npow(const infnum&x,infnum::t1elem e);

std::ostream&operator << (std::ostream&a,infnum x);
std::istream&operator >> (std::istream&a,infnum&x);

inline void infnum::setsign(int s)
{
  sign = ( ( (s>0) | !N() )<<1 ) - 1;	
}


//inlined functions 
/** compares two numbers
 *\param a, \param b : number to compare
 *\return: <0 if a<b,>0 if a==b, =0 if a==b
 */
inline int vgl(const infnum&a,const infnum&b)
{
  if(a.sgn()!=b.sgn())return a.sgn()-b.sgn();

  return absvgl(a,b)*a.sgn();
}

/** compares number to 0
 *\return: <0 if this<0,  >0 if this >0, =0 if this=0
 */
inline int infnum::vglz() const
{
  corsiz();
  if(!N())return 0;
  return sign;
}


infnum::infnum()
{ p.resize(0);sign=1;sh=0;}

infnum::~infnum()
{
}

/** constructor from string */
infnum::infnum(const char*s)
{
  p.resize(0);
  val(s);
}


inline infnum::t1elem infnum::operator[](tindex i) const
{
  i-=sh;
  if(i<0||i>=N())return 0;
  return p[i];
}

inline infnum operator/(infnum a,const infnum&b)
{infnum tmp;a.quot(b,tmp);return tmp;}


#define INFNUM_BINARY_FROM_UNARY(op,righttype)\
inline infnum operator op(infnum i,righttype t) \
{ i op##= t ; return i; }


INFNUM_BINARY_FROM_UNARY(>>,infnum::tindex);
INFNUM_BINARY_FROM_UNARY(<<,infnum::tindex);
INFNUM_BINARY_FROM_UNARY(*,infnum::t1elem);
INFNUM_BINARY_FROM_UNARY(%,const infnum&);
INFNUM_BINARY_FROM_UNARY(+,const infnum&);
//INFNUM_BINARY_FROM_UNARY(+,infnum::tfloat);
INFNUM_BINARY_FROM_UNARY(-,const infnum&);
//INFNUM_BINARY_FROM_UNARY(-,infnum::tfloat);

#undef INFNUM_BINARY_FROM_UNARY

inline infnum operator/(infnum a,infnum::t1elem b)
{a.t1elemdiv(b);return a;}

inline infnum abs(infnum a)
{a.setsign(1);return a;}

#define INFNUM_VGL_DERIVE(op) \
inline int operator op (const infnum&a,const infnum&b) \
{return vgl(a,b) op 0;}\
inline int operator op (infnum::tfloat a,const infnum&b) \
{return vgl(infnum(a),b) op 0;}\
inline int operator op (const infnum&a,infnum::tfloat b) \
{return vgl(a,infnum(b)) op 0;}\
inline int operator op (int a,const infnum&b) \
{return vgl(infnum((infnum::tfloat)a),b) op 0;}\
inline int operator op (const infnum&a,int b) \
{return vgl(a,infnum((infnum::tfloat)b)) op 0;}\


INFNUM_VGL_DERIVE(<);
INFNUM_VGL_DERIVE(>);
INFNUM_VGL_DERIVE(<=);
INFNUM_VGL_DERIVE(>=);
INFNUM_VGL_DERIVE(!=);
INFNUM_VGL_DERIVE(==);

#undef INFNUM_VGL_DERIVE

#endif

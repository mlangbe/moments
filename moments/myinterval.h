/*
    library for interval arithmetics.

	This file is part of the source code used for the calculation of the moment invariants
	as described in the dissertation of Max Langbein
	https://nbn-resolving.org/urn:nbn:de:hbz:386-kluedo-38558

	Copyright (C) 2020 TU Kaiserslautern,
	              Prof.Dr. Hans Hagen (AG Computergraphik und HCI) (hagen@cs.uni-kl.de)

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

#ifndef GNUCC
#pragma once
#endif

#ifndef MYINTERVAL_H
#define MYINTERVAL_H

#define _USE_MATH_DEFINES
#include <math.h>
/*
* a class for myintervalval arithmetics.
* if floating point unit is set 
* so it rounds down, 
* (and no optimization used)
* outward rounding is implemented.
*/

template<class T>
class myinterval
{

public:
  T a,b;

  typedef T value_type;

  inline myinterval(const T& a,const T&b):a(a),b(b){}
  inline myinterval(const T& a):a(a),b(a){}
  inline myinterval(){}

  template<class S>
  explicit myinterval(const myinterval<S>&x)
  {a=x.a;b=x.b;}

  inline myinterval sqr() const
  { 
    myinterval ret;
    if(a>=0)
    { ret.a=a*a;ret.b=-b;ret.b*=b; ret.b=-ret.b; }
    else if(b>0)
      if(-a>b)
      { ret.a=0;ret.b=-a;ret.b*=a; ret.b=-ret.b; }
      else
      { ret.a=0;ret.b=-b;ret.b*=b; ret.b=-ret.b; }
    else
    { ret.a=b*b;ret.b=-a;ret.b*=a; ret.b=-ret.b; }

    return ret;
  }

  // PI is a function so it grasps the current rounding model.
  inline static myinterval PI()
  { 
    myinterval ret(M_PI,-M_PI);
    ret.b+=T( -1e-35 ) ;
    ret.b= -ret.b;
    return ret; 
  }

  inline static myinterval EMPTY()
  { 
    return myinterval(HUGE_VAL,-HUGE_VAL);
  }

  inline T mid() const {return .5*(a+b); }
  inline T width() const { return -(a-b); }

  inline myinterval operator - () const {return myinterval(-b,-a);}

  inline void intersect(const myinterval&x)
  {
    if( x.a > a )a=x.a;
    if( x.b < b) b=x.b;  
    if(a>b)
    { a = +HUGE_VAL; b=-HUGE_VAL; }
  }

  inline void join(const myinterval&x)
  {
    if( x.a < a )a=x.a;
    if( x.b > b) b=x.b;  
  }

  inline myinterval operator + (const myinterval&x) const
  { return myinterval(x.a+a,-(-x.b-b) ); }

  inline void operator += (const T&x) 
  { a+=x; b=-b-x;b=-b; }

  inline void operator += (const myinterval&x) 
  { a+=x.a; b=-b-x.b;b=-b; }

  inline void operator -= (const myinterval&x) 
  { a-=x.b; b=-b+x.a;b=-b; }

  inline void operator -= (const T&x) 
  { a-=x; b=-b+x;b=-b; }

  inline myinterval operator - (const myinterval&x) const
  { myinterval r;r.a=a-x.b;r.b=-b+x.a;r.b=-r.b;return r; }

  inline myinterval& operator *= (const T&x)
  { 
    if(x>0)
    {a*=x; b*=-x;b=-b;}
    else
    {T newb = -a*x; a = b*x; b =-newb;}
    return *this;
  }

  inline bool contains(T x) const
  {return (a<=x) & (x<=b);}

 
  enum FLAGS{
    POS2=0, NULL2=2, NEG2=3,EMPTY2=1,
    POS1  = POS2<<2,NULL1 = NULL2<<2, NEG1  = NEG2<<2,EMPTY1= EMPTY2<<2
  };

  inline myinterval operator * (const myinterval&x) const
  { 
    int f=((a<=0)<<3)|((b<0)<<2)|((x.a<=0)<<1)|(x.b<0) ;           

    myinterval r;
    switch(f)
    {
    case POS1|POS2:   r.a= a * x.a; r.b=( b * -x.b); r.b= -r.b; return r; 
    case POS1|NULL2:  r.a= b * x.a; r.b=( b * -x.b); r.b= -r.b; return r; 
    case POS1|NEG2:   r.a= b * x.a; r.b=( a * -x.b); r.b= -r.b; return r; 
    case NULL1|POS2:  r.a= a * x.b; r.b=( b * -x.b); r.b= -r.b; return r;
    case NULL1|NULL2: { 
        T ab;
        ab = a*x.b; r.a = b*x.a;
        if(ab<r.a)r.a=ab;
        ab = a*-x.a; 
        r.b = b*-x.b;
        if(ab<r.b) r.b=ab;
        r.b=-r.b;
        return r;
    } 
    case NULL1|NEG2:  r.a= b * x.a ; r.b=( a * -x.a); r.b= -r.b; return r;
    case NEG1|POS2:   r.a= a * x.b ; r.b=( b * -x.a); r.b= -r.b; return r; 
    case NEG1|NULL2:  r.a= a * x.b ; r.b=( a * -x.a); r.b= -r.b; return r;   
    case NEG1|NEG2:   r.a= b * x.b ; r.b=( a * -x.a); r.b= -r.b; return r; 
    default: r.a=HUGE_VAL,r.b=-HUGE_VAL; return r;
      //one or both intervals empty
    }
  }

  inline void operator *= (const myinterval&x)
  {*this = *this * x ;}

  inline myinterval operator / (const T&x) const
  {
    myinterval ret;
    if(x>0)
    {ret.a = a/x;ret.b = b/-x;}
    else
    {ret.a = b/x;ret.b = a/-x;}
    ret.b=-ret.b;
    return ret;
  }

  inline myinterval operator / (const myinterval&x) const
  {     
    int f=((a<=0)<<3)|((b<0)<<2)|((x.a<=0)<<1)|(x.b<0) ;           

    myinterval ret;
    switch(f)
    {
    case POS1|POS2:   ret.a = a / x.b; ret.b = b / -x.a; ret.b=-ret.b; return ret;
    case POS1|NEG2:   ret.a = b / x.b; ret.b = a / -x.a; ret.b=-ret.b; return ret;
    case NULL1|POS2:  ret.a = a / x.a; ret.b = b / -x.a; ret.b=-ret.b; return ret;
    case NULL1|NEG2:  ret.a = b / x.b; ret.b = a / -x.b; ret.b=-ret.b; return ret;
    case NEG1|POS2:   ret.a = a / x.a; ret.b = b / -x.b; ret.b=-ret.b; return ret;
    case NEG1|NEG2:   ret.a = b / x.a; ret.b = a / -x.b; ret.b=-ret.b; return ret;
    case POS1|NULL2:  //division by an interval overlapping zero: infinite interval
    case NEG1|NULL2:
    case NULL1|NULL2: 
      ret.a=-HUGE_VAL;ret.b=HUGE_VAL; return ret;
    default:  // enumerator or denominator are empty intervals:
      ret.a =HUGE_VAL;ret.b=-HUGE_VAL; return ret; //empty interval
    }    
  }

  inline void operator /= (const myinterval&x) 
  { *this= *this / x;}
 


  //normalize angle to interval [-pi,pi];
  //round down (always to the negative direction)
  //(helper function for cosine)
  inline static T normalizeAngle(T rot)
  {    
    myinterval pi = PI();

    //rot >0: subtract as much as possible
    //rot <0: add as less as possible
    if(rot>pi.a)
    {
      pi *= 2; 
      T f = -rot;
      f /= pi.a; 
      f = floor(f + .5); 
      f *= pi.b;
      return rot + f;
    }    
    if(rot<-pi.b)
    {
      pi *= 2;
      T f = -rot;
      f /= pi.b;  
      f = floor(f + .5 );      
      f *= pi.a;
      return rot + f;
    }
    return rot;
  }

  //the minimum of cos([x.a,x.b])
  //(it is assumed x.b-x.a < 2*PI )
  //(helper function for cos)
  inline T minCos() const
  {
    myinterval pi=pi.PI();
    myinterval w( normalizeAngle(a),  - normalizeAngle(-b) );

    if( w.a==1 && w.b==1 )
      return 1;

  // if interval overlaps PI or -PI
    if( w.a>w.b || (w.a < -pi.a ) ||  (pi.a < w.b ) )
      return -1;

    T ret;

    // the one with the greater absolute value (<PI) gives the smaller cosine 
    if(fabs(w.a)>fabs(w.b))
      ret=cos(w.a);
    else
      ret=cos(w.b);

    //printf("\n conv.Interval: %.6g %.6g",w.a/M_PI,w.b/M_PI);
    //make sure it is rounded down  (if floatpoint rounding set to down)
    ret += T(-1e-35) ;   
  
    return ret;
  }

    inline myinterval<bool> operator == (const myinterval&x) const
	{
		if(a==b && x.a==x.b && a==x.a)
			return true;

		if(b<x.a || x.b<a)
			return false;

		return myinterval<bool>(false,true);
	}

   inline myinterval<bool> operator < (const myinterval&x) const
   {
	   if(b<x.a)
		   return true;
	   if(a>x.b)
		   return false;

	   return myinterval<bool>(false,true);
   }


};

template<class T>
 inline myinterval<T> operator * (const T& a,const myinterval<T>&x) 
 {
	myinterval<T> ret(x);
	ret*=a;
	return ret;
 }

inline bool  containsTrue(const myinterval<bool>&x)
{
	return x.b;
}


template<class T>
myinterval<T> mini(const myinterval<T>&a,const myinterval<T>&b)
{
	return myinterval<T>(a.a<b.a?a.a:b.a,  a.b<b.b?a.b:b.b);
}

template<class T>
myinterval<T> maxi(const myinterval<T>&a,const myinterval<T>&b)
{
	return myinterval<T>(a.a>b.a?a.a:b.a,  a.b>b.b?a.b:b.b);
}


template<class T>
myinterval<T> sqrt(const myinterval<T>& x)
{
	return myinterval<T>(sqrt(x.a),sqrt(x.b));
}

template<class T>
myinterval<T> exp(const myinterval<T>& x)
{
	return myinterval<T>(exp(x.a),exp(x.b));
}

template<class T>
myinterval<T> log(const myinterval<T>& x)
{
	return myinterval<T>(log(x.a),log(x.b));
}

template<class T>
myinterval<T> pow(const myinterval<T>& x, const myinterval<T>&ex)
{
	return exp(log(x)*ex);
}

template<class T>
myinterval<T> pow(const myinterval<T>& x, const myinterval<int>&ex)
{
	myinterval<T> ret = pow(x,ex.a);
	for(int i=ex.a+1;i<ex.b;++i)
		ret.join(pow(x,i));
	return ret;
}

template<class T>
myinterval<T> pow(const myinterval<T>& x, int i)
{
	myinterval<T> ret=1,pw=x;
	
	if(i<0)
	{
		i=-i;pw = myinterval<T>(1)/pw;
	}

	for(;;)
	{
		if(i&1)
			ret*=pw;
	
		i>>=1;

		if(!i)
			return ret;

		pw=pw.sqr();
	}

}



template<class T>
myinterval<T> cos(const myinterval<T>& x)
{
  myinterval<T> y, pi = myinterval<T>::PI();
  if(x.width() > pi.a)
  {y.a=-1;y.b=1;}
  else{
    y.a = x.minCos();
    y.b = (pi-x).minCos();
    y.b = -y.b;
  }
  return y;
}

template<class T>
inline myinterval<T> sin(const myinterval<T>& x){
  return cos(x-x.PI()*0.5);
}

template<class T>
inline myinterval<T> atan(const myinterval<T>& y)
{
	return myinterval<T>(atan(y.a),atan(y.b));	
}


template<class T>
inline void joinperiodic( myinterval<T>&a,const myinterval<T>&b, const myinterval<T>& period)
{
	//a contains b
	if(a.a<=b.a&&b.b<=a.b)
		return;

	//b contains a
	if(b.a<=a.a&&a.b<=b.b)
	{
		a=b;
		return;
	}
	
	T dmid=a.mid()-b.mid();

	int nj= (int) floor( dmid / period.b + 0.5 );	
	a-=period*nj;	

	a.join(b);
}


template<class T>
inline myinterval<T> atan2(const myinterval<T>& y,const myinterval<T>& x){

	static const myinterval<T> PI=myinterval<T>::PI(),_2PI=PI*2,PI_2=PI*0.5;

	myinterval<T> ret=myinterval<T>::EMPTY();

	if(x.contains(0) && y.contains(0))
		return myinterval<T>(-PI.b,PI.b);

	myinterval<bool> lt = x<y;
	if(lt.b) //x<y
	{
		myinterval<bool> lt2 = -x<y;
		if(lt2.b) // -x<y (quadrant oben)
		{
			ret=PI_2 -atan(x/y);
		}
		if(!lt2.a)// -x>=y (quadrant links) 
		{
			ret.join( atan(y/x)+PI);
		}
	}
	if(!lt.a) //x>y possible
	{
		myinterval<bool> lt2 = -y<x;
		if(lt2.b)  //x>-y  (quadrant rechts)
		{
			joinperiodic(ret,atan(y/x),_2PI);
		}
		if(!lt2.a) //x<-y (quadrant unten)
		{
			joinperiodic(ret,-PI_2 - atan(x/y),_2PI);
		}	
	}

	if( ret.width() > _2PI.b )
	{
		return myinterval<T>(-PI.b,PI.b);
	}
	

	//shift interval center to 0
	T dmid=ret.mid();
	int nj= (int) floor( dmid / _2PI.b + 0.5 );
	ret-=_2PI*nj;

	return ret;

}

template<class T>
myinterval<T> floor(const myinterval<T>&x)
{
	return myinterval<T>(floor(x.a),floor(x.b));
}


template<class T>
inline myinterval<T> abs(const myinterval<T>& x){
	if(x.a>=0)
		return x;
	if(x.b<=0)
		return -x;

	return myinterval<T>(0, -x.a>x.b?-x.a:x.b);
}


inline myinterval<int> operator % (const myinterval<int>&a,int x)
{
	if(a.width()<x)
	{
		myinterval<int> mod(a.a%x,a.b%x);
		if(mod.a <= mod.b)
			return mod;
	}
	if(a.a>=0)
		return myinterval<int>(0,x-1);
	if(a.b<=0)
		return myinterval<int>(-x+1,0);
	
	return myinterval<int>(-x+1,x-1);
}


template<class T>
inline std::ostream &operator<<(std::ostream&o,const myinterval<T>& x)
{
  o<<"["<<x.a<<","<<x.b<<"]";
  return o;
}


#endif

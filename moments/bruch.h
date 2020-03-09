/*
 	A class implementing a rational number represented with enumerator and denominator.

	This file is part of the source code used for the calculation of the moment invariants
	as described in the dissertation of Max Langbein
	https://nbn-resolving.org/urn:nbn:de:hbz:386-kluedo-38558

	Copyright (C) 2020 TU Kaiserslautern, Prof.Dr. Hans Hagen (AG Computergraphik und HCI)

	Author: Max Langbein

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
#ifndef BRUCH_H
#define BRUCH_H

#ifndef GNUCC
#pragma once
#endif

template<class T,T maxT>
struct bruch{
			public:
			T z,n;
      static inline T abs(T x){return x<0?-x:x;}

			/*
#ifdef GNUCC
  static const T maxT=(T)0x3fffffffffffffffllu;
#else
  static const T maxT= 0x3fffffffffffffff;
#endif
*/
			inline bruch(){z=(T)0.0;n=(T)1.0;}
			inline bruch(long double zz,long double nn){z=(T)zz;n=(T)nn;}
			inline bruch(long double zz){z=(T)zz;n=1;}
			//inline bruch(const T& zz){z=zz;n=1;}
			void operator+=(const bruch&a);
			void operator-=(const bruch&a);
			void operator*=(const bruch&a);
			void operator/=(const bruch&a);

			bruch operator - () { bruch x; x.z=-x.z; return x;}

            bool operator < (const bruch&b) const
            { return z * b.n < b.z * n;  }

            bool operator ! () const 
            { return z==0; }

            bool operator == (long double x) const 
            { return z == T(x)*n; }

						bool operator != (long double x) const 
            { return z != T(x)*n; }

			void kuerze();

			
			template<class TT>
			operator TT () const
			{
				return ((TT)z)/((TT)n);
			}
			

};






template<class T>
inline T ggt(T a,T b)
{
  for(;;){
    if(b==0)return a;
    a%=b;
    if(a==0)return b;
    b%=a;
  }
}


template<class T>
inline T kgv(const T&a,const T&b)
{return a/ggt(a,b)*b;}

template<class T,T maxT>
inline void bruch<T,maxT>::kuerze()
{
T g=ggt(z,n);
if(g<0)g=-g;
z/=g;n/=g;
}


template<class T,T maxT>
inline void bruch<T,maxT>::operator*=(const bruch<T,maxT>&a)
{
 assert( maxT / (abs(z)+1) > abs(a.z)           );
 assert( maxT / abs(n) > abs(a.n)            );
z*=a.z;n*=a.n;
kuerze();
}


template<class T,T maxT>
inline void bruch<T,maxT>::operator/=(const bruch<T,maxT>&a)
{

 assert( maxT / abs(z) > abs(a.n)     );
 assert( maxT / abs(a.z) > abs(n)     );
 z*=a.n;n*=a.z;

 if(n<0)n*=-1,z*=-1;

 kuerze();
}

template<class T,T maxT>
inline void bruch<T,maxT>::operator+=(const bruch<T,maxT>&a)
{
  assert( maxT / abs(a.n) > abs(z) );
  z*=a.n;
  assert( maxT / abs(n) > abs(a.z) );
  z+=a.z*n;
  assert( maxT / abs(a.n) > abs(n) );
  n*=a.n;
  kuerze();
}

template<class T,T maxT>
inline void bruch<T,maxT>::operator -=(const bruch<T,maxT>&a)
{
  z=-z;(*this)+=a;z=-z;
}


#define derivesingle(op) \
template<class T,T maxT> \
inline bruch<T,maxT> operator op (bruch<T,maxT> a, const bruch<T,maxT>&b) \
{	a op ## = b ; return a; }

derivesingle(*)
derivesingle(/)
derivesingle(-)
derivesingle(+)


#define derivevgl(op,nop) \
template<class T,T maxT> \
inline bool operator op (const bruch<T,maxT>& a, const bruch<T,maxT>&b) \
{	return nop; }

derivevgl(>,(b<a));
derivevgl(<=,!(b<a));
derivevgl(>=,!(a<b));


template<class T,T maxT>
inline std::ostream& operator<<(std::ostream&o,const bruch<T,maxT>&b)
{ o<<b.z; if(b.n!=1)o<<'/'<<b.n;return o;}

#endif

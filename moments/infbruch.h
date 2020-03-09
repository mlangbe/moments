/*
	implementations for rational numbers with infinite denominator/enumerator size

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

#ifndef INFBRUCH_H
#define INFBRUCH_H
#pragma once

#include"infnum.hh"
struct infbruch{
			public:
			infnum z,n;
			inline infbruch(){n=infnum::one;}
			inline infbruch(const infnum& zz,const infnum& nn){z=zz;n=nn;}
			inline infbruch(const infnum& zz){z=zz;n=infnum::one;}
			inline infbruch(const infbruch& zz){z=zz.z;n=zz.n;}
			inline infbruch& operator =(const infbruch& zz){z=zz.z;n=zz.n;return *this;}
			inline infbruch& operator =(const infnum& zz){z=zz;n=infnum::one;return *this;}
			inline infbruch& operator =(int x){z.fromint(x);n=infnum::one;return *this;}
			void operator+=(const infbruch&a);
			void operator-=(const infbruch&a);
			void operator*=(const infbruch&a);
			void operator/=(const infbruch&a);

			infbruch operator -()
			{
				infbruch ret=*this;
				ret.z.setsign(-ret.z.sgn());
				return ret;
			}

            bool operator < (const infbruch&b) const
            { return z * b.n < b.z * n;  }

            bool operator ! () const 
            { return z==0; }

            bool operator == (const infbruch& x) const 
            { return z * x.n == x.z * n; }

			bool operator != (const infbruch& x) const 
            { return !(*this==x);; }


			bool operator == (long double x) const 
            { return z == infnum(x)*n; }

			bool operator != (long double x) const 
            { return !(*this==x);; }

			void kuerze();
};

infnum ggt(infnum a,infnum b);


inline std::ostream& operator<<(std::ostream&o,const infbruch&b)
{ o<<b.z; if(b.n!=1)o<<'/'<<b.n;return o;}

#define BINARY_FROM_UNARY(op) \
	inline infbruch operator op (infbruch a, const infbruch& b) \
{ a op##= b; return a; }

BINARY_FROM_UNARY(+)
BINARY_FROM_UNARY(-)
BINARY_FROM_UNARY(*)
BINARY_FROM_UNARY(/)

#undef BINARY_FROM_UNARY


#endif

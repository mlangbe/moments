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


#include"infbruch.h"

infnum kgv(const infnum&a,const infnum&b)
{return a/ggt(a,b)*b;}

void infbruch::kuerze()
{
infnum g=ggt(z,n);
g.setsign(1);
z/=g;n/=g;
}


void infbruch::operator*=(const infbruch&a)
{
z*=a.z;n*=a.n;
kuerze();
}

void infbruch::operator/=(const infbruch&a)
{
a.z.corsiz();
1/a.z.N();//pr�fen auf division durch null
z*=a.n;n*=a.z;

z.setsign(z.sgn()* n.sgn());
n.setsign(1);

kuerze();
}

void infbruch::operator+=(const infbruch&a)
{
z*=a.n;
z+=a.z*n;
n*=a.n;
kuerze();
}

void infbruch::operator -=(const infbruch&a)
{z.setsign(-z.sgn());(*this)+=a;z.setsign(-z.sgn());}


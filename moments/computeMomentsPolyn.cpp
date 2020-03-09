/*
	moment-invariant computation for polynomials, creating the polynomials overt the moment tensor components
	corresponding to the  moment invariant

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

#include "polynom.h"
#include "computeMoments.h"

//we don't want to wait forever!
#pragma optimize("", off)

#define MOMTYP polyn
#include "momenteFunctions.h"
#undef MOMTYP




#if 0 
#include<iostream>
#include "poly.h"

typedef poly<32,int> T;
typedef T mypoly;

inline mypoly operator * (int i,mypoly b)
{ b*=mypoly(i);return b;}



#include "integralPolynomes.h"



void testinteg()
{
	 mypoly A[32];
	for(int i=0;i<32;++i)
		A[i]=mypoly(1,i);

	std::cout
		<<  calcIntegralPoly(A)-calcIntegralPoly2(A)<<std::endl;


}
#endif

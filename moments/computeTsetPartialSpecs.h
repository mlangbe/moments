/*
	specializations for moment-creating functions for specific orders in 3d

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


/** partial specializations of add3dcoords2<N>
*created by create3dcoords2s.
* this file will be included only in pointCloud.cpp
*/

#ifndef OLD_VC
#define _TEMPLATE_ template<>
#else
#define _TEMPLATE_
#endif



_TEMPLATE_
add3dcoords2<1>::add3dcoords2(mytype *A,const mytype *x){
 mytype a=1,b,c;
 b=a;
 c=b;A[0]+=c;c=x[2];A[1]+=c;b=x[1];
 c=b;A[2]+=c;a=x[0];
 b=a;
 c=b;A[3]+=c;
}

_TEMPLATE_
add3dcoords2<2>::add3dcoords2(mytype *A,const mytype *x){
 mytype a=1,b,c;
 b=a;
 c=b;A[0]+=c;c=x[2];A[1]+=c;c*=x[2];A[2]+=c;b=x[1];
 c=b;A[3]+=c;c*=x[2];A[4]+=c;b*=x[1];
 c=b;A[5]+=c;a=x[0];
 b=a;
 c=b;A[6]+=c;c*=x[2];A[7]+=c;b*=x[1];
 c=b;A[8]+=c;a*=x[0];
 b=a;
 c=b;A[9]+=c;
}

_TEMPLATE_
add3dcoords2<3>::add3dcoords2(mytype *A,const mytype *x){
 mytype a=1,b,c;
 b=a;
 c=b;A[0]+=c;c=x[2];A[1]+=c;c*=x[2];A[2]+=c;c*=x[2];A[3]+=c;b=x[1];
 c=b;A[4]+=c;c*=x[2];A[5]+=c;c*=x[2];A[6]+=c;b*=x[1];
 c=b;A[7]+=c;c*=x[2];A[8]+=c;b*=x[1];
 c=b;A[9]+=c;a=x[0];
 b=a;
 c=b;A[10]+=c;c*=x[2];A[11]+=c;c*=x[2];A[12]+=c;b*=x[1];
 c=b;A[13]+=c;c*=x[2];A[14]+=c;b*=x[1];
 c=b;A[15]+=c;a*=x[0];
 b=a;
 c=b;A[16]+=c;c*=x[2];A[17]+=c;b*=x[1];
 c=b;A[18]+=c;a*=x[0];
 b=a;
 c=b;A[19]+=c;
}

_TEMPLATE_
add3dcoords2<4>::add3dcoords2(mytype *A,const mytype *x){
 mytype a=1,b,c;
 b=a;
 c=b;A[0]+=c;c=x[2];A[1]+=c;c*=x[2];A[2]+=c;c*=x[2];A[3]+=c;c*=x[2];A[4]+=c;b=x[1];
 c=b;A[5]+=c;c*=x[2];A[6]+=c;c*=x[2];A[7]+=c;c*=x[2];A[8]+=c;b*=x[1];
 c=b;A[9]+=c;c*=x[2];A[10]+=c;c*=x[2];A[11]+=c;b*=x[1];
 c=b;A[12]+=c;c*=x[2];A[13]+=c;b*=x[1];
 c=b;A[14]+=c;a=x[0];
 b=a;
 c=b;A[15]+=c;c*=x[2];A[16]+=c;c*=x[2];A[17]+=c;c*=x[2];A[18]+=c;b*=x[1];
 c=b;A[19]+=c;c*=x[2];A[20]+=c;c*=x[2];A[21]+=c;b*=x[1];
 c=b;A[22]+=c;c*=x[2];A[23]+=c;b*=x[1];
 c=b;A[24]+=c;a*=x[0];
 b=a;
 c=b;A[25]+=c;c*=x[2];A[26]+=c;c*=x[2];A[27]+=c;b*=x[1];
 c=b;A[28]+=c;c*=x[2];A[29]+=c;b*=x[1];
 c=b;A[30]+=c;a*=x[0];
 b=a;
 c=b;A[31]+=c;c*=x[2];A[32]+=c;b*=x[1];
 c=b;A[33]+=c;a*=x[0];
 b=a;
 c=b;A[34]+=c;
}


#undef _TEMPLATE_

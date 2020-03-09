/*
	a file containing the partial specializations of the constructors of
	translate<N>
	it will be included in pointCloud.cc


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

/* a file containing the partial specializations of the constructors of 
 * translate<N>
 * it will be included in pointCloud.cc
 */


/** partial specializations created by createtranslate() */

#ifndef OLD_VC
#define _TEMPLATE_ template<>
#else
#define _TEMPLATE_
#endif


_TEMPLATE_
translate<1>::translate(mytype*newA, const mytype*A, const mytype *t)
{
mytype at[4];
{
 mytype a=1,b,c;
 b=a;
 c=b;at[0]=c;c=t[2];at[1]=c;b=t[1];
 c=b;at[2]=c;a=t[0];
 b=a;
 c=b;at[3]=c;
}
newA[0]= + A[0];
newA[1]= + A[0]*at[1] + A[1];
newA[2]= + A[0]*at[2] + A[2];
newA[3]= + A[0]*at[3] + A[3];
}

_TEMPLATE_
translate<2>::translate(mytype*newA, const mytype*A, const mytype *t)
{
mytype at[10];
{
 mytype a=1,b,c;
 b=a;
 c=b;at[0]=c;c=t[2];at[1]=c;c*=t[2];at[2]=c;b=t[1];
 c=b;at[3]=c;c*=t[2];at[4]=c;b*=t[1];
 c=b;at[5]=c;a=t[0];
 b=a;
 c=b;at[6]=c;c*=t[2];at[7]=c;b*=t[1];
 c=b;at[8]=c;a*=t[0];
 b=a;
 c=b;at[9]=c;
}
newA[0]= + A[0];
newA[1]= + A[0]*at[1] + A[1];
newA[2]= + A[0]*at[2] + 2*A[1]*at[1] + A[2];
newA[3]= + A[0]*at[3] + A[3];
newA[4]= + A[0]*at[4] + A[1]*at[3] + A[3]*at[1] + A[4];
newA[5]= + A[0]*at[5] + 2*A[3]*at[3] + A[5];
newA[6]= + A[0]*at[6] + A[6];
newA[7]= + A[0]*at[7] + A[1]*at[6] + A[6]*at[1] + A[7];
newA[8]= + A[0]*at[8] + A[3]*at[6] + A[6]*at[3] + A[8];
newA[9]= + A[0]*at[9] + 2*A[6]*at[6] + A[9];
}

_TEMPLATE_
translate<3>::translate(mytype*newA, const mytype*A, const mytype *t)
{
mytype at[20];
{
 mytype a=1,b,c;
 b=a;
 c=b;at[0]=c;c=t[2];at[1]=c;c*=t[2];at[2]=c;c*=t[2];at[3]=c;b=t[1];
 c=b;at[4]=c;c*=t[2];at[5]=c;c*=t[2];at[6]=c;b*=t[1];
 c=b;at[7]=c;c*=t[2];at[8]=c;b*=t[1];
 c=b;at[9]=c;a=t[0];
 b=a;
 c=b;at[10]=c;c*=t[2];at[11]=c;c*=t[2];at[12]=c;b*=t[1];
 c=b;at[13]=c;c*=t[2];at[14]=c;b*=t[1];
 c=b;at[15]=c;a*=t[0];
 b=a;
 c=b;at[16]=c;c*=t[2];at[17]=c;b*=t[1];
 c=b;at[18]=c;a*=t[0];
 b=a;
 c=b;at[19]=c;
}
newA[0]= + A[0];
newA[1]= + A[0]*at[1] + A[1];
newA[2]= + A[0]*at[2] + 2*A[1]*at[1] + A[2];
newA[3]= + A[0]*at[3] + 3*A[1]*at[2] + 3*A[2]*at[1] + A[3];
newA[4]= + A[0]*at[4] + A[4];
newA[5]= + A[0]*at[5] + A[1]*at[4] + A[4]*at[1] + A[5];
newA[6]= + A[0]*at[6] + 2*A[1]*at[5] + A[2]*at[4] + A[4]*at[2] + 2*A[5]*at[1] + A[6];
newA[7]= + A[0]*at[7] + 2*A[4]*at[4] + A[7];
newA[8]= + A[0]*at[8] + A[1]*at[7] + 2*A[4]*at[5] + 2*A[5]*at[4] + A[7]*at[1] + A[8];
newA[9]= + A[0]*at[9] + 3*A[4]*at[7] + 3*A[7]*at[4] + A[9];
newA[10]= + A[0]*at[10] + A[10];
newA[11]= + A[0]*at[11] + A[1]*at[10] + A[10]*at[1] + A[11];
newA[12]= + A[0]*at[12] + 2*A[1]*at[11] + A[2]*at[10] + A[10]*at[2] + 2*A[11]*at[1] + A[12];
newA[13]= + A[0]*at[13] + A[4]*at[10] + A[10]*at[4] + A[13];
newA[14]= + A[0]*at[14] + A[1]*at[13] + A[4]*at[11] + A[5]*at[10] + A[10]*at[5] + A[11]*at[4] + A[13]*at[1] + A[14];
newA[15]= + A[0]*at[15] + 2*A[4]*at[13] + A[7]*at[10] + A[10]*at[7] + 2*A[13]*at[4] + A[15];
newA[16]= + A[0]*at[16] + 2*A[10]*at[10] + A[16];
newA[17]= + A[0]*at[17] + A[1]*at[16] + 2*A[10]*at[11] + 2*A[11]*at[10] + A[16]*at[1] + A[17];
newA[18]= + A[0]*at[18] + A[4]*at[16] + 2*A[10]*at[13] + 2*A[13]*at[10] + A[16]*at[4] + A[18];
newA[19]= + A[0]*at[19] + 3*A[10]*at[16] + 3*A[16]*at[10] + A[19];
}

_TEMPLATE_
translate<4>::translate(mytype*newA, const mytype*A, const mytype *t)
{
mytype at[35];
{
 mytype a=1,b,c;
 b=a;
 c=b;at[0]=c;c=t[2];at[1]=c;c*=t[2];at[2]=c;c*=t[2];at[3]=c;c*=t[2];at[4]=c;b=t[1];
 c=b;at[5]=c;c*=t[2];at[6]=c;c*=t[2];at[7]=c;c*=t[2];at[8]=c;b*=t[1];
 c=b;at[9]=c;c*=t[2];at[10]=c;c*=t[2];at[11]=c;b*=t[1];
 c=b;at[12]=c;c*=t[2];at[13]=c;b*=t[1];
 c=b;at[14]=c;a=t[0];
 b=a;
 c=b;at[15]=c;c*=t[2];at[16]=c;c*=t[2];at[17]=c;c*=t[2];at[18]=c;b*=t[1];
 c=b;at[19]=c;c*=t[2];at[20]=c;c*=t[2];at[21]=c;b*=t[1];
 c=b;at[22]=c;c*=t[2];at[23]=c;b*=t[1];
 c=b;at[24]=c;a*=t[0];
 b=a;
 c=b;at[25]=c;c*=t[2];at[26]=c;c*=t[2];at[27]=c;b*=t[1];
 c=b;at[28]=c;c*=t[2];at[29]=c;b*=t[1];
 c=b;at[30]=c;a*=t[0];
 b=a;
 c=b;at[31]=c;c*=t[2];at[32]=c;b*=t[1];
 c=b;at[33]=c;a*=t[0];
 b=a;
 c=b;at[34]=c;
}
newA[0]= + A[0];
newA[1]= + A[0]*at[1] + A[1];
newA[2]= + A[0]*at[2] + 2*A[1]*at[1] + A[2];
newA[3]= + A[0]*at[3] + 3*A[1]*at[2] + 3*A[2]*at[1] + A[3];
newA[4]= + A[0]*at[4] + 4*A[1]*at[3] + 6*A[2]*at[2] + 4*A[3]*at[1] + A[4];
newA[5]= + A[0]*at[5] + A[5];
newA[6]= + A[0]*at[6] + A[1]*at[5] + A[5]*at[1] + A[6];
newA[7]= + A[0]*at[7] + 2*A[1]*at[6] + A[2]*at[5] + A[5]*at[2] + 2*A[6]*at[1] + A[7];
newA[8]= + A[0]*at[8] + 3*A[1]*at[7] + 3*A[2]*at[6] + A[3]*at[5] + A[5]*at[3] + 3*A[6]*at[2] + 3*A[7]*at[1] + A[8];
newA[9]= + A[0]*at[9] + 2*A[5]*at[5] + A[9];
newA[10]= + A[0]*at[10] + A[1]*at[9] + 2*A[5]*at[6] + 2*A[6]*at[5] + A[9]*at[1] + A[10];
newA[11]= + A[0]*at[11] + 2*A[1]*at[10] + A[2]*at[9] + 2*A[5]*at[7] + 4*A[6]*at[6] + 2*A[7]*at[5] + A[9]*at[2] + 2*A[10]*at[1] + A[11];
newA[12]= + A[0]*at[12] + 3*A[5]*at[9] + 3*A[9]*at[5] + A[12];
newA[13]= + A[0]*at[13] + A[1]*at[12] + 3*A[5]*at[10] + 3*A[6]*at[9] + 3*A[9]*at[6] + 3*A[10]*at[5] + A[12]*at[1] + A[13];
newA[14]= + A[0]*at[14] + 4*A[5]*at[12] + 6*A[9]*at[9] + 4*A[12]*at[5] + A[14];
newA[15]= + A[0]*at[15] + A[15];
newA[16]= + A[0]*at[16] + A[1]*at[15] + A[15]*at[1] + A[16];
newA[17]= + A[0]*at[17] + 2*A[1]*at[16] + A[2]*at[15] + A[15]*at[2] + 2*A[16]*at[1] + A[17];
newA[18]= + A[0]*at[18] + 3*A[1]*at[17] + 3*A[2]*at[16] + A[3]*at[15] + A[15]*at[3] + 3*A[16]*at[2] + 3*A[17]*at[1] + A[18];
newA[19]= + A[0]*at[19] + A[5]*at[15] + A[15]*at[5] + A[19];
newA[20]= + A[0]*at[20] + A[1]*at[19] + A[5]*at[16] + A[6]*at[15] + A[15]*at[6] + A[16]*at[5] + A[19]*at[1] + A[20];
newA[21]= + A[0]*at[21] + 2*A[1]*at[20] + A[2]*at[19] + A[5]*at[17] + 2*A[6]*at[16] + A[7]*at[15] + A[15]*at[7] + 2*A[16]*at[6] + A[17]*at[5] + A[19]*at[2] + 2*A[20]*at[1] + A[21];
newA[22]= + A[0]*at[22] + 2*A[5]*at[19] + A[9]*at[15] + A[15]*at[9] + 2*A[19]*at[5] + A[22];
newA[23]= + A[0]*at[23] + A[1]*at[22] + 2*A[5]*at[20] + 2*A[6]*at[19] + A[9]*at[16] + A[10]*at[15] + A[15]*at[10] + A[16]*at[9] + 2*A[19]*at[6] + 2*A[20]*at[5] + A[22]*at[1] + A[23];
newA[24]= + A[0]*at[24] + 3*A[5]*at[22] + 3*A[9]*at[19] + A[12]*at[15] + A[15]*at[12] + 3*A[19]*at[9] + 3*A[22]*at[5] + A[24];
newA[25]= + A[0]*at[25] + 2*A[15]*at[15] + A[25];
newA[26]= + A[0]*at[26] + A[1]*at[25] + 2*A[15]*at[16] + 2*A[16]*at[15] + A[25]*at[1] + A[26];
newA[27]= + A[0]*at[27] + 2*A[1]*at[26] + A[2]*at[25] + 2*A[15]*at[17] + 4*A[16]*at[16] + 2*A[17]*at[15] + A[25]*at[2] + 2*A[26]*at[1] + A[27];
newA[28]= + A[0]*at[28] + A[5]*at[25] + 2*A[15]*at[19] + 2*A[19]*at[15] + A[25]*at[5] + A[28];
newA[29]= + A[0]*at[29] + A[1]*at[28] + A[5]*at[26] + A[6]*at[25] + 2*A[15]*at[20] + 2*A[16]*at[19] + 2*A[19]*at[16] + 2*A[20]*at[15] + A[25]*at[6] + A[26]*at[5] + A[28]*at[1] + A[29];
newA[30]= + A[0]*at[30] + 2*A[5]*at[28] + A[9]*at[25] + 2*A[15]*at[22] + 4*A[19]*at[19] + 2*A[22]*at[15] + A[25]*at[9] + 2*A[28]*at[5] + A[30];
newA[31]= + A[0]*at[31] + 3*A[15]*at[25] + 3*A[25]*at[15] + A[31];
newA[32]= + A[0]*at[32] + A[1]*at[31] + 3*A[15]*at[26] + 3*A[16]*at[25] + 3*A[25]*at[16] + 3*A[26]*at[15] + A[31]*at[1] + A[32];
newA[33]= + A[0]*at[33] + A[5]*at[31] + 3*A[15]*at[28] + 3*A[19]*at[25] + 3*A[25]*at[19] + 3*A[28]*at[15] + A[31]*at[5] + A[33];
newA[34]= + A[0]*at[34] + 4*A[15]*at[31] + 6*A[25]*at[25] + 4*A[31]*at[15] + A[34];
}
/********************************************************************************/
#undef _TEMPLATE_

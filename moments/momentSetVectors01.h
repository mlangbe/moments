/*
	implementations of invariant computing methods from moment tensors.
	(generated code)

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

//this file is only included in momentSetVectors.cpp.

//MOMTYPE is a define which has to be set in the including class.
//it denotes the omponent type of the tensors
//e.g. polynomes, double values, etc.

#define MYMOMENTSET momentSet<MOMTYPE,27,0,1>

#ifndef OLD_VC
template<>
#endif
const char*const MYMOMENTSET::comment
=
"computation optimized by using the tensor graph\n "
"and splitting it into multiple tensor products and folds.\n"
"here, only one contraction is performed per tensor computation.\n"
"In momentSetVectors02 also multiple contractions are performed";

#ifndef OLD_VC
template<>
#endif
void MYMOMENTSET::compute(valuetype*M,const valuetype*A)
{

//
valuetype a[3];
/* 0*/ a[0]= ( A[20] ) ;
/* 1*/ a[1]= ( A[10] ) ;
/* 2*/ a[2]= ( A[0] ) ;

//
valuetype b[9];
/* 0 0*/ b[0]= ( A[26] ) ;
/* 0 1*/ b[1]= ( A[16] ) ;
/* 0 2*/ b[2]= ( A[6] ) ;
/* 1 0*/ b[3]= ( A[23] ) ;
/* 1 1*/ b[4]= ( A[13] ) ;
/* 1 2*/ b[5]= ( A[3] ) ;
/* 2 0*/ b[6]= ( A[21] ) ;
/* 2 1*/ b[7]= ( A[11] ) ;
/* 2 2*/ b[8]= ( A[1] ) ;

//
valuetype c[18];
/* 0 0 0*/ c[0]= ( A[29] ) ;
/* 0 0 1*/ c[1]= ( A[19] ) ;
/* 0 0 2*/ c[2]= ( A[9] ) ;
/* 1 0 0*/ c[3]= ( A[28] ) ;
/* 1 0 1*/ c[4]= ( A[18] ) ;
/* 1 0 2*/ c[5]= ( A[8] ) ;
/* 2 0 0*/ c[6]= ( A[27] ) ;
/* 2 0 1*/ c[7]= ( A[17] ) ;
/* 2 0 2*/ c[8]= ( A[7] ) ;
/* 1 1 0*/ c[9]= ( A[25] ) ;
/* 1 1 1*/ c[10]= ( A[15] ) ;
/* 1 1 2*/ c[11]= ( A[5] ) ;
/* 2 1 0*/ c[12]= ( A[24] ) ;
/* 2 1 1*/ c[13]= ( A[14] ) ;
/* 2 1 2*/ c[14]= ( A[4] ) ;
/* 2 2 0*/ c[15]= ( A[22] ) ;
/* 2 2 1*/ c[16]= ( A[12] ) ;
/* 2 2 2*/ c[17]= ( A[2] ) ;

/*  sum_{(1 \; 0)\\} ( c )  */
/* ( =  )  */
valuetype d[3];
/* 0*/ d[0]=c[0] + c[9] + c[15];
/* 1*/ d[1]=c[1] + c[10] + c[16];
/* 2*/ d[2]=c[2] + c[11] + c[17];

/*  sum_{(0 \; 3)\\} a c  */
/*  (=sum_{(0 \; 3)\\}()())  */
valuetype f[6];
/* 0 0*/ f[0]=a[0]*c[0] + a[1]*c[1] + a[2]*c[2];
/* 1 0*/ f[1]=a[0]*c[3] + a[1]*c[4] + a[2]*c[5];
/* 2 0*/ f[2]=a[0]*c[6] + a[1]*c[7] + a[2]*c[8];
/* 1 1*/ f[3]=a[0]*c[9] + a[1]*c[10] + a[2]*c[11];
/* 2 1*/ f[4]=a[0]*c[12] + a[1]*c[13] + a[2]*c[14];
/* 2 2*/ f[5]=a[0]*c[15] + a[1]*c[16] + a[2]*c[17];

/*  sum_{(0 \; 2)\\} a b  */
/*  (=sum_{(0 \; 2)\\}()())  */
valuetype h[3];
/* 0*/ h[0]=a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
/* 1*/ h[1]=a[0]*b[3] + a[1]*b[4] + a[2]*b[5];
/* 2*/ h[2]=a[0]*b[6] + a[1]*b[7] + a[2]*b[8];

/*  sum_{(1 \; 4)\\} c c  */
/*  (=sum_{(1 \; 4)\\}()())  */
valuetype i[45];
/* 0 0 0 0*/ i[0]=c[0]*c[0] + c[3]*c[3] + c[6]*c[6];
/* 0 1 0 0*/ i[1]=c[0]*c[1] + c[3]*c[4] + c[6]*c[7];
/* 0 2 0 0*/ i[2]=c[0]*c[2] + c[3]*c[5] + c[6]*c[8];
/* 1 0 0 0*/ i[3]=c[0]*c[3] + c[3]*c[9] + c[6]*c[12];
/* 1 1 0 0*/ i[4]=c[0]*c[4] + c[3]*c[10] + c[6]*c[13];
/* 1 2 0 0*/ i[5]=c[0]*c[5] + c[3]*c[11] + c[6]*c[14];
/* 2 0 0 0*/ i[6]=c[0]*c[6] + c[3]*c[12] + c[6]*c[15];
/* 2 1 0 0*/ i[7]=c[0]*c[7] + c[3]*c[13] + c[6]*c[16];
/* 2 2 0 0*/ i[8]=c[0]*c[8] + c[3]*c[14] + c[6]*c[17];
/* 0 1 0 1*/ i[9]=c[1]*c[1] + c[4]*c[4] + c[7]*c[7];
/* 0 2 0 1*/ i[10]=c[1]*c[2] + c[4]*c[5] + c[7]*c[8];
/* 1 0 0 1*/ i[11]=c[1]*c[3] + c[4]*c[9] + c[7]*c[12];
/* 1 1 0 1*/ i[12]=c[1]*c[4] + c[4]*c[10] + c[7]*c[13];
/* 1 2 0 1*/ i[13]=c[1]*c[5] + c[4]*c[11] + c[7]*c[14];
/* 2 0 0 1*/ i[14]=c[1]*c[6] + c[4]*c[12] + c[7]*c[15];
/* 2 1 0 1*/ i[15]=c[1]*c[7] + c[4]*c[13] + c[7]*c[16];
/* 2 2 0 1*/ i[16]=c[1]*c[8] + c[4]*c[14] + c[7]*c[17];
/* 0 2 0 2*/ i[17]=c[2]*c[2] + c[5]*c[5] + c[8]*c[8];
/* 1 0 0 2*/ i[18]=c[2]*c[3] + c[5]*c[9] + c[8]*c[12];
/* 1 1 0 2*/ i[19]=c[2]*c[4] + c[5]*c[10] + c[8]*c[13];
/* 1 2 0 2*/ i[20]=c[2]*c[5] + c[5]*c[11] + c[8]*c[14];
/* 2 0 0 2*/ i[21]=c[2]*c[6] + c[5]*c[12] + c[8]*c[15];
/* 2 1 0 2*/ i[22]=c[2]*c[7] + c[5]*c[13] + c[8]*c[16];
/* 2 2 0 2*/ i[23]=c[2]*c[8] + c[5]*c[14] + c[8]*c[17];
/* 1 0 1 0*/ i[24]=c[3]*c[3] + c[9]*c[9] + c[12]*c[12];
/* 1 1 1 0*/ i[25]=c[3]*c[4] + c[9]*c[10] + c[12]*c[13];
/* 1 2 1 0*/ i[26]=c[3]*c[5] + c[9]*c[11] + c[12]*c[14];
/* 2 0 1 0*/ i[27]=c[3]*c[6] + c[9]*c[12] + c[12]*c[15];
/* 2 1 1 0*/ i[28]=c[3]*c[7] + c[9]*c[13] + c[12]*c[16];
/* 2 2 1 0*/ i[29]=c[3]*c[8] + c[9]*c[14] + c[12]*c[17];
/* 1 1 1 1*/ i[30]=c[4]*c[4] + c[10]*c[10] + c[13]*c[13];
/* 1 2 1 1*/ i[31]=c[4]*c[5] + c[10]*c[11] + c[13]*c[14];
/* 2 0 1 1*/ i[32]=c[4]*c[6] + c[10]*c[12] + c[13]*c[15];
/* 2 1 1 1*/ i[33]=c[4]*c[7] + c[10]*c[13] + c[13]*c[16];
/* 2 2 1 1*/ i[34]=c[4]*c[8] + c[10]*c[14] + c[13]*c[17];
/* 1 2 1 2*/ i[35]=c[5]*c[5] + c[11]*c[11] + c[14]*c[14];
/* 2 0 1 2*/ i[36]=c[5]*c[6] + c[11]*c[12] + c[14]*c[15];
/* 2 1 1 2*/ i[37]=c[5]*c[7] + c[11]*c[13] + c[14]*c[16];
/* 2 2 1 2*/ i[38]=c[5]*c[8] + c[11]*c[14] + c[14]*c[17];
/* 2 0 2 0*/ i[39]=c[6]*c[6] + c[12]*c[12] + c[15]*c[15];
/* 2 1 2 0*/ i[40]=c[6]*c[7] + c[12]*c[13] + c[15]*c[16];
/* 2 2 2 0*/ i[41]=c[6]*c[8] + c[12]*c[14] + c[15]*c[17];
/* 2 1 2 1*/ i[42]=c[7]*c[7] + c[13]*c[13] + c[16]*c[16];
/* 2 2 2 1*/ i[43]=c[7]*c[8] + c[13]*c[14] + c[16]*c[17];
/* 2 2 2 2*/ i[44]=c[8]*c[8] + c[14]*c[14] + c[17]*c[17];

/*  sum_{(1 \; 2)\\} ( c )  */
/* ( =  )  */
valuetype j[3];
/* 0*/ j[0]=c[0] + c[4] + c[8];
/* 1*/ j[1]=c[3] + c[10] + c[14];
/* 2*/ j[2]=c[6] + c[13] + c[17];

/*  sum_{(1 \; 2)\\} b b  */
/*  (=sum_{(1 \; 2)\\}()())  */
valuetype k[9];
/* 0 0*/ k[0]=b[0]*b[0] + b[1]*b[3] + b[2]*b[6];
/* 0 1*/ k[1]=b[0]*b[1] + b[1]*b[4] + b[2]*b[7];
/* 0 2*/ k[2]=b[0]*b[2] + b[1]*b[5] + b[2]*b[8];
/* 1 0*/ k[3]=b[3]*b[0] + b[4]*b[3] + b[5]*b[6];
/* 1 1*/ k[4]=b[3]*b[1] + b[4]*b[4] + b[5]*b[7];
/* 1 2*/ k[5]=b[3]*b[2] + b[4]*b[5] + b[5]*b[8];
/* 2 0*/ k[6]=b[6]*b[0] + b[7]*b[3] + b[8]*b[6];
/* 2 1*/ k[7]=b[6]*b[1] + b[7]*b[4] + b[8]*b[7];
/* 2 2*/ k[8]=b[6]*b[2] + b[7]*b[5] + b[8]*b[8];

/*  sum_{(2 \; 4)\\} c c  */
/*  (=sum_{(2 \; 4)\\}()())  */
valuetype l[54];
/* 0 0 0 0*/ l[0]=c[0]*c[0] + c[1]*c[3] + c[2]*c[6];
/* 0 0 0 1*/ l[1]=c[0]*c[1] + c[1]*c[4] + c[2]*c[7];
/* 0 0 0 2*/ l[2]=c[0]*c[2] + c[1]*c[5] + c[2]*c[8];
/* 0 0 1 0*/ l[3]=c[0]*c[3] + c[1]*c[9] + c[2]*c[12];
/* 0 0 1 1*/ l[4]=c[0]*c[4] + c[1]*c[10] + c[2]*c[13];
/* 0 0 1 2*/ l[5]=c[0]*c[5] + c[1]*c[11] + c[2]*c[14];
/* 0 0 2 0*/ l[6]=c[0]*c[6] + c[1]*c[12] + c[2]*c[15];
/* 0 0 2 1*/ l[7]=c[0]*c[7] + c[1]*c[13] + c[2]*c[16];
/* 0 0 2 2*/ l[8]=c[0]*c[8] + c[1]*c[14] + c[2]*c[17];
/* 1 0 0 0*/ l[9]=c[3]*c[0] + c[4]*c[3] + c[5]*c[6];
/* 1 0 0 1*/ l[10]=c[3]*c[1] + c[4]*c[4] + c[5]*c[7];
/* 1 0 0 2*/ l[11]=c[3]*c[2] + c[4]*c[5] + c[5]*c[8];
/* 1 0 1 0*/ l[12]=c[3]*c[3] + c[4]*c[9] + c[5]*c[12];
/* 1 0 1 1*/ l[13]=c[3]*c[4] + c[4]*c[10] + c[5]*c[13];
/* 1 0 1 2*/ l[14]=c[3]*c[5] + c[4]*c[11] + c[5]*c[14];
/* 1 0 2 0*/ l[15]=c[3]*c[6] + c[4]*c[12] + c[5]*c[15];
/* 1 0 2 1*/ l[16]=c[3]*c[7] + c[4]*c[13] + c[5]*c[16];
/* 1 0 2 2*/ l[17]=c[3]*c[8] + c[4]*c[14] + c[5]*c[17];
/* 2 0 0 0*/ l[18]=c[6]*c[0] + c[7]*c[3] + c[8]*c[6];
/* 2 0 0 1*/ l[19]=c[6]*c[1] + c[7]*c[4] + c[8]*c[7];
/* 2 0 0 2*/ l[20]=c[6]*c[2] + c[7]*c[5] + c[8]*c[8];
/* 2 0 1 0*/ l[21]=c[6]*c[3] + c[7]*c[9] + c[8]*c[12];
/* 2 0 1 1*/ l[22]=c[6]*c[4] + c[7]*c[10] + c[8]*c[13];
/* 2 0 1 2*/ l[23]=c[6]*c[5] + c[7]*c[11] + c[8]*c[14];
/* 2 0 2 0*/ l[24]=c[6]*c[6] + c[7]*c[12] + c[8]*c[15];
/* 2 0 2 1*/ l[25]=c[6]*c[7] + c[7]*c[13] + c[8]*c[16];
/* 2 0 2 2*/ l[26]=c[6]*c[8] + c[7]*c[14] + c[8]*c[17];
/* 1 1 0 0*/ l[27]=c[9]*c[0] + c[10]*c[3] + c[11]*c[6];
/* 1 1 0 1*/ l[28]=c[9]*c[1] + c[10]*c[4] + c[11]*c[7];
/* 1 1 0 2*/ l[29]=c[9]*c[2] + c[10]*c[5] + c[11]*c[8];
/* 1 1 1 0*/ l[30]=c[9]*c[3] + c[10]*c[9] + c[11]*c[12];
/* 1 1 1 1*/ l[31]=c[9]*c[4] + c[10]*c[10] + c[11]*c[13];
/* 1 1 1 2*/ l[32]=c[9]*c[5] + c[10]*c[11] + c[11]*c[14];
/* 1 1 2 0*/ l[33]=c[9]*c[6] + c[10]*c[12] + c[11]*c[15];
/* 1 1 2 1*/ l[34]=c[9]*c[7] + c[10]*c[13] + c[11]*c[16];
/* 1 1 2 2*/ l[35]=c[9]*c[8] + c[10]*c[14] + c[11]*c[17];
/* 2 1 0 0*/ l[36]=c[12]*c[0] + c[13]*c[3] + c[14]*c[6];
/* 2 1 0 1*/ l[37]=c[12]*c[1] + c[13]*c[4] + c[14]*c[7];
/* 2 1 0 2*/ l[38]=c[12]*c[2] + c[13]*c[5] + c[14]*c[8];
/* 2 1 1 0*/ l[39]=c[12]*c[3] + c[13]*c[9] + c[14]*c[12];
/* 2 1 1 1*/ l[40]=c[12]*c[4] + c[13]*c[10] + c[14]*c[13];
/* 2 1 1 2*/ l[41]=c[12]*c[5] + c[13]*c[11] + c[14]*c[14];
/* 2 1 2 0*/ l[42]=c[12]*c[6] + c[13]*c[12] + c[14]*c[15];
/* 2 1 2 1*/ l[43]=c[12]*c[7] + c[13]*c[13] + c[14]*c[16];
/* 2 1 2 2*/ l[44]=c[12]*c[8] + c[13]*c[14] + c[14]*c[17];
/* 2 2 0 0*/ l[45]=c[15]*c[0] + c[16]*c[3] + c[17]*c[6];
/* 2 2 0 1*/ l[46]=c[15]*c[1] + c[16]*c[4] + c[17]*c[7];
/* 2 2 0 2*/ l[47]=c[15]*c[2] + c[16]*c[5] + c[17]*c[8];
/* 2 2 1 0*/ l[48]=c[15]*c[3] + c[16]*c[9] + c[17]*c[12];
/* 2 2 1 1*/ l[49]=c[15]*c[4] + c[16]*c[10] + c[17]*c[13];
/* 2 2 1 2*/ l[50]=c[15]*c[5] + c[16]*c[11] + c[17]*c[14];
/* 2 2 2 0*/ l[51]=c[15]*c[6] + c[16]*c[12] + c[17]*c[15];
/* 2 2 2 1*/ l[52]=c[15]*c[7] + c[16]*c[13] + c[17]*c[16];
/* 2 2 2 2*/ l[53]=c[15]*c[8] + c[16]*c[14] + c[17]*c[17];

/*  sum_{(1 \; 3)\\} b b  */
/*  (=sum_{(1 \; 3)\\}()())  */
valuetype m[6];
/* 0 0*/ m[0]=b[0]*b[0] + b[1]*b[1] + b[2]*b[2];
/* 1 0*/ m[1]=b[0]*b[3] + b[1]*b[4] + b[2]*b[5];
/* 2 0*/ m[2]=b[0]*b[6] + b[1]*b[7] + b[2]*b[8];
/* 1 1*/ m[3]=b[3]*b[3] + b[4]*b[4] + b[5]*b[5];
/* 2 1*/ m[4]=b[3]*b[6] + b[4]*b[7] + b[5]*b[8];
/* 2 2*/ m[5]=b[6]*b[6] + b[7]*b[7] + b[8]*b[8];

/*  sum_{(2 \; 5)\\} c c  */
/*  (=sum_{(2 \; 5)\\}()())  */
valuetype n[21];
/* 0 0 0 0*/ n[0]=c[0]*c[0] + c[1]*c[1] + c[2]*c[2];
/* 1 0 0 0*/ n[1]=c[0]*c[3] + c[1]*c[4] + c[2]*c[5];
/* 2 0 0 0*/ n[2]=c[0]*c[6] + c[1]*c[7] + c[2]*c[8];
/* 1 1 0 0*/ n[3]=c[0]*c[9] + c[1]*c[10] + c[2]*c[11];
/* 2 1 0 0*/ n[4]=c[0]*c[12] + c[1]*c[13] + c[2]*c[14];
/* 2 2 0 0*/ n[5]=c[0]*c[15] + c[1]*c[16] + c[2]*c[17];
/* 1 0 1 0*/ n[6]=c[3]*c[3] + c[4]*c[4] + c[5]*c[5];
/* 2 0 1 0*/ n[7]=c[3]*c[6] + c[4]*c[7] + c[5]*c[8];
/* 1 1 1 0*/ n[8]=c[3]*c[9] + c[4]*c[10] + c[5]*c[11];
/* 2 1 1 0*/ n[9]=c[3]*c[12] + c[4]*c[13] + c[5]*c[14];
/* 2 2 1 0*/ n[10]=c[3]*c[15] + c[4]*c[16] + c[5]*c[17];
/* 2 0 2 0*/ n[11]=c[6]*c[6] + c[7]*c[7] + c[8]*c[8];
/* 2 0 1 1*/ n[12]=c[6]*c[9] + c[7]*c[10] + c[8]*c[11];
/* 2 1 2 0*/ n[13]=c[6]*c[12] + c[7]*c[13] + c[8]*c[14];
/* 2 2 2 0*/ n[14]=c[6]*c[15] + c[7]*c[16] + c[8]*c[17];
/* 1 1 1 1*/ n[15]=c[9]*c[9] + c[10]*c[10] + c[11]*c[11];
/* 2 1 1 1*/ n[16]=c[9]*c[12] + c[10]*c[13] + c[11]*c[14];
/* 2 2 1 1*/ n[17]=c[9]*c[15] + c[10]*c[16] + c[11]*c[17];
/* 2 1 2 1*/ n[18]=c[12]*c[12] + c[13]*c[13] + c[14]*c[14];
/* 2 2 2 1*/ n[19]=c[12]*c[15] + c[13]*c[16] + c[14]*c[17];
/* 2 2 2 2*/ n[20]=c[15]*c[15] + c[16]*c[16] + c[17]*c[17];

/*  sum_{(1 \; 2)\\} f d  */
/*  (=sum_{(1 \; 2)\\}(sum_{(0 \; 3)\\}()())())  */
valuetype o[3];
/* 0*/ o[0]=f[0]*d[0] + f[1]*d[1] + f[2]*d[2];
/* 1*/ o[1]=f[1]*d[0] + f[3]*d[1] + f[4]*d[2];
/* 2*/ o[2]=f[2]*d[0] + f[4]*d[1] + f[5]*d[2];

/*  sum_{(0 \; 2)\\} ( i )  */
/* ( =  )  */
valuetype q[6];
/* 0 0*/ q[0]=i[0] + i[24] + i[39];
/* 1 0*/ q[1]=i[1] + i[25] + i[40];
/* 2 0*/ q[2]=i[2] + i[26] + i[41];
/* 1 1*/ q[3]=i[9] + i[30] + i[42];
/* 2 1*/ q[4]=i[10] + i[31] + i[43];
/* 2 2*/ q[5]=i[17] + i[35] + i[44];

/*  sum_{(0 \; 2)\\} b d  */
/*  (=sum_{(0 \; 2)\\}()())  */
valuetype s[3];
/* 0*/ s[0]=b[0]*d[0] + b[3]*d[1] + b[6]*d[2];
/* 1*/ s[1]=b[1]*d[0] + b[4]*d[1] + b[7]*d[2];
/* 2*/ s[2]=b[2]*d[0] + b[5]*d[1] + b[8]*d[2];

/*  sum_{(1 \; 3)\\} ( n )  */
/* ( =  )  */
valuetype t[6];
/* 0 0*/ t[0]=n[0] + n[6] + n[11];
/* 1 0*/ t[1]=n[1] + n[8] + n[13];
/* 2 0*/ t[2]=n[2] + n[9] + n[14];
/* 1 1*/ t[3]=n[6] + n[15] + n[18];
/* 2 1*/ t[4]=n[7] + n[16] + n[19];
/* 2 2*/ t[5]=n[11] + n[18] + n[20];

/*  sum_{(1 \; 2)\\} b d  */
/*  (=sum_{(1 \; 2)\\}()())  */
valuetype u[3];
/* 0*/ u[0]=b[0]*d[0] + b[1]*d[1] + b[2]*d[2];
/* 1*/ u[1]=b[3]*d[0] + b[4]*d[1] + b[5]*d[2];
/* 2*/ u[2]=b[6]*d[0] + b[7]*d[1] + b[8]*d[2];

/*  sum_{(1 \; 0)\\} ( n )  */
/* ( =  )  */
valuetype v[6];
/* 0 0*/ v[0]=n[0] + n[3] + n[5];
/* 1 0*/ v[1]=n[1] + n[8] + n[10];
/* 2 0*/ v[2]=n[2] + n[12] + n[14];
/* 1 1*/ v[3]=n[3] + n[15] + n[17];
/* 2 1*/ v[4]=n[4] + n[16] + n[19];
/* 2 2*/ v[5]=n[5] + n[17] + n[20];

/*  sum_{(1 \; 2)\\} b k  */
/*  (=sum_{(1 \; 2)\\}()(sum_{(1 \; 2)\\}()()))  */
valuetype x[9];
/* 0 0*/ x[0]=b[0]*k[0] + b[1]*k[3] + b[2]*k[6];
/* 0 1*/ x[1]=b[0]*k[1] + b[1]*k[4] + b[2]*k[7];
/* 0 2*/ x[2]=b[0]*k[2] + b[1]*k[5] + b[2]*k[8];
/* 1 0*/ x[3]=b[3]*k[0] + b[4]*k[3] + b[5]*k[6];
/* 1 1*/ x[4]=b[3]*k[1] + b[4]*k[4] + b[5]*k[7];
/* 1 2*/ x[5]=b[3]*k[2] + b[4]*k[5] + b[5]*k[8];
/* 2 0*/ x[6]=b[6]*k[0] + b[7]*k[3] + b[8]*k[6];
/* 2 1*/ x[7]=b[6]*k[1] + b[7]*k[4] + b[8]*k[7];
/* 2 2*/ x[8]=b[6]*k[2] + b[7]*k[5] + b[8]*k[8];

/*  sum_{(1 \; 3)\\} ( l )  */
/* ( =  )  */
valuetype z[6];
/* 0 0*/ z[0]=l[0] + l[10] + l[20];
/* 1 0*/ z[1]=l[3] + l[13] + l[23];
/* 2 0*/ z[2]=l[6] + l[16] + l[26];
/* 1 1*/ z[3]=l[12] + l[31] + l[41];
/* 2 1*/ z[4]=l[15] + l[34] + l[44];
/* 2 2*/ z[5]=l[24] + l[43] + l[53];

/*  sum_{(1 \; 2)\\} ( l )  */
/* ( =  )  */
valuetype ab[9];
/* 0 0*/ ab[0]=l[0] + l[12] + l[24];
/* 0 1*/ ab[1]=l[1] + l[13] + l[25];
/* 0 2*/ ab[2]=l[2] + l[14] + l[26];
/* 1 0*/ ab[3]=l[9] + l[30] + l[42];
/* 1 1*/ ab[4]=l[10] + l[31] + l[43];
/* 1 2*/ ab[5]=l[11] + l[32] + l[44];
/* 2 0*/ ab[6]=l[18] + l[39] + l[51];
/* 2 1*/ ab[7]=l[19] + l[40] + l[52];
/* 2 2*/ ab[8]=l[20] + l[41] + l[53];

/*  sum_{(1 \; 2)\\} b m  */
/*  (=sum_{(1 \; 2)\\}()(sum_{(1 \; 3)\\}()()))  */
valuetype bb[9];
/* 0 0*/ bb[0]=b[0]*m[0] + b[1]*m[1] + b[2]*m[2];
/* 0 1*/ bb[1]=b[0]*m[1] + b[1]*m[3] + b[2]*m[4];
/* 0 2*/ bb[2]=b[0]*m[2] + b[1]*m[4] + b[2]*m[5];
/* 1 0*/ bb[3]=b[3]*m[0] + b[4]*m[1] + b[5]*m[2];
/* 1 1*/ bb[4]=b[3]*m[1] + b[4]*m[3] + b[5]*m[4];
/* 1 2*/ bb[5]=b[3]*m[2] + b[4]*m[4] + b[5]*m[5];
/* 2 0*/ bb[6]=b[6]*m[0] + b[7]*m[1] + b[8]*m[2];
/* 2 1*/ bb[7]=b[6]*m[1] + b[7]*m[3] + b[8]*m[4];
/* 2 2*/ bb[8]=b[6]*m[2] + b[7]*m[4] + b[8]*m[5];

/*  sum_{(2 \; 3)\\} c j  */
/*  (=sum_{(2 \; 3)\\}()())  */
valuetype eb[6];
/* 0 0*/ eb[0]=c[0]*j[0] + c[1]*j[1] + c[2]*j[2];
/* 1 0*/ eb[1]=c[3]*j[0] + c[4]*j[1] + c[5]*j[2];
/* 2 0*/ eb[2]=c[6]*j[0] + c[7]*j[1] + c[8]*j[2];
/* 1 1*/ eb[3]=c[9]*j[0] + c[10]*j[1] + c[11]*j[2];
/* 2 1*/ eb[4]=c[12]*j[0] + c[13]*j[1] + c[14]*j[2];
/* 2 2*/ eb[5]=c[15]*j[0] + c[16]*j[1] + c[17]*j[2];

/*  sum_{(0 \; 2)\\} o f  */
/*  (=sum_{(0 \; 2)\\}(sum_{(1 \; 2)\\}(sum_{(0 \; 3)\\}()())())(sum_{(0 \; 3)\\}()()))  */
valuetype fb[3];
/* 0*/ fb[0]=o[0]*f[0] + o[1]*f[1] + o[2]*f[2];
/* 1*/ fb[1]=o[0]*f[1] + o[1]*f[3] + o[2]*f[4];
/* 2*/ fb[2]=o[0]*f[2] + o[1]*f[4] + o[2]*f[5];

/*  sum_{(1 \; 2)\\} v t  */
/*  (=sum_{(1 \; 2)\\}()())  */
valuetype gb[9];
/* 0 0*/ gb[0]=v[0]*t[0] + v[1]*t[1] + v[2]*t[2];
/* 0 1*/ gb[1]=v[0]*t[1] + v[1]*t[3] + v[2]*t[4];
/* 0 2*/ gb[2]=v[0]*t[2] + v[1]*t[4] + v[2]*t[5];
/* 1 0*/ gb[3]=v[1]*t[0] + v[3]*t[1] + v[4]*t[2];
/* 1 1*/ gb[4]=v[1]*t[1] + v[3]*t[3] + v[4]*t[4];
/* 1 2*/ gb[5]=v[1]*t[2] + v[3]*t[4] + v[4]*t[5];
/* 2 0*/ gb[6]=v[2]*t[0] + v[4]*t[1] + v[5]*t[2];
/* 2 1*/ gb[7]=v[2]*t[1] + v[4]*t[3] + v[5]*t[4];
/* 2 2*/ gb[8]=v[2]*t[2] + v[4]*t[4] + v[5]*t[5];

/*  sum_{(1 \; 2)\\} b bb  */
/*  (=sum_{(1 \; 2)\\}()(sum_{(1 \; 2)\\}()(sum_{(1 \; 3)\\}()())))  */
valuetype lb[6];
/* 0 0*/ lb[0]=b[0]*bb[0] + b[1]*bb[1] + b[2]*bb[2];
/* 1 0*/ lb[1]=b[0]*bb[3] + b[1]*bb[4] + b[2]*bb[5];
/* 2 0*/ lb[2]=b[0]*bb[6] + b[1]*bb[7] + b[2]*bb[8];
/* 1 1*/ lb[3]=b[3]*bb[3] + b[4]*bb[4] + b[5]*bb[5];
/* 2 1*/ lb[4]=b[3]*bb[6] + b[4]*bb[7] + b[5]*bb[8];
/* 2 2*/ lb[5]=b[6]*bb[6] + b[7]*bb[7] + b[8]*bb[8];

/*  sum_{(1 \; 3)\\} v v  */
/*  (=sum_{(1 \; 3)\\}()())  */
valuetype mb[6];
/* 0 0*/ mb[0]=v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
/* 1 0*/ mb[1]=v[0]*v[1] + v[1]*v[3] + v[2]*v[4];
/* 2 0*/ mb[2]=v[0]*v[2] + v[1]*v[4] + v[2]*v[5];
/* 1 1*/ mb[3]=v[1]*v[1] + v[3]*v[3] + v[4]*v[4];
/* 2 1*/ mb[4]=v[1]*v[2] + v[3]*v[4] + v[4]*v[5];
/* 2 2*/ mb[5]=v[2]*v[2] + v[4]*v[4] + v[5]*v[5];

/*  sum_{(0 \; 2)\\} d v  */
/*  (=sum_{(0 \; 2)\\}()())  */
valuetype pb[3];
/* 0*/ pb[0]=d[0]*v[0] + d[1]*v[1] + d[2]*v[2];
/* 1*/ pb[1]=d[0]*v[1] + d[1]*v[3] + d[2]*v[4];
/* 2*/ pb[2]=d[0]*v[2] + d[1]*v[4] + d[2]*v[5];

/*  sum_{(0 \; 2)\\} d ab  */
/*  (=sum_{(0 \; 2)\\}()())  */
valuetype qb[3];
/* 0*/ qb[0]=d[0]*ab[0] + d[1]*ab[1] + d[2]*ab[2];
/* 1*/ qb[1]=d[0]*ab[3] + d[1]*ab[4] + d[2]*ab[5];
/* 2*/ qb[2]=d[0]*ab[6] + d[1]*ab[7] + d[2]*ab[8];

/*  sum_{(1 \; 2)\\} v z  */
/*  (=sum_{(1 \; 2)\\}()())  */
valuetype sb[9];
/* 0 0*/ sb[0]=v[0]*z[0] + v[1]*z[1] + v[2]*z[2];
/* 0 1*/ sb[1]=v[0]*z[1] + v[1]*z[3] + v[2]*z[4];
/* 0 2*/ sb[2]=v[0]*z[2] + v[1]*z[4] + v[2]*z[5];
/* 1 0*/ sb[3]=v[1]*z[0] + v[3]*z[1] + v[4]*z[2];
/* 1 1*/ sb[4]=v[1]*z[1] + v[3]*z[3] + v[4]*z[4];
/* 1 2*/ sb[5]=v[1]*z[2] + v[3]*z[4] + v[4]*z[5];
/* 2 0*/ sb[6]=v[2]*z[0] + v[4]*z[1] + v[5]*z[2];
/* 2 1*/ sb[7]=v[2]*z[1] + v[4]*z[3] + v[5]*z[4];
/* 2 2*/ sb[8]=v[2]*z[2] + v[4]*z[4] + v[5]*z[5];

/*  sum_{(1 \; 3)\\} v eb  */
/*  (=sum_{(1 \; 3)\\}()(sum_{(2 \; 3)\\}()()))  */
valuetype tb[9];
/* 0 0*/ tb[0]=v[0]*eb[0] + v[1]*eb[1] + v[2]*eb[2];
/* 0 1*/ tb[1]=v[0]*eb[1] + v[1]*eb[3] + v[2]*eb[4];
/* 0 2*/ tb[2]=v[0]*eb[2] + v[1]*eb[4] + v[2]*eb[5];
/* 1 0*/ tb[3]=v[1]*eb[0] + v[3]*eb[1] + v[4]*eb[2];
/* 1 1*/ tb[4]=v[1]*eb[1] + v[3]*eb[3] + v[4]*eb[4];
/* 1 2*/ tb[5]=v[1]*eb[2] + v[3]*eb[4] + v[4]*eb[5];
/* 2 0*/ tb[6]=v[2]*eb[0] + v[4]*eb[1] + v[5]*eb[2];
/* 2 1*/ tb[7]=v[2]*eb[1] + v[4]*eb[3] + v[5]*eb[4];
/* 2 2*/ tb[8]=v[2]*eb[2] + v[4]*eb[4] + v[5]*eb[5];

/*  sum_{(1 \; 3)\\} v ab  */
/*  (=sum_{(1 \; 3)\\}()())  */
valuetype ub[9];
/* 0 0*/ ub[0]=v[0]*ab[0] + v[1]*ab[1] + v[2]*ab[2];
/* 0 1*/ ub[1]=v[0]*ab[3] + v[1]*ab[4] + v[2]*ab[5];
/* 0 2*/ ub[2]=v[0]*ab[6] + v[1]*ab[7] + v[2]*ab[8];
/* 1 0*/ ub[3]=v[1]*ab[0] + v[3]*ab[1] + v[4]*ab[2];
/* 1 1*/ ub[4]=v[1]*ab[3] + v[3]*ab[4] + v[4]*ab[5];
/* 1 2*/ ub[5]=v[1]*ab[6] + v[3]*ab[7] + v[4]*ab[8];
/* 2 0*/ ub[6]=v[2]*ab[0] + v[4]*ab[1] + v[5]*ab[2];
/* 2 1*/ ub[7]=v[2]*ab[3] + v[4]*ab[4] + v[5]*ab[5];
/* 2 2*/ ub[8]=v[2]*ab[6] + v[4]*ab[7] + v[5]*ab[8];

/*  sum_{(1 \; 2)\\} v q  */
/*  (=sum_{(1 \; 2)\\}()())  */
valuetype vb[9];
/* 0 0*/ vb[0]=v[0]*q[0] + v[1]*q[1] + v[2]*q[2];
/* 0 1*/ vb[1]=v[0]*q[1] + v[1]*q[3] + v[2]*q[4];
/* 0 2*/ vb[2]=v[0]*q[2] + v[1]*q[4] + v[2]*q[5];
/* 1 0*/ vb[3]=v[1]*q[0] + v[3]*q[1] + v[4]*q[2];
/* 1 1*/ vb[4]=v[1]*q[1] + v[3]*q[3] + v[4]*q[4];
/* 1 2*/ vb[5]=v[1]*q[2] + v[3]*q[4] + v[4]*q[5];
/* 2 0*/ vb[6]=v[2]*q[0] + v[4]*q[1] + v[5]*q[2];
/* 2 1*/ vb[7]=v[2]*q[1] + v[4]*q[3] + v[5]*q[4];
/* 2 2*/ vb[8]=v[2]*q[2] + v[4]*q[4] + v[5]*q[5];

/*  sum_{(0 \; 1)\\} a a  */
/*  (=sum_{(0 \; 1)\\}()())  */
/**/ M[0]=a[0]*a[0] + a[1]*a[1] + a[2]*a[2];

/*  sum_{(0 \; 1)\\} a fb  */
/*  (=sum_{(0 \; 1)\\}()(sum_{(0 \; 2)\\}(sum_{(1 \; 2)\\}(sum_{(0 \; 3)\\}()())())(sum_{(0 \; 3)\\}()())))  */
/**/ M[1]=a[0]*fb[0] + a[1]*fb[1] + a[2]*fb[2];

/*  sum_{(0 \; 1)\\} a h  */
/*  (=sum_{(0 \; 1)\\}()(sum_{(0 \; 2)\\}()()))  */
/**/ M[2]=a[0]*h[0] + a[1]*h[1] + a[2]*h[2];

/*  sum_{(0 \; 1)\\} ( b )  */
/* ( =  )  */
/**/ M[3]=b[0] + b[4] + b[8];

/*  sum_{(0 \; 1)\\} ( m )  */
/* ( =  )  */
/**/ M[4]=m[0] + m[3] + m[5];

/*  sum_{(0 \; 1)\\} ( k )  */
/* ( =  )  */
/**/ M[5]=k[0] + k[4] + k[8];

/*  sum_{(0 \; 1)\\} ( bb )  */
/* ( =  )  */
/**/ M[6]=bb[0] + bb[4] + bb[8];

/*  sum_{(0 \; 1)\\} ( x )  */
/* ( =  )  */
/**/ M[7]=x[0] + x[4] + x[8];

/*  sum_{(0 \; 1)\\} ( lb )  */
/* ( =  )  */
/**/ M[8]=lb[0] + lb[3] + lb[5];

/*  sum_{(0 \; 1)\\} u d  */
/*  (=sum_{(0 \; 1)\\}(sum_{(1 \; 2)\\}()())())  */
/**/ M[9]=u[0]*d[0] + u[1]*d[1] + u[2]*d[2];

/*  sum_{(0 \; 1)\\} u j  */
/*  (=sum_{(0 \; 1)\\}(sum_{(1 \; 2)\\}()())())  */
/**/ M[10]=u[0]*j[0] + u[1]*j[1] + u[2]*j[2];

/*  sum_{(0 \; 1)\\} s j  */
/*  (=sum_{(0 \; 1)\\}(sum_{(0 \; 2)\\}()())())  */
/**/ M[11]=s[0]*j[0] + s[1]*j[1] + s[2]*j[2];

/*  sum_{(0 \; 1)\\} d d  */
/*  (=sum_{(0 \; 1)\\}()())  */
/**/ M[12]=d[0]*d[0] + d[1]*d[1] + d[2]*d[2];

/*  sum_{(0 \; 1)\\} d j  */
/*  (=sum_{(0 \; 1)\\}()())  */
/**/ M[13]=d[0]*j[0] + d[1]*j[1] + d[2]*j[2];

/*  sum_{(0 \; 1)\\} j j  */
/*  (=sum_{(0 \; 1)\\}()())  */
/**/ M[14]=j[0]*j[0] + j[1]*j[1] + j[2]*j[2];

/*  sum_{(0 \; 1)\\} ( t )  */
/* ( =  )  */
/**/ M[15]=t[0] + t[3] + t[5];

/*  sum_{(0 \; 1)\\} ( z )  */
/* ( =  )  */
/**/ M[16]=z[0] + z[3] + z[5];

/*  sum_{(0 \; 1)\\} d pb  */
/*  (=sum_{(0 \; 1)\\}()(sum_{(0 \; 2)\\}()()))  */
/**/ M[17]=d[0]*pb[0] + d[1]*pb[1] + d[2]*pb[2];

/*  sum_{(0 \; 1)\\} pb j  */
/*  (=sum_{(0 \; 1)\\}(sum_{(0 \; 2)\\}()())())  */
/**/ M[18]=pb[0]*j[0] + pb[1]*j[1] + pb[2]*j[2];

/*  sum_{(0 \; 1)\\} ( mb )  */
/* ( =  )  */
/**/ M[19]=mb[0] + mb[3] + mb[5];

/*  sum_{(0 \; 1)\\} d qb  */
/*  (=sum_{(0 \; 1)\\}()(sum_{(0 \; 2)\\}()()))  */
/**/ M[20]=d[0]*qb[0] + d[1]*qb[1] + d[2]*qb[2];

/*  sum_{(0 \; 1)\\} qb j  */
/*  (=sum_{(0 \; 1)\\}(sum_{(0 \; 2)\\}()())())  */
/**/ M[21]=qb[0]*j[0] + qb[1]*j[1] + qb[2]*j[2];

/*  sum_{(0 \; 1)\\} ( tb )  */
/* ( =  )  */
/**/ M[22]=tb[0] + tb[4] + tb[8];

/*  sum_{(0 \; 1)\\} ( gb )  */
/* ( =  )  */
/**/ M[23]=gb[0] + gb[4] + gb[8];

/*  sum_{(0 \; 1)\\} ( vb )  */
/* ( =  )  */
/**/ M[24]=vb[0] + vb[4] + vb[8];

/*  sum_{(0 \; 1)\\} ( ub )  */
/* ( =  )  */
/**/ M[25]=ub[0] + ub[4] + ub[8];

/*  sum_{(0 \; 1)\\} ( sb )  */
/* ( =  )  */
/**/ M[26]=sb[0] + sb[4] + sb[8];

}





#undef MYMOMENTSET

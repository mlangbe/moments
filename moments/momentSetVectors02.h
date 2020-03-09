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

//this file is only included in momentSetVectors2.cpp.

//MOMTYPE is a define which has to be set in the including class.


#define MYMOMENTSET momentSet<MOMTYPE,27,0,2>

#ifndef OLD_VC
template<>
#endif
const char*const MYMOMENTSET::comment
=
"computation optimized by using the tensor graph\n "
"and splitting it into multiple tensor products and folds.\n"
"also multiple contractions are performed";

#ifndef OLD_VC
template<>
#endif
void MYMOMENTSET::compute(valuetype*M,
												 const valuetype*A)

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
valuetype ic[3];
/* 0*/ ic[0]=c[0] + c[9] + c[15];
/* 1*/ ic[1]=c[1] + c[10] + c[16];
/* 2*/ ic[2]=c[2] + c[11] + c[17];

/*  sum_{(0 \; 3)\\} a c  */
/*  (=sum_{(0 \; 3)\\}()())  */
valuetype kc[6];
/* 0 0*/ kc[0]=a[0]*c[0] + a[1]*c[1] + a[2]*c[2];
/* 1 0*/ kc[1]=a[0]*c[3] + a[1]*c[4] + a[2]*c[5];
/* 2 0*/ kc[2]=a[0]*c[6] + a[1]*c[7] + a[2]*c[8];
/* 1 1*/ kc[3]=a[0]*c[9] + a[1]*c[10] + a[2]*c[11];
/* 2 1*/ kc[4]=a[0]*c[12] + a[1]*c[13] + a[2]*c[14];
/* 2 2*/ kc[5]=a[0]*c[15] + a[1]*c[16] + a[2]*c[17];

/*  sum_{(0 \; 2)\\} a b  */
/*  (=sum_{(0 \; 2)\\}()())  */
valuetype mc[3];
/* 0*/ mc[0]=a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
/* 1*/ mc[1]=a[0]*b[3] + a[1]*b[4] + a[2]*b[5];
/* 2*/ mc[2]=a[0]*b[6] + a[1]*b[7] + a[2]*b[8];

/*  sum_{(1 \; 4)(0 \; 3)\\} c c  */
/*  (=sum_{(1 \; 4)(0 \; 3)\\}()())  */
valuetype nc[6];
/* 0 0*/ nc[0]=c[0]*c[0] + c[9]*c[9] + c[15]*c[15] + 2 * ( c[3]*c[3] + c[6]*c[6] + c[12]*c[12] );
/* 1 0*/ nc[1]=c[0]*c[1] + c[9]*c[10] + c[15]*c[16] + 2 * ( c[3]*c[4] + c[6]*c[7] + c[12]*c[13] );
/* 2 0*/ nc[2]=c[0]*c[2] + c[9]*c[11] + c[15]*c[17] + 2 * ( c[3]*c[5] + c[6]*c[8] + c[12]*c[14] );
/* 1 1*/ nc[3]=c[1]*c[1] + c[10]*c[10] + c[16]*c[16] + 2 * ( c[4]*c[4] + c[7]*c[7] + c[13]*c[13] );
/* 2 1*/ nc[4]=c[1]*c[2] + c[10]*c[11] + c[16]*c[17] + 2 * ( c[4]*c[5] + c[7]*c[8] + c[13]*c[14] );
/* 2 2*/ nc[5]=c[2]*c[2] + c[11]*c[11] + c[17]*c[17] + 2 * ( c[5]*c[5] + c[8]*c[8] + c[14]*c[14] );

/*  sum_{(1 \; 2)\\} ( c )  */
/* ( =  )  */
valuetype oc[3];
/* 0*/ oc[0]=c[0] + c[4] + c[8];
/* 1*/ oc[1]=c[3] + c[10] + c[14];
/* 2*/ oc[2]=c[6] + c[13] + c[17];

/*  sum_{(1 \; 2)\\} b b  */
/*  (=sum_{(1 \; 2)\\}()())  */
valuetype pc[9];
/* 0 0*/ pc[0]=b[0]*b[0] + b[1]*b[3] + b[2]*b[6];
/* 0 1*/ pc[1]=b[0]*b[1] + b[1]*b[4] + b[2]*b[7];
/* 0 2*/ pc[2]=b[0]*b[2] + b[1]*b[5] + b[2]*b[8];
/* 1 0*/ pc[3]=b[3]*b[0] + b[4]*b[3] + b[5]*b[6];
/* 1 1*/ pc[4]=b[3]*b[1] + b[4]*b[4] + b[5]*b[7];
/* 1 2*/ pc[5]=b[3]*b[2] + b[4]*b[5] + b[5]*b[8];
/* 2 0*/ pc[6]=b[6]*b[0] + b[7]*b[3] + b[8]*b[6];
/* 2 1*/ pc[7]=b[6]*b[1] + b[7]*b[4] + b[8]*b[7];
/* 2 2*/ pc[8]=b[6]*b[2] + b[7]*b[5] + b[8]*b[8];

/*  sum_{(2 \; 4)\\} c c  */
/*  (=sum_{(2 \; 4)\\}()())  */
valuetype qc[54];
/* 0 0 0 0*/ qc[0]=c[0]*c[0] + c[1]*c[3] + c[2]*c[6];
/* 0 0 0 1*/ qc[1]=c[0]*c[1] + c[1]*c[4] + c[2]*c[7];
/* 0 0 0 2*/ qc[2]=c[0]*c[2] + c[1]*c[5] + c[2]*c[8];
/* 0 0 1 0*/ qc[3]=c[0]*c[3] + c[1]*c[9] + c[2]*c[12];
/* 0 0 1 1*/ qc[4]=c[0]*c[4] + c[1]*c[10] + c[2]*c[13];
/* 0 0 1 2*/ qc[5]=c[0]*c[5] + c[1]*c[11] + c[2]*c[14];
/* 0 0 2 0*/ qc[6]=c[0]*c[6] + c[1]*c[12] + c[2]*c[15];
/* 0 0 2 1*/ qc[7]=c[0]*c[7] + c[1]*c[13] + c[2]*c[16];
/* 0 0 2 2*/ qc[8]=c[0]*c[8] + c[1]*c[14] + c[2]*c[17];
/* 1 0 0 0*/ qc[9]=c[3]*c[0] + c[4]*c[3] + c[5]*c[6];
/* 1 0 0 1*/ qc[10]=c[3]*c[1] + c[4]*c[4] + c[5]*c[7];
/* 1 0 0 2*/ qc[11]=c[3]*c[2] + c[4]*c[5] + c[5]*c[8];
/* 1 0 1 0*/ qc[12]=c[3]*c[3] + c[4]*c[9] + c[5]*c[12];
/* 1 0 1 1*/ qc[13]=c[3]*c[4] + c[4]*c[10] + c[5]*c[13];
/* 1 0 1 2*/ qc[14]=c[3]*c[5] + c[4]*c[11] + c[5]*c[14];
/* 1 0 2 0*/ qc[15]=c[3]*c[6] + c[4]*c[12] + c[5]*c[15];
/* 1 0 2 1*/ qc[16]=c[3]*c[7] + c[4]*c[13] + c[5]*c[16];
/* 1 0 2 2*/ qc[17]=c[3]*c[8] + c[4]*c[14] + c[5]*c[17];
/* 2 0 0 0*/ qc[18]=c[6]*c[0] + c[7]*c[3] + c[8]*c[6];
/* 2 0 0 1*/ qc[19]=c[6]*c[1] + c[7]*c[4] + c[8]*c[7];
/* 2 0 0 2*/ qc[20]=c[6]*c[2] + c[7]*c[5] + c[8]*c[8];
/* 2 0 1 0*/ qc[21]=c[6]*c[3] + c[7]*c[9] + c[8]*c[12];
/* 2 0 1 1*/ qc[22]=c[6]*c[4] + c[7]*c[10] + c[8]*c[13];
/* 2 0 1 2*/ qc[23]=c[6]*c[5] + c[7]*c[11] + c[8]*c[14];
/* 2 0 2 0*/ qc[24]=c[6]*c[6] + c[7]*c[12] + c[8]*c[15];
/* 2 0 2 1*/ qc[25]=c[6]*c[7] + c[7]*c[13] + c[8]*c[16];
/* 2 0 2 2*/ qc[26]=c[6]*c[8] + c[7]*c[14] + c[8]*c[17];
/* 1 1 0 0*/ qc[27]=c[9]*c[0] + c[10]*c[3] + c[11]*c[6];
/* 1 1 0 1*/ qc[28]=c[9]*c[1] + c[10]*c[4] + c[11]*c[7];
/* 1 1 0 2*/ qc[29]=c[9]*c[2] + c[10]*c[5] + c[11]*c[8];
/* 1 1 1 0*/ qc[30]=c[9]*c[3] + c[10]*c[9] + c[11]*c[12];
/* 1 1 1 1*/ qc[31]=c[9]*c[4] + c[10]*c[10] + c[11]*c[13];
/* 1 1 1 2*/ qc[32]=c[9]*c[5] + c[10]*c[11] + c[11]*c[14];
/* 1 1 2 0*/ qc[33]=c[9]*c[6] + c[10]*c[12] + c[11]*c[15];
/* 1 1 2 1*/ qc[34]=c[9]*c[7] + c[10]*c[13] + c[11]*c[16];
/* 1 1 2 2*/ qc[35]=c[9]*c[8] + c[10]*c[14] + c[11]*c[17];
/* 2 1 0 0*/ qc[36]=c[12]*c[0] + c[13]*c[3] + c[14]*c[6];
/* 2 1 0 1*/ qc[37]=c[12]*c[1] + c[13]*c[4] + c[14]*c[7];
/* 2 1 0 2*/ qc[38]=c[12]*c[2] + c[13]*c[5] + c[14]*c[8];
/* 2 1 1 0*/ qc[39]=c[12]*c[3] + c[13]*c[9] + c[14]*c[12];
/* 2 1 1 1*/ qc[40]=c[12]*c[4] + c[13]*c[10] + c[14]*c[13];
/* 2 1 1 2*/ qc[41]=c[12]*c[5] + c[13]*c[11] + c[14]*c[14];
/* 2 1 2 0*/ qc[42]=c[12]*c[6] + c[13]*c[12] + c[14]*c[15];
/* 2 1 2 1*/ qc[43]=c[12]*c[7] + c[13]*c[13] + c[14]*c[16];
/* 2 1 2 2*/ qc[44]=c[12]*c[8] + c[13]*c[14] + c[14]*c[17];
/* 2 2 0 0*/ qc[45]=c[15]*c[0] + c[16]*c[3] + c[17]*c[6];
/* 2 2 0 1*/ qc[46]=c[15]*c[1] + c[16]*c[4] + c[17]*c[7];
/* 2 2 0 2*/ qc[47]=c[15]*c[2] + c[16]*c[5] + c[17]*c[8];
/* 2 2 1 0*/ qc[48]=c[15]*c[3] + c[16]*c[9] + c[17]*c[12];
/* 2 2 1 1*/ qc[49]=c[15]*c[4] + c[16]*c[10] + c[17]*c[13];
/* 2 2 1 2*/ qc[50]=c[15]*c[5] + c[16]*c[11] + c[17]*c[14];
/* 2 2 2 0*/ qc[51]=c[15]*c[6] + c[16]*c[12] + c[17]*c[15];
/* 2 2 2 1*/ qc[52]=c[15]*c[7] + c[16]*c[13] + c[17]*c[16];
/* 2 2 2 2*/ qc[53]=c[15]*c[8] + c[16]*c[14] + c[17]*c[17];

/*  sum_{(1 \; 3)\\} b b  */
/*  (=sum_{(1 \; 3)\\}()())  */
valuetype rc[6];
/* 0 0*/ rc[0]=b[0]*b[0] + b[1]*b[1] + b[2]*b[2];
/* 1 0*/ rc[1]=b[0]*b[3] + b[1]*b[4] + b[2]*b[5];
/* 2 0*/ rc[2]=b[0]*b[6] + b[1]*b[7] + b[2]*b[8];
/* 1 1*/ rc[3]=b[3]*b[3] + b[4]*b[4] + b[5]*b[5];
/* 2 1*/ rc[4]=b[3]*b[6] + b[4]*b[7] + b[5]*b[8];
/* 2 2*/ rc[5]=b[6]*b[6] + b[7]*b[7] + b[8]*b[8];

/*  sum_{(2 \; 5)\\} c c  */
/*  (=sum_{(2 \; 5)\\}()())  */
valuetype sc[21];
/* 0 0 0 0*/ sc[0]=c[0]*c[0] + c[1]*c[1] + c[2]*c[2];
/* 1 0 0 0*/ sc[1]=c[0]*c[3] + c[1]*c[4] + c[2]*c[5];
/* 2 0 0 0*/ sc[2]=c[0]*c[6] + c[1]*c[7] + c[2]*c[8];
/* 1 1 0 0*/ sc[3]=c[0]*c[9] + c[1]*c[10] + c[2]*c[11];
/* 2 1 0 0*/ sc[4]=c[0]*c[12] + c[1]*c[13] + c[2]*c[14];
/* 2 2 0 0*/ sc[5]=c[0]*c[15] + c[1]*c[16] + c[2]*c[17];
/* 1 0 1 0*/ sc[6]=c[3]*c[3] + c[4]*c[4] + c[5]*c[5];
/* 2 0 1 0*/ sc[7]=c[3]*c[6] + c[4]*c[7] + c[5]*c[8];
/* 1 1 1 0*/ sc[8]=c[3]*c[9] + c[4]*c[10] + c[5]*c[11];
/* 2 1 1 0*/ sc[9]=c[3]*c[12] + c[4]*c[13] + c[5]*c[14];
/* 2 2 1 0*/ sc[10]=c[3]*c[15] + c[4]*c[16] + c[5]*c[17];
/* 2 0 2 0*/ sc[11]=c[6]*c[6] + c[7]*c[7] + c[8]*c[8];
/* 2 0 1 1*/ sc[12]=c[6]*c[9] + c[7]*c[10] + c[8]*c[11];
/* 2 1 2 0*/ sc[13]=c[6]*c[12] + c[7]*c[13] + c[8]*c[14];
/* 2 2 2 0*/ sc[14]=c[6]*c[15] + c[7]*c[16] + c[8]*c[17];
/* 1 1 1 1*/ sc[15]=c[9]*c[9] + c[10]*c[10] + c[11]*c[11];
/* 2 1 1 1*/ sc[16]=c[9]*c[12] + c[10]*c[13] + c[11]*c[14];
/* 2 2 1 1*/ sc[17]=c[9]*c[15] + c[10]*c[16] + c[11]*c[17];
/* 2 1 2 1*/ sc[18]=c[12]*c[12] + c[13]*c[13] + c[14]*c[14];
/* 2 2 2 1*/ sc[19]=c[12]*c[15] + c[13]*c[16] + c[14]*c[17];
/* 2 2 2 2*/ sc[20]=c[15]*c[15] + c[16]*c[16] + c[17]*c[17];

/*  sum_{(1 \; 2)\\} kc ic  */
/*  (=sum_{(1 \; 2)\\}(sum_{(0 \; 3)\\}()())())  */
valuetype tc[3];
/* 0*/ tc[0]=kc[0]*ic[0] + kc[1]*ic[1] + kc[2]*ic[2];
/* 1*/ tc[1]=kc[1]*ic[0] + kc[3]*ic[1] + kc[4]*ic[2];
/* 2*/ tc[2]=kc[2]*ic[0] + kc[4]*ic[1] + kc[5]*ic[2];

/*  sum_{(0 \; 2)\\} b ic  */
/*  (=sum_{(0 \; 2)\\}()())  */
valuetype wc[3];
/* 0*/ wc[0]=b[0]*ic[0] + b[3]*ic[1] + b[6]*ic[2];
/* 1*/ wc[1]=b[1]*ic[0] + b[4]*ic[1] + b[7]*ic[2];
/* 2*/ wc[2]=b[2]*ic[0] + b[5]*ic[1] + b[8]*ic[2];

/*  sum_{(1 \; 3)\\} ( sc )  */
/* ( =  )  */
valuetype xc[6];
/* 0 0*/ xc[0]=sc[0] + sc[6] + sc[11];
/* 1 0*/ xc[1]=sc[1] + sc[8] + sc[13];
/* 2 0*/ xc[2]=sc[2] + sc[9] + sc[14];
/* 1 1*/ xc[3]=sc[6] + sc[15] + sc[18];
/* 2 1*/ xc[4]=sc[7] + sc[16] + sc[19];
/* 2 2*/ xc[5]=sc[11] + sc[18] + sc[20];

/*  sum_{(1 \; 2)\\} b ic  */
/*  (=sum_{(1 \; 2)\\}()())  */
valuetype yc[3];
/* 0*/ yc[0]=b[0]*ic[0] + b[1]*ic[1] + b[2]*ic[2];
/* 1*/ yc[1]=b[3]*ic[0] + b[4]*ic[1] + b[5]*ic[2];
/* 2*/ yc[2]=b[6]*ic[0] + b[7]*ic[1] + b[8]*ic[2];

/*  sum_{(1 \; 0)\\} ( sc )  */
/* ( =  )  */
valuetype zc[6];
/* 0 0*/ zc[0]=sc[0] + sc[3] + sc[5];
/* 1 0*/ zc[1]=sc[1] + sc[8] + sc[10];
/* 2 0*/ zc[2]=sc[2] + sc[12] + sc[14];
/* 1 1*/ zc[3]=sc[3] + sc[15] + sc[17];
/* 2 1*/ zc[4]=sc[4] + sc[16] + sc[19];
/* 2 2*/ zc[5]=sc[5] + sc[17] + sc[20];

/*  sum_{(1 \; 3)\\} ( qc )  */
/* ( =  )  */
valuetype dd[6];
/* 0 0*/ dd[0]=qc[0] + qc[10] + qc[20];
/* 1 0*/ dd[1]=qc[3] + qc[13] + qc[23];
/* 2 0*/ dd[2]=qc[6] + qc[16] + qc[26];
/* 1 1*/ dd[3]=qc[12] + qc[31] + qc[41];
/* 2 1*/ dd[4]=qc[15] + qc[34] + qc[44];
/* 2 2*/ dd[5]=qc[24] + qc[43] + qc[53];

/*  sum_{(1 \; 2)\\} ( qc )  */
/* ( =  )  */
valuetype ed[9];
/* 0 0*/ ed[0]=qc[0] + qc[12] + qc[24];
/* 0 1*/ ed[1]=qc[1] + qc[13] + qc[25];
/* 0 2*/ ed[2]=qc[2] + qc[14] + qc[26];
/* 1 0*/ ed[3]=qc[9] + qc[30] + qc[42];
/* 1 1*/ ed[4]=qc[10] + qc[31] + qc[43];
/* 1 2*/ ed[5]=qc[11] + qc[32] + qc[44];
/* 2 0*/ ed[6]=qc[18] + qc[39] + qc[51];
/* 2 1*/ ed[7]=qc[19] + qc[40] + qc[52];
/* 2 2*/ ed[8]=qc[20] + qc[41] + qc[53];

/*  sum_{(1 \; 2)\\} b rc  */
/*  (=sum_{(1 \; 2)\\}()(sum_{(1 \; 3)\\}()()))  */
valuetype fd[9];
/* 0 0*/ fd[0]=b[0]*rc[0] + b[1]*rc[1] + b[2]*rc[2];
/* 0 1*/ fd[1]=b[0]*rc[1] + b[1]*rc[3] + b[2]*rc[4];
/* 0 2*/ fd[2]=b[0]*rc[2] + b[1]*rc[4] + b[2]*rc[5];
/* 1 0*/ fd[3]=b[3]*rc[0] + b[4]*rc[1] + b[5]*rc[2];
/* 1 1*/ fd[4]=b[3]*rc[1] + b[4]*rc[3] + b[5]*rc[4];
/* 1 2*/ fd[5]=b[3]*rc[2] + b[4]*rc[4] + b[5]*rc[5];
/* 2 0*/ fd[6]=b[6]*rc[0] + b[7]*rc[1] + b[8]*rc[2];
/* 2 1*/ fd[7]=b[6]*rc[1] + b[7]*rc[3] + b[8]*rc[4];
/* 2 2*/ fd[8]=b[6]*rc[2] + b[7]*rc[4] + b[8]*rc[5];

/*  sum_{(2 \; 3)\\} c oc  */
/*  (=sum_{(2 \; 3)\\}()())  */
valuetype id[6];
/* 0 0*/ id[0]=c[0]*oc[0] + c[1]*oc[1] + c[2]*oc[2];
/* 1 0*/ id[1]=c[3]*oc[0] + c[4]*oc[1] + c[5]*oc[2];
/* 2 0*/ id[2]=c[6]*oc[0] + c[7]*oc[1] + c[8]*oc[2];
/* 1 1*/ id[3]=c[9]*oc[0] + c[10]*oc[1] + c[11]*oc[2];
/* 2 1*/ id[4]=c[12]*oc[0] + c[13]*oc[1] + c[14]*oc[2];
/* 2 2*/ id[5]=c[15]*oc[0] + c[16]*oc[1] + c[17]*oc[2];

/*  sum_{(0 \; 2)\\} tc kc  */
/*  (=sum_{(0 \; 2)\\}(sum_{(1 \; 2)\\}(sum_{(0 \; 3)\\}()())())(sum_{(0 \; 3)\\}()()))  */
valuetype jd[3];
/* 0*/ jd[0]=tc[0]*kc[0] + tc[1]*kc[1] + tc[2]*kc[2];
/* 1*/ jd[1]=tc[0]*kc[1] + tc[1]*kc[3] + tc[2]*kc[4];
/* 2*/ jd[2]=tc[0]*kc[2] + tc[1]*kc[4] + tc[2]*kc[5];

/*  sum_{(0 \; 2)\\} ic zc  */
/*  (=sum_{(0 \; 2)\\}()())  */
valuetype sd[3];
/* 0*/ sd[0]=ic[0]*zc[0] + ic[1]*zc[1] + ic[2]*zc[2];
/* 1*/ sd[1]=ic[0]*zc[1] + ic[1]*zc[3] + ic[2]*zc[4];
/* 2*/ sd[2]=ic[0]*zc[2] + ic[1]*zc[4] + ic[2]*zc[5];

/*  sum_{(0 \; 2)\\} ic ed  */
/*  (=sum_{(0 \; 2)\\}()())  */
valuetype td[3];
/* 0*/ td[0]=ic[0]*ed[0] + ic[1]*ed[1] + ic[2]*ed[2];
/* 1*/ td[1]=ic[0]*ed[3] + ic[1]*ed[4] + ic[2]*ed[5];
/* 2*/ td[2]=ic[0]*ed[6] + ic[1]*ed[7] + ic[2]*ed[8];

/*  sum_{(0 \; 1)\\} a a  */
/*  (=sum_{(0 \; 1)\\}()())  */
/**/ M[0]=a[0]*a[0] + a[1]*a[1] + a[2]*a[2];

/*  sum_{(0 \; 1)\\} a jd  */
/*  (=sum_{(0 \; 1)\\}()(sum_{(0 \; 2)\\}(sum_{(1 \; 2)\\}(sum_{(0 \; 3)\\}()())())(sum_{(0 \; 3)\\}()())))  */
/**/ M[1]=a[0]*jd[0] + a[1]*jd[1] + a[2]*jd[2];

/*  sum_{(0 \; 1)\\} a mc  */
/*  (=sum_{(0 \; 1)\\}()(sum_{(0 \; 2)\\}()()))  */
/**/ M[2]=a[0]*mc[0] + a[1]*mc[1] + a[2]*mc[2];

/*  sum_{(0 \; 1)\\} ( b )  */
/* ( =  )  */
/**/ M[3]=b[0] + b[4] + b[8];

/*  sum_{(0 \; 1)\\} ( rc )  */
/* ( =  )  */
/**/ M[4]=rc[0] + rc[3] + rc[5];

/*  sum_{(0 \; 1)\\} ( pc )  */
/* ( =  )  */
/**/ M[5]=pc[0] + pc[4] + pc[8];

/*  sum_{(0 \; 1)\\} ( fd )  */
/* ( =  )  */
/**/ M[6]=fd[0] + fd[4] + fd[8];

/*  sum_{(0 \; 3)(1 \; 2)\\} b pc  */
/*  (=sum_{(0 \; 3)(1 \; 2)\\}()(sum_{(1 \; 2)\\}()()))  */
/**/ M[7]=b[0]*pc[0] + b[1]*pc[3] + b[2]*pc[6] + b[3]*pc[1] + b[4]*pc[4] + b[5]*pc[7] + b[6]*pc[2] + b[7]*pc[5] + b[8]*pc[8];

/*  sum_{(0 \; 3)(1 \; 2)\\} b fd  */
/*  (=sum_{(0 \; 3)(1 \; 2)\\}()(sum_{(1 \; 2)\\}()(sum_{(1 \; 3)\\}()())))  */
/**/ M[8]=b[0]*fd[0] + b[1]*fd[1] + b[2]*fd[2] + b[3]*fd[3] + b[4]*fd[4] + b[5]*fd[5] + b[6]*fd[6] + b[7]*fd[7] + b[8]*fd[8];

/*  sum_{(0 \; 1)\\} yc ic  */
/*  (=sum_{(0 \; 1)\\}(sum_{(1 \; 2)\\}()())())  */
/**/ M[9]=yc[0]*ic[0] + yc[1]*ic[1] + yc[2]*ic[2];

/*  sum_{(0 \; 1)\\} yc oc  */
/*  (=sum_{(0 \; 1)\\}(sum_{(1 \; 2)\\}()())())  */
/**/ M[10]=yc[0]*oc[0] + yc[1]*oc[1] + yc[2]*oc[2];

/*  sum_{(0 \; 1)\\} wc oc  */
/*  (=sum_{(0 \; 1)\\}(sum_{(0 \; 2)\\}()())())  */
/**/ M[11]=wc[0]*oc[0] + wc[1]*oc[1] + wc[2]*oc[2];

/*  sum_{(0 \; 1)\\} ic ic  */
/*  (=sum_{(0 \; 1)\\}()())  */
/**/ M[12]=ic[0]*ic[0] + ic[1]*ic[1] + ic[2]*ic[2];

/*  sum_{(0 \; 1)\\} ic oc  */
/*  (=sum_{(0 \; 1)\\}()())  */
/**/ M[13]=ic[0]*oc[0] + ic[1]*oc[1] + ic[2]*oc[2];

/*  sum_{(0 \; 1)\\} oc oc  */
/*  (=sum_{(0 \; 1)\\}()())  */
/**/ M[14]=oc[0]*oc[0] + oc[1]*oc[1] + oc[2]*oc[2];

/*  sum_{(0 \; 1)\\} ( xc )  */
/* ( =  )  */
/**/ M[15]=xc[0] + xc[3] + xc[5];

/*  sum_{(0 \; 1)\\} ( dd )  */
/* ( =  )  */
/**/ M[16]=dd[0] + dd[3] + dd[5];

/*  sum_{(0 \; 1)\\} ic sd  */
/*  (=sum_{(0 \; 1)\\}()(sum_{(0 \; 2)\\}()()))  */
/**/ M[17]=ic[0]*sd[0] + ic[1]*sd[1] + ic[2]*sd[2];

/*  sum_{(0 \; 1)\\} sd oc  */
/*  (=sum_{(0 \; 1)\\}(sum_{(0 \; 2)\\}()())())  */
/**/ M[18]=sd[0]*oc[0] + sd[1]*oc[1] + sd[2]*oc[2];

/*  sum_{(1 \; 3)(0 \; 2)\\} zc zc  */
/*  (=sum_{(1 \; 3)(0 \; 2)\\}()())  */
/**/ M[19]=zc[0]*zc[0] + zc[3]*zc[3] + zc[5]*zc[5] + 2 * ( zc[1]*zc[1] + zc[2]*zc[2] + zc[4]*zc[4] );

/*  sum_{(0 \; 1)\\} ic td  */
/*  (=sum_{(0 \; 1)\\}()(sum_{(0 \; 2)\\}()()))  */
/**/ M[20]=ic[0]*td[0] + ic[1]*td[1] + ic[2]*td[2];

/*  sum_{(0 \; 1)\\} td oc  */
/*  (=sum_{(0 \; 1)\\}(sum_{(0 \; 2)\\}()())())  */
/**/ M[21]=td[0]*oc[0] + td[1]*oc[1] + td[2]*oc[2];

/*  sum_{(1 \; 3)(0 \; 2)\\} zc id  */
/*  (=sum_{(1 \; 3)(0 \; 2)\\}()(sum_{(2 \; 3)\\}()()))  */
/**/ M[22]=zc[0]*id[0] + zc[3]*id[3] + zc[5]*id[5] + 2 * ( zc[1]*id[1] + zc[2]*id[2] + zc[4]*id[4] );

/*  sum_{(1 \; 3)(0 \; 2)\\} zc xc  */
/*  (=sum_{(1 \; 3)(0 \; 2)\\}()())  */
/**/ M[23]=zc[0]*xc[0] + zc[3]*xc[3] + zc[5]*xc[5] + 2 * ( zc[1]*xc[1] + zc[2]*xc[2] + zc[4]*xc[4] );

/*  sum_{(1 \; 3)(0 \; 2)\\} zc nc  */
/*  (=sum_{(1 \; 3)(0 \; 2)\\}()(sum_{(1 \; 4)(0 \; 3)\\}()()))  */
/**/ M[24]=zc[0]*nc[0] + zc[3]*nc[3] + zc[5]*nc[5] + 2 * ( zc[1]*nc[1] + zc[2]*nc[2] + zc[4]*nc[4] );

/*  sum_{(1 \; 2)(0 \; 3)\\} zc ed  */
/*  (=sum_{(1 \; 2)(0 \; 3)\\}()())  */
/**/ M[25]=zc[0]*ed[0] + zc[1]*ed[1] + zc[1]*ed[3] + zc[2]*ed[2] + zc[2]*ed[6] + zc[3]*ed[4] + zc[4]*ed[5] + zc[4]*ed[7] + zc[5]*ed[8];

/*  sum_{(1 \; 3)(0 \; 2)\\} zc dd  */
/*  (=sum_{(1 \; 3)(0 \; 2)\\}()())  */
/**/ M[26]=zc[0]*dd[0] + zc[3]*dd[3] + zc[5]*dd[5] + 2 * ( zc[1]*dd[1] + zc[2]*dd[2] + zc[4]*dd[4] );

}

#undef MYMOMENTSET

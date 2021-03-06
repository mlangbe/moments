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

#define MYMOMENTSET momentSet<MOMTYPE,28,0,2>

#ifndef OLD_VC
template<>
#endif
const char*const MYMOMENTSET::comment
="computation optimized by using an optimized tensor Graph";

#ifndef OLD_VC
template<>
#endif
void MYMOMENTSET::compute(valuetype*M,const valuetype*A){

//
valuetype a[3];
/* 0*/ a[0]= ( A[15] ) ;
/* 1*/ a[1]= ( A[5] ) ;
/* 2*/ a[2]= ( A[1] ) ;

//
valuetype b[6];
/* 0 0*/ b[0]= ( A[25] ) ;
/* 1 0*/ b[1]= ( A[19] ) ;
/* 2 0*/ b[2]= ( A[16] ) ;
/* 1 1*/ b[3]= ( A[9] ) ;
/* 2 1*/ b[4]= ( A[6] ) ;
/* 2 2*/ b[5]= ( A[2] ) ;

//
valuetype c[10];
/* 0 0 0*/ c[0]= ( A[31] ) ;
/* 1 0 0*/ c[1]= ( A[28] ) ;
/* 2 0 0*/ c[2]= ( A[26] ) ;
/* 1 1 0*/ c[3]= ( A[22] ) ;
/* 2 1 0*/ c[4]= ( A[20] ) ;
/* 2 2 0*/ c[5]= ( A[17] ) ;
/* 1 1 1*/ c[6]= ( A[12] ) ;
/* 2 1 1*/ c[7]= ( A[10] ) ;
/* 2 2 1*/ c[8]= ( A[7] ) ;
/* 2 2 2*/ c[9]= ( A[3] ) ;

//
valuetype d[15];
/* 0 0 0 0*/ d[0]= ( A[34] ) ;
/* 1 0 0 0*/ d[1]= ( A[33] ) ;
/* 2 0 0 0*/ d[2]= ( A[32] ) ;
/* 1 1 0 0*/ d[3]= ( A[30] ) ;
/* 2 1 0 0*/ d[4]= ( A[29] ) ;
/* 2 2 0 0*/ d[5]= ( A[27] ) ;
/* 1 1 1 0*/ d[6]= ( A[24] ) ;
/* 2 1 1 0*/ d[7]= ( A[23] ) ;
/* 2 2 1 0*/ d[8]= ( A[21] ) ;
/* 2 2 2 0*/ d[9]= ( A[18] ) ;
/* 1 1 1 1*/ d[10]= ( A[14] ) ;
/* 2 1 1 1*/ d[11]= ( A[13] ) ;
/* 2 2 1 1*/ d[12]= ( A[11] ) ;
/* 2 2 2 1*/ d[13]= ( A[8] ) ;
/* 2 2 2 2*/ d[14]= ( A[4] ) ;

/*  sum_{(1 \; 3)\\} b b  */
/*  (=sum_{(1 \; 3)\\}()())  */
valuetype zc[6];
/* 0 0*/ zc[0]=b[0]*b[0] + b[1]*b[1] + b[2]*b[2];
/* 1 0*/ zc[1]=b[0]*b[1] + b[1]*b[3] + b[2]*b[4];
/* 2 0*/ zc[2]=b[0]*b[2] + b[1]*b[4] + b[2]*b[5];
/* 1 1*/ zc[3]=b[1]*b[1] + b[3]*b[3] + b[4]*b[4];
/* 2 1*/ zc[4]=b[1]*b[2] + b[3]*b[4] + b[4]*b[5];
/* 2 2*/ zc[5]=b[2]*b[2] + b[4]*b[4] + b[5]*b[5];

/*  sum_{(1 \; 5)\\} b d  */
/*  (=sum_{(1 \; 5)\\}()())  */
valuetype ad[30];
/* 0 0 0 0*/ ad[0]=b[0]*d[0] + b[1]*d[1] + b[2]*d[2];
/* 0 1 0 0*/ ad[1]=b[0]*d[1] + b[1]*d[3] + b[2]*d[4];
/* 0 2 0 0*/ ad[2]=b[0]*d[2] + b[1]*d[4] + b[2]*d[5];
/* 0 1 1 0*/ ad[3]=b[0]*d[3] + b[1]*d[6] + b[2]*d[7];
/* 0 2 1 0*/ ad[4]=b[0]*d[4] + b[1]*d[7] + b[2]*d[8];
/* 0 2 2 0*/ ad[5]=b[0]*d[5] + b[1]*d[8] + b[2]*d[9];
/* 0 1 1 1*/ ad[6]=b[0]*d[6] + b[1]*d[10] + b[2]*d[11];
/* 0 2 1 1*/ ad[7]=b[0]*d[7] + b[1]*d[11] + b[2]*d[12];
/* 0 2 2 1*/ ad[8]=b[0]*d[8] + b[1]*d[12] + b[2]*d[13];
/* 0 2 2 2*/ ad[9]=b[0]*d[9] + b[1]*d[13] + b[2]*d[14];
/* 1 0 0 0*/ ad[10]=b[1]*d[0] + b[3]*d[1] + b[4]*d[2];
/* 1 1 0 0*/ ad[11]=b[1]*d[1] + b[3]*d[3] + b[4]*d[4];
/* 1 2 0 0*/ ad[12]=b[1]*d[2] + b[3]*d[4] + b[4]*d[5];
/* 1 1 1 0*/ ad[13]=b[1]*d[3] + b[3]*d[6] + b[4]*d[7];
/* 1 2 1 0*/ ad[14]=b[1]*d[4] + b[3]*d[7] + b[4]*d[8];
/* 1 2 2 0*/ ad[15]=b[1]*d[5] + b[3]*d[8] + b[4]*d[9];
/* 1 1 1 1*/ ad[16]=b[1]*d[6] + b[3]*d[10] + b[4]*d[11];
/* 1 2 1 1*/ ad[17]=b[1]*d[7] + b[3]*d[11] + b[4]*d[12];
/* 1 2 2 1*/ ad[18]=b[1]*d[8] + b[3]*d[12] + b[4]*d[13];
/* 1 2 2 2*/ ad[19]=b[1]*d[9] + b[3]*d[13] + b[4]*d[14];
/* 2 0 0 0*/ ad[20]=b[2]*d[0] + b[4]*d[1] + b[5]*d[2];
/* 2 1 0 0*/ ad[21]=b[2]*d[1] + b[4]*d[3] + b[5]*d[4];
/* 2 2 0 0*/ ad[22]=b[2]*d[2] + b[4]*d[4] + b[5]*d[5];
/* 2 1 1 0*/ ad[23]=b[2]*d[3] + b[4]*d[6] + b[5]*d[7];
/* 2 2 1 0*/ ad[24]=b[2]*d[4] + b[4]*d[7] + b[5]*d[8];
/* 2 2 2 0*/ ad[25]=b[2]*d[5] + b[4]*d[8] + b[5]*d[9];
/* 2 1 1 1*/ ad[26]=b[2]*d[6] + b[4]*d[10] + b[5]*d[11];
/* 2 2 1 1*/ ad[27]=b[2]*d[7] + b[4]*d[11] + b[5]*d[12];
/* 2 2 2 1*/ ad[28]=b[2]*d[8] + b[4]*d[12] + b[5]*d[13];
/* 2 2 2 2*/ ad[29]=b[2]*d[9] + b[4]*d[13] + b[5]*d[14];

/*  sum_{(2 \; 1)\\} ( c )  */
/* ( =  )  */
valuetype bd[3];
/* 0*/ bd[0]=c[0] + c[3] + c[5];
/* 1*/ bd[1]=c[1] + c[6] + c[8];
/* 2*/ bd[2]=c[2] + c[7] + c[9];

/*  sum_{(2 \; 5)\\} c c  */
/*  (=sum_{(2 \; 5)\\}()())  */
valuetype cd[21];
/* 0 0 0 0*/ cd[0]=c[0]*c[0] + c[1]*c[1] + c[2]*c[2];
/* 1 0 0 0*/ cd[1]=c[0]*c[1] + c[1]*c[3] + c[2]*c[4];
/* 2 0 0 0*/ cd[2]=c[0]*c[2] + c[1]*c[4] + c[2]*c[5];
/* 1 1 0 0*/ cd[3]=c[0]*c[3] + c[1]*c[6] + c[2]*c[7];
/* 2 1 0 0*/ cd[4]=c[0]*c[4] + c[1]*c[7] + c[2]*c[8];
/* 2 2 0 0*/ cd[5]=c[0]*c[5] + c[1]*c[8] + c[2]*c[9];
/* 1 0 1 0*/ cd[6]=c[1]*c[1] + c[3]*c[3] + c[4]*c[4];
/* 2 0 1 0*/ cd[7]=c[1]*c[2] + c[3]*c[4] + c[4]*c[5];
/* 1 1 1 0*/ cd[8]=c[1]*c[3] + c[3]*c[6] + c[4]*c[7];
/* 2 1 1 0*/ cd[9]=c[1]*c[4] + c[3]*c[7] + c[4]*c[8];
/* 2 2 1 0*/ cd[10]=c[1]*c[5] + c[3]*c[8] + c[4]*c[9];
/* 2 0 2 0*/ cd[11]=c[2]*c[2] + c[4]*c[4] + c[5]*c[5];
/* 2 0 1 1*/ cd[12]=c[2]*c[3] + c[4]*c[6] + c[5]*c[7];
/* 2 1 2 0*/ cd[13]=c[2]*c[4] + c[4]*c[7] + c[5]*c[8];
/* 2 2 2 0*/ cd[14]=c[2]*c[5] + c[4]*c[8] + c[5]*c[9];
/* 1 1 1 1*/ cd[15]=c[3]*c[3] + c[6]*c[6] + c[7]*c[7];
/* 2 1 1 1*/ cd[16]=c[3]*c[4] + c[6]*c[7] + c[7]*c[8];
/* 2 2 1 1*/ cd[17]=c[3]*c[5] + c[6]*c[8] + c[7]*c[9];
/* 2 1 2 1*/ cd[18]=c[4]*c[4] + c[7]*c[7] + c[8]*c[8];
/* 2 2 2 1*/ cd[19]=c[4]*c[5] + c[7]*c[8] + c[8]*c[9];
/* 2 2 2 2*/ cd[20]=c[5]*c[5] + c[8]*c[8] + c[9]*c[9];

/*  sum_{(3 \; 2)\\} ( d )  */
/* ( =  )  */
valuetype dd[6];
/* 0 0*/ dd[0]=d[0] + d[3] + d[5];
/* 1 0*/ dd[1]=d[1] + d[6] + d[8];
/* 2 0*/ dd[2]=d[2] + d[7] + d[9];
/* 1 1*/ dd[3]=d[3] + d[10] + d[12];
/* 2 1*/ dd[4]=d[4] + d[11] + d[13];
/* 2 2*/ dd[5]=d[5] + d[12] + d[14];

/*  sum_{(3 \; 7)(2 \; 6)\\} d d  */
/*  (=sum_{(3 \; 7)(2 \; 6)\\}()())  */
valuetype ed[21];
/* 0 0 0 0*/ ed[0]=d[0]*d[0] + d[3]*d[3] + d[5]*d[5] + 2 * ( d[1]*d[1] + d[2]*d[2] + d[4]*d[4] );
/* 1 0 0 0*/ ed[1]=d[0]*d[1] + d[3]*d[6] + d[5]*d[8] + 2 * ( d[1]*d[3] + d[2]*d[4] + d[4]*d[7] );
/* 2 0 0 0*/ ed[2]=d[0]*d[2] + d[3]*d[7] + d[5]*d[9] + 2 * ( d[1]*d[4] + d[2]*d[5] + d[4]*d[8] );
/* 1 1 0 0*/ ed[3]=d[0]*d[3] + d[3]*d[10] + d[5]*d[12] + 2 * ( d[1]*d[6] + d[2]*d[7] + d[4]*d[11] );
/* 2 1 0 0*/ ed[4]=d[0]*d[4] + d[3]*d[11] + d[5]*d[13] + 2 * ( d[1]*d[7] + d[2]*d[8] + d[4]*d[12] );
/* 2 2 0 0*/ ed[5]=d[0]*d[5] + d[3]*d[12] + d[5]*d[14] + 2 * ( d[1]*d[8] + d[2]*d[9] + d[4]*d[13] );
/* 1 0 1 0*/ ed[6]=d[1]*d[1] + d[6]*d[6] + d[8]*d[8] + 2 * ( d[3]*d[3] + d[4]*d[4] + d[7]*d[7] );
/* 2 0 1 0*/ ed[7]=d[1]*d[2] + d[6]*d[7] + d[8]*d[9] + 2 * ( d[3]*d[4] + d[4]*d[5] + d[7]*d[8] );
/* 1 1 1 0*/ ed[8]=d[1]*d[3] + d[6]*d[10] + d[8]*d[12] + 2 * ( d[3]*d[6] + d[4]*d[7] + d[7]*d[11] );
/* 2 1 1 0*/ ed[9]=d[1]*d[4] + d[6]*d[11] + d[8]*d[13] + 2 * ( d[3]*d[7] + d[4]*d[8] + d[7]*d[12] );
/* 2 2 1 0*/ ed[10]=d[1]*d[5] + d[6]*d[12] + d[8]*d[14] + 2 * ( d[3]*d[8] + d[4]*d[9] + d[7]*d[13] );
/* 2 0 2 0*/ ed[11]=d[2]*d[2] + d[7]*d[7] + d[9]*d[9] + 2 * ( d[4]*d[4] + d[5]*d[5] + d[8]*d[8] );
/* 2 0 1 1*/ ed[12]=d[2]*d[3] + d[7]*d[10] + d[9]*d[12] + 2 * ( d[4]*d[6] + d[5]*d[7] + d[8]*d[11] );
/* 2 1 2 0*/ ed[13]=d[2]*d[4] + d[7]*d[11] + d[9]*d[13] + 2 * ( d[4]*d[7] + d[5]*d[8] + d[8]*d[12] );
/* 2 2 2 0*/ ed[14]=d[2]*d[5] + d[7]*d[12] + d[9]*d[14] + 2 * ( d[4]*d[8] + d[5]*d[9] + d[8]*d[13] );
/* 1 1 1 1*/ ed[15]=d[3]*d[3] + d[10]*d[10] + d[12]*d[12] + 2 * ( d[6]*d[6] + d[7]*d[7] + d[11]*d[11] );
/* 2 1 1 1*/ ed[16]=d[3]*d[4] + d[10]*d[11] + d[12]*d[13] + 2 * ( d[6]*d[7] + d[7]*d[8] + d[11]*d[12] );
/* 2 2 1 1*/ ed[17]=d[3]*d[5] + d[10]*d[12] + d[12]*d[14] + 2 * ( d[6]*d[8] + d[7]*d[9] + d[11]*d[13] );
/* 2 1 2 1*/ ed[18]=d[4]*d[4] + d[11]*d[11] + d[13]*d[13] + 2 * ( d[7]*d[7] + d[8]*d[8] + d[12]*d[12] );
/* 2 2 2 1*/ ed[19]=d[4]*d[5] + d[11]*d[12] + d[13]*d[14] + 2 * ( d[7]*d[8] + d[8]*d[9] + d[12]*d[13] );
/* 2 2 2 2*/ ed[20]=d[5]*d[5] + d[12]*d[12] + d[14]*d[14] + 2 * ( d[8]*d[8] + d[9]*d[9] + d[13]*d[13] );

/*  sum_{(0 \; 3)\\} ( ad )  */
/* ( =  )  */
valuetype hd[6];
/* 0 0*/ hd[0]=ad[0] + ad[11] + ad[22];
/* 1 0*/ hd[1]=ad[1] + ad[13] + ad[24];
/* 2 0*/ hd[2]=ad[2] + ad[14] + ad[25];
/* 1 1*/ hd[3]=ad[3] + ad[16] + ad[27];
/* 2 1*/ hd[4]=ad[4] + ad[17] + ad[28];
/* 2 2*/ hd[5]=ad[5] + ad[18] + ad[29];

/*  sum_{(1 \; 2)\\} b bd  */
/*  (=sum_{(1 \; 2)\\}()())  */
valuetype id[3];
/* 0*/ id[0]=b[0]*bd[0] + b[1]*bd[1] + b[2]*bd[2];
/* 1*/ id[1]=b[1]*bd[0] + b[3]*bd[1] + b[4]*bd[2];
/* 2*/ id[2]=b[2]*bd[0] + b[4]*bd[1] + b[5]*bd[2];

/*  sum_{(1 \; 3)\\} ( cd )  */
/* ( =  )  */
valuetype jd[6];
/* 0 0*/ jd[0]=cd[0] + cd[6] + cd[11];
/* 1 0*/ jd[1]=cd[1] + cd[8] + cd[13];
/* 2 0*/ jd[2]=cd[2] + cd[9] + cd[14];
/* 1 1*/ jd[3]=cd[6] + cd[15] + cd[18];
/* 2 1*/ jd[4]=cd[7] + cd[16] + cd[19];
/* 2 2*/ jd[5]=cd[11] + cd[18] + cd[20];

/*  sum_{(3 \; 2)\\} ( ad )  */
/* ( =  )  */
valuetype kd[9];
/* 0 0*/ kd[0]=ad[0] + ad[3] + ad[5];
/* 0 1*/ kd[1]=ad[1] + ad[6] + ad[8];
/* 0 2*/ kd[2]=ad[2] + ad[7] + ad[9];
/* 1 0*/ kd[3]=ad[10] + ad[13] + ad[15];
/* 1 1*/ kd[4]=ad[11] + ad[16] + ad[18];
/* 1 2*/ kd[5]=ad[12] + ad[17] + ad[19];
/* 2 0*/ kd[6]=ad[20] + ad[23] + ad[25];
/* 2 1*/ kd[7]=ad[21] + ad[26] + ad[28];
/* 2 2*/ kd[8]=ad[22] + ad[27] + ad[29];

/*  sum_{(0 \; 3)\\} bd c  */
/*  (=sum_{(0 \; 3)\\}()())  */
valuetype ld[6];
/* 0 0*/ ld[0]=bd[0]*c[0] + bd[1]*c[1] + bd[2]*c[2];
/* 1 0*/ ld[1]=bd[0]*c[1] + bd[1]*c[3] + bd[2]*c[4];
/* 2 0*/ ld[2]=bd[0]*c[2] + bd[1]*c[4] + bd[2]*c[5];
/* 1 1*/ ld[3]=bd[0]*c[3] + bd[1]*c[6] + bd[2]*c[7];
/* 2 1*/ ld[4]=bd[0]*c[4] + bd[1]*c[7] + bd[2]*c[8];
/* 2 2*/ ld[5]=bd[0]*c[5] + bd[1]*c[8] + bd[2]*c[9];

/*  sum_{(1 \; 3)\\} dd dd  */
/*  (=sum_{(1 \; 3)\\}()())  */
valuetype od[6];
/* 0 0*/ od[0]=dd[0]*dd[0] + dd[1]*dd[1] + dd[2]*dd[2];
/* 1 0*/ od[1]=dd[0]*dd[1] + dd[1]*dd[3] + dd[2]*dd[4];
/* 2 0*/ od[2]=dd[0]*dd[2] + dd[1]*dd[4] + dd[2]*dd[5];
/* 1 1*/ od[3]=dd[1]*dd[1] + dd[3]*dd[3] + dd[4]*dd[4];
/* 2 1*/ od[4]=dd[1]*dd[2] + dd[3]*dd[4] + dd[4]*dd[5];
/* 2 2*/ od[5]=dd[2]*dd[2] + dd[4]*dd[4] + dd[5]*dd[5];

/*  sum_{(3 \; 6)(1 \; 5)\\} cd c  */
/*  (=sum_{(3 \; 6)(1 \; 5)\\}(sum_{(2 \; 5)\\}()())())  */
valuetype pd[10];
/* 0 0 0*/ pd[0]=cd[0]*c[0] + cd[6]*c[3] + cd[11]*c[5] + 2 * ( cd[1]*c[1] + cd[2]*c[2] + cd[7]*c[4] );
/* 1 0 0*/ pd[1]=cd[0]*c[1] + cd[6]*c[6] + cd[11]*c[8] + 2 * ( cd[1]*c[3] + cd[2]*c[4] + cd[7]*c[7] );
/* 2 0 0*/ pd[2]=cd[0]*c[2] + cd[6]*c[7] + cd[11]*c[9] + 2 * ( cd[1]*c[4] + cd[2]*c[5] + cd[7]*c[8] );
/* 1 1 0*/ pd[3]=cd[1]*c[1] + cd[3]*c[3] + cd[4]*c[4] + cd[6]*c[3] + cd[7]*c[4] + cd[8]*c[6] + cd[9]*c[7] + cd[12]*c[7] + cd[13]*c[8];
/* 2 1 0*/ pd[4]=cd[1]*c[2] + cd[3]*c[4] + cd[4]*c[5] + cd[6]*c[4] + cd[7]*c[5] + cd[8]*c[7] + cd[9]*c[8] + cd[12]*c[8] + cd[13]*c[9];
/* 2 2 0*/ pd[5]=cd[2]*c[2] + cd[4]*c[4] + cd[5]*c[5] + cd[7]*c[4] + cd[9]*c[7] + cd[10]*c[8] + cd[11]*c[5] + cd[13]*c[8] + cd[14]*c[9];
/* 1 1 1*/ pd[6]=cd[6]*c[1] + cd[15]*c[6] + cd[18]*c[8] + 2 * ( cd[8]*c[3] + cd[9]*c[4] + cd[16]*c[7] );
/* 2 1 1*/ pd[7]=cd[6]*c[2] + cd[15]*c[7] + cd[18]*c[9] + 2 * ( cd[8]*c[4] + cd[9]*c[5] + cd[16]*c[8] );
/* 2 2 1*/ pd[8]=cd[7]*c[2] + cd[9]*c[4] + cd[10]*c[5] + cd[12]*c[4] + cd[13]*c[5] + cd[16]*c[7] + cd[17]*c[8] + cd[18]*c[8] + cd[19]*c[9];
/* 2 2 2*/ pd[9]=cd[11]*c[2] + cd[18]*c[7] + cd[20]*c[9] + 2 * ( cd[13]*c[4] + cd[14]*c[5] + cd[19]*c[8] );

/*  sum_{(0 \; 4)\\} bd d  */
/*  (=sum_{(0 \; 4)\\}()())  */
valuetype qd[10];
/* 0 0 0*/ qd[0]=bd[0]*d[0] + bd[1]*d[1] + bd[2]*d[2];
/* 1 0 0*/ qd[1]=bd[0]*d[1] + bd[1]*d[3] + bd[2]*d[4];
/* 2 0 0*/ qd[2]=bd[0]*d[2] + bd[1]*d[4] + bd[2]*d[5];
/* 1 1 0*/ qd[3]=bd[0]*d[3] + bd[1]*d[6] + bd[2]*d[7];
/* 2 1 0*/ qd[4]=bd[0]*d[4] + bd[1]*d[7] + bd[2]*d[8];
/* 2 2 0*/ qd[5]=bd[0]*d[5] + bd[1]*d[8] + bd[2]*d[9];
/* 1 1 1*/ qd[6]=bd[0]*d[6] + bd[1]*d[10] + bd[2]*d[11];
/* 2 1 1*/ qd[7]=bd[0]*d[7] + bd[1]*d[11] + bd[2]*d[12];
/* 2 2 1*/ qd[8]=bd[0]*d[8] + bd[1]*d[12] + bd[2]*d[13];
/* 2 2 2*/ qd[9]=bd[0]*d[9] + bd[1]*d[13] + bd[2]*d[14];

/*  sum_{(1 \; 3)\\} ( ed )  */
/* ( =  )  */
valuetype rd[6];
/* 0 0*/ rd[0]=ed[0] + ed[6] + ed[11];
/* 1 0*/ rd[1]=ed[1] + ed[8] + ed[13];
/* 2 0*/ rd[2]=ed[2] + ed[9] + ed[14];
/* 1 1*/ rd[3]=ed[6] + ed[15] + ed[18];
/* 2 1*/ rd[4]=ed[7] + ed[16] + ed[19];
/* 2 2*/ rd[5]=ed[11] + ed[18] + ed[20];

/*  sum_{(1 \; 0)\\} ( ed )  */
/* ( =  )  */
valuetype sd[6];
/* 0 0*/ sd[0]=ed[0] + ed[3] + ed[5];
/* 1 0*/ sd[1]=ed[1] + ed[8] + ed[10];
/* 2 0*/ sd[2]=ed[2] + ed[12] + ed[14];
/* 1 1*/ sd[3]=ed[3] + ed[15] + ed[17];
/* 2 1*/ sd[4]=ed[4] + ed[16] + ed[19];
/* 2 2*/ sd[5]=ed[5] + ed[17] + ed[20];

/*  sum_{(0 \; 3)\\} id c  */
/*  (=sum_{(0 \; 3)\\}(sum_{(1 \; 2)\\}()())())  */
valuetype vd[6];
/* 0 0*/ vd[0]=id[0]*c[0] + id[1]*c[1] + id[2]*c[2];
/* 1 0*/ vd[1]=id[0]*c[1] + id[1]*c[3] + id[2]*c[4];
/* 2 0*/ vd[2]=id[0]*c[2] + id[1]*c[4] + id[2]*c[5];
/* 1 1*/ vd[3]=id[0]*c[3] + id[1]*c[6] + id[2]*c[7];
/* 2 1*/ vd[4]=id[0]*c[4] + id[1]*c[7] + id[2]*c[8];
/* 2 2*/ vd[5]=id[0]*c[5] + id[1]*c[8] + id[2]*c[9];

/*  sum_{(0 \; 2)\\} bd ld  */
/*  (=sum_{(0 \; 2)\\}()(sum_{(0 \; 3)\\}()()))  */
valuetype ge[3];
/* 0*/ ge[0]=bd[0]*ld[0] + bd[1]*ld[1] + bd[2]*ld[2];
/* 1*/ ge[1]=bd[0]*ld[1] + bd[1]*ld[3] + bd[2]*ld[4];
/* 2*/ ge[2]=bd[0]*ld[2] + bd[1]*ld[4] + bd[2]*ld[5];

/*  sum_{(1 \; 0)\\} ( b )  */
/* ( =  )  */
/**/ M[0]=b[0] + b[3] + b[5];

/*  sum_{(0 \; 1)\\} ( zc )  */
/* ( =  )  */
/**/ M[1]=zc[0] + zc[3] + zc[5];

/*  sum_{(1 \; 3)(0 \; 2)\\} zc b  */
/*  (=sum_{(1 \; 3)(0 \; 2)\\}(sum_{(1 \; 3)\\}()())())  */
/**/ M[2]=zc[0]*b[0] + zc[3]*b[3] + zc[5]*b[5] + 2 * ( zc[1]*b[1] + zc[2]*b[2] + zc[4]*b[4] );

/*  sum_{(1 \; 3)(0 \; 2)\\} zc hd  */
/*  (=sum_{(1 \; 3)(0 \; 2)\\}(sum_{(1 \; 3)\\}()())())  */
/**/ M[3]=zc[0]*hd[0] + zc[3]*hd[3] + zc[5]*hd[5] + 2 * ( zc[1]*hd[1] + zc[2]*hd[2] + zc[4]*hd[4] );

/*  sum_{(1 \; 3)(0 \; 2)\\} vd b  */
/*  (=sum_{(1 \; 3)(0 \; 2)\\}(sum_{(0 \; 3)\\}(sum_{(1 \; 2)\\}()())())())  */
/**/ M[4]=vd[0]*b[0] + vd[3]*b[3] + vd[5]*b[5] + 2 * ( vd[1]*b[1] + vd[2]*b[2] + vd[4]*b[4] );

/*  sum_{(1 \; 2)(0 \; 3)\\} b kd  */
/*  (=sum_{(1 \; 2)(0 \; 3)\\}()())  */
/**/ M[5]=b[0]*kd[0] + b[1]*kd[1] + b[1]*kd[3] + b[2]*kd[2] + b[2]*kd[6] + b[3]*kd[4] + b[4]*kd[5] + b[4]*kd[7] + b[5]*kd[8];

/*  sum_{(1 \; 3)(0 \; 2)\\} hd b  */
/*  (=sum_{(1 \; 3)(0 \; 2)\\}()())  */
/**/ M[6]=hd[0]*b[0] + hd[3]*b[3] + hd[5]*b[5] + 2 * ( hd[1]*b[1] + hd[2]*b[2] + hd[4]*b[4] );

/*  sum_{(0 \; 1)\\} id bd  */
/*  (=sum_{(0 \; 1)\\}(sum_{(1 \; 2)\\}()())())  */
/**/ M[7]=id[0]*bd[0] + id[1]*bd[1] + id[2]*bd[2];

/*  sum_{(1 \; 3)(0 \; 2)\\} b ld  */
/*  (=sum_{(1 \; 3)(0 \; 2)\\}()(sum_{(0 \; 3)\\}()()))  */
/**/ M[8]=b[0]*ld[0] + b[3]*ld[3] + b[5]*ld[5] + 2 * ( b[1]*ld[1] + b[2]*ld[2] + b[4]*ld[4] );

/*  sum_{(1 \; 3)(0 \; 2)\\} b jd  */
/*  (=sum_{(1 \; 3)(0 \; 2)\\}()())  */
/**/ M[9]=b[0]*jd[0] + b[3]*jd[3] + b[5]*jd[5] + 2 * ( b[1]*jd[1] + b[2]*jd[2] + b[4]*jd[4] );

/*  sum_{(0 \; 1)\\} ( kd )  */
/* ( =  )  */
/**/ M[10]=kd[0] + kd[4] + kd[8];

/*  sum_{(0 \; 3)(1 \; 2)\\} kd dd  */
/*  (=sum_{(0 \; 3)(1 \; 2)\\}()())  */
/**/ M[11]=kd[0]*dd[0] + kd[1]*dd[1] + kd[2]*dd[2] + kd[3]*dd[1] + kd[4]*dd[3] + kd[5]*dd[4] + kd[6]*dd[2] + kd[7]*dd[4] + kd[8]*dd[5];

/*  sum_{(1 \; 3)(0 \; 2)\\} hd dd  */
/*  (=sum_{(1 \; 3)(0 \; 2)\\}()())  */
/**/ M[12]=hd[0]*dd[0] + hd[3]*dd[3] + hd[5]*dd[5] + 2 * ( hd[1]*dd[1] + hd[2]*dd[2] + hd[4]*dd[4] );

/*  sum_{(0 \; 7)(3 \; 6)(2 \; 5)(1 \; 4)\\} ad d  */
/*  (=sum_{(0 \; 7)(3 \; 6)(2 \; 5)(1 \; 4)\\}(sum_{(1 \; 5)\\}()())())  */
/**/ M[13]=ad[0]*d[0] + ad[6]*d[6] + ad[9]*d[9] + ad[10]*d[1] + ad[16]*d[10] + ad[19]*d[13] + ad[20]*d[2] + ad[26]*d[11] + ad[29]*d[14] + 3 * ( ad[1]*d[1] + ad[2]*d[2] + ad[3]*d[3] + ad[5]*d[5] + ad[7]*d[7] + ad[8]*d[8] + ad[11]*d[3] + ad[12]*d[4] + ad[13]*d[6] + ad[15]*d[8] + ad[17]*d[11] + ad[18]*d[12] + ad[21]*d[4] + ad[22]*d[5] + ad[23]*d[7] + ad[25]*d[9] + ad[27]*d[12] + ad[28]*d[13] ) + 6 * ( ad[4]*d[4] + ad[14]*d[7] + ad[24]*d[8] );

/*  sum_{(1 \; 0)\\} ( ld )  */
/* ( =  )  */
/**/ M[14]=ld[0] + ld[3] + ld[5];

/*  sum_{(0 \; 1)\\} ( jd )  */
/* ( =  )  */
/**/ M[15]=jd[0] + jd[3] + jd[5];

/*  sum_{(0 \; 1)\\} bd ge  */
/*  (=sum_{(0 \; 1)\\}()(sum_{(0 \; 2)\\}()(sum_{(0 \; 3)\\}()())))  */
/**/ M[16]=bd[0]*ge[0] + bd[1]*ge[1] + bd[2]*ge[2];

/*  sum_{(1 \; 3)(0 \; 2)\\} ld jd  */
/*  (=sum_{(1 \; 3)(0 \; 2)\\}(sum_{(0 \; 3)\\}()())())  */
/**/ M[17]=ld[0]*jd[0] + ld[3]*jd[3] + ld[5]*jd[5] + 2 * ( ld[1]*jd[1] + ld[2]*jd[2] + ld[4]*jd[4] );

/*  sum_{(1 \; 3)(0 \; 2)\\} jd jd  */
/*  (=sum_{(1 \; 3)(0 \; 2)\\}()())  */
/**/ M[18]=jd[0]*jd[0] + jd[3]*jd[3] + jd[5]*jd[5] + 2 * ( jd[1]*jd[1] + jd[2]*jd[2] + jd[4]*jd[4] );

/*  sum_{(2 \; 5)(1 \; 4)(0 \; 3)\\} pd c  */
/*  (=sum_{(2 \; 5)(1 \; 4)(0 \; 3)\\}(sum_{(3 \; 6)(1 \; 5)\\}(sum_{(2 \; 5)\\}()())())())  */
/**/ M[19]=pd[0]*c[0] + pd[6]*c[6] + pd[9]*c[9] + 3 * ( pd[1]*c[1] + pd[2]*c[2] + pd[3]*c[3] + pd[5]*c[5] + pd[7]*c[7] + pd[8]*c[8] ) + 6*pd[4]*c[4];

/*  sum_{(2 \; 5)(1 \; 4)(0 \; 3)\\} qd c  */
/*  (=sum_{(2 \; 5)(1 \; 4)(0 \; 3)\\}(sum_{(0 \; 4)\\}()())())  */
/**/ M[20]=qd[0]*c[0] + qd[6]*c[6] + qd[9]*c[9] + 3 * ( qd[1]*c[1] + qd[2]*c[2] + qd[3]*c[3] + qd[5]*c[5] + qd[7]*c[7] + qd[8]*c[8] ) + 6*qd[4]*c[4];

/*  sum_{(1 \; 0)\\} ( dd )  */
/* ( =  )  */
/**/ M[21]=dd[0] + dd[3] + dd[5];

/*  sum_{(0 \; 1)\\} ( od )  */
/* ( =  )  */
/**/ M[22]=od[0] + od[3] + od[5];

/*  sum_{(0 \; 1)\\} ( rd )  */
/* ( =  )  */
/**/ M[23]=rd[0] + rd[3] + rd[5];

/*  sum_{(1 \; 3)(0 \; 2)\\} od dd  */
/*  (=sum_{(1 \; 3)(0 \; 2)\\}(sum_{(1 \; 3)\\}()())())  */
/**/ M[24]=od[0]*dd[0] + od[3]*dd[3] + od[5]*dd[5] + 2 * ( od[1]*dd[1] + od[2]*dd[2] + od[4]*dd[4] );

/*  sum_{(1 \; 3)(0 \; 2)\\} sd dd  */
/*  (=sum_{(1 \; 3)(0 \; 2)\\}()())  */
/**/ M[25]=sd[0]*dd[0] + sd[3]*dd[3] + sd[5]*dd[5] + 2 * ( sd[1]*dd[1] + sd[2]*dd[2] + sd[4]*dd[4] );

/*  sum_{(1 \; 3)(0 \; 2)\\} dd rd  */
/*  (=sum_{(1 \; 3)(0 \; 2)\\}()())  */
/**/ M[26]=dd[0]*rd[0] + dd[3]*rd[3] + dd[5]*rd[5] + 2 * ( dd[1]*rd[1] + dd[2]*rd[2] + dd[4]*rd[4] );

/*  sum_{(1 \; 7)(3 \; 6)(2 \; 5)(0 \; 4)\\} ed d  */
/*  (=sum_{(1 \; 7)(3 \; 6)(2 \; 5)(0 \; 4)\\}(sum_{(3 \; 7)(2 \; 6)\\}()())())  */
/**/ M[27]=ed[0]*d[0] + ed[15]*d[10] + ed[20]*d[14] + 2 * ( ed[3]*d[3] + ed[5]*d[5] + ed[17]*d[12] ) + 4 * ( ed[1]*d[1] + ed[2]*d[2] + ed[4]*d[4] + ed[6]*d[3] + ed[8]*d[6] + ed[10]*d[8] + ed[11]*d[5] + ed[12]*d[7] + ed[14]*d[9] + ed[16]*d[11] + ed[18]*d[12] + ed[19]*d[13] ) + 8 * ( ed[7]*d[4] + ed[9]*d[7] + ed[13]*d[8] );

}

#undef MYMOMENTSET

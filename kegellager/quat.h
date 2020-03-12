/*
	This file is part of example code using interval arithmetics for implicit surface triangulation

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

#ifndef QUAT_H
#define QUAT_H

#pragma once

template<class T,class S,class U>
void mulQuat(T*a,const S*b,const U *c)
{
	a[3]= b[3]*c[3] -d3prod(b,c);
	d3kreuz(a,b,c);
	d3op2(,a,+=,c,*b[3]);
	d3op2(,a,+=,b,*c[3]);
}



template<class T,class S>
void rotateUsingQuat(T* ret,const S*q)
{
	T
		buf[4]={ret[0],ret[1],ret[2],0},
		buf2[4];
	float
		qbar[4]={-q[0],-q[1],-q[2],q[3]};
	
	mulQuat(buf2,q,buf);
	mulQuat(buf,buf2,qbar);
	d3op2(,ret,=,buf,);
}
template<class T>
void normQuat(T* quat)
{
	float sc= 1.0/sqrt( d3prod(quat,quat)+quat[3]*quat[3]);
	d3op1(,quat,*=sc);
	quat[3]*=sc;
}
#endif

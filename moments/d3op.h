/*
 * a simple macro library for the most used operations on 3d vectors.
 *
 *
 *  the operations
 *  d3op<n>(opsym0,a,... )
 *  are given so one can but together expressions three times with the arguments between the commas,
 *  where <n> gives the number of vector arguments.
 *  <br>
 *  Example: <br>
 *  d3op2(,a, *=  ,b,  *5)<br>
 *  expands to
 *  a[0] *= b[0]*5;
 *  a[1] *= b[1]*5;
 *  a[2] *= b[2]*5;
 *
 *
 *
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


#ifndef D3OP_H
#define D3OP_H
#pragma once

#define d3IsInside(max,min,p)\
 (min[0]<=p[0]&&p[0]<=max[0]  \
  &&min[1]<=p[1]&&p[1]<=max[1] \
  &&min[2]<=p[2]&&p[2]<=max[2])

#define d3kreuz(a,b,c)\
a[0]=b[1]*c[2]-b[2]*c[1];\
a[1]=b[2]*c[0]-b[0]*c[2];\
a[2]=b[0]*c[1]-b[1]*c[0];

#define d3prod(a,b)\
(a[0]*b[0]+a[1]*b[1]+a[2]*b[2])

#define d3out(p)\
'('<<p[0]<<' '<<p[1]<<' '<<p[2]<<')'

#define d3norm(x)\
sqrt(d3prod(x,x))

#define d3op5(opsym0,a,opsym1,b,opsym2,c,opsym3,d,opsym4,e,opsym5)\
{opsym0 a[0] opsym1 b[0] opsym2 c[0] opsym3 d[0] opsym4 e[0] opsym5;\
opsym0 a[1] opsym1 b[1] opsym2 c[1] opsym3 d[1] opsym4 e[1] opsym5;\
opsym0 a[2] opsym1 b[2] opsym2 c[2] opsym3 d[2] opsym4 e[2] opsym5;}

#define d3op4(opsym0,a,opsym1,b,opsym2,c,opsym3,d,opsym4)\
{opsym0 a[0] opsym1 b[0] opsym2 c[0] opsym3 d[0] opsym4;\
opsym0 a[1] opsym1 b[1] opsym2 c[1] opsym3 d[1] opsym4;\
opsym0 a[2] opsym1 b[2] opsym2 c[2] opsym3 d[2] opsym4;}

#define d3op3(opsym0,a,opsym1,b,opsym2,c,opsym3)\
{opsym0 a[0] opsym1 b[0] opsym2 c[0] opsym3;\
opsym0 a[1] opsym1 b[1] opsym2 c[1] opsym3;\
opsym0 a[2] opsym1 b[2] opsym2 c[2] opsym3;}

#define d3op2(opsym0,a,opsym1,b,opsym2)\
{opsym0 a[0] opsym1 b[0] opsym2;\
opsym0 a[1] opsym1 b[1] opsym2;\
opsym0 a[2] opsym1 b[2] opsym2;}

#define d3op1(opsym0,a,opsym1)\
{opsym0 a[0] opsym1;\
opsym0 a[1] opsym1;\
opsym0 a[2] opsym1;}

#endif

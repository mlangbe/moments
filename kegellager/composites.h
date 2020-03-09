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

#ifndef COMPOSITES_H
#define COMPOSITES_H

#include "Objects.h"

struct KegelKegel
{

	Kegel kegel;
	

	Cylinder subTop;
	Cylinder subBottom;
	Intersection s1,s2;
	Rotated x;
	ArrangedInCircle aic;
	Translated ret;
	
	KegelKegel(float biga,float smalla,int n,float minz,float maxz,float cz,float off=0,float tanaend=0)
		:kegel(0,tan(smalla)),
		subBottom(minz*tan(smalla)*.25, minz*tan(smalla)*.25+tanaend*(maxz-minz),minz,maxz),
		subTop(maxz*tan(smalla)*.25 +tanaend*(maxz-minz),maxz*tan(smalla)*.25,minz,maxz),
		s1(kegel,subBottom),s2(s1,subTop),
		x(s2,0,sin(biga*.5),0,cos(biga*.5)),aic(x,n,off),ret(aic,0,0,cz)
	{
	}

	operator Implicit& (){return ret;}
};

#endif
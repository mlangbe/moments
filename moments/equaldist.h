/*
	equal distribution over triangulated surface.

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

#pragma once
#ifndef EQUALDIST_H
#define EQUALDIST_H
#include<vector>
#include<cmath>
#include "d3op.h"
#include "momentSet.h"
typedef faceMomSet::valuetype mynum;



class equaldist
{
public:
	std::vector<mynum> points;
	std::vector<int> tris;
	std::vector<int> quads;
	std::vector<std::vector<int> > polygons;




	inline void extract(int id,mynum x[3])
	{
		int id3=id*3;	
		x[0]=points[id3];x[1]=points[id3+1];x[2]=points[id3+2];

	}

	inline mynum area(int i1,int i2,int i3)
	{
		mynum g[3];
		areaNormal(g,i1,i2,i3);
		return sqrt(d3prod(g,g))/2;
	}

	inline void areaNormal(mynum g[3],int i1,int i2,int i3)
	{
		mynum a[3],b[3],c[3];
		extract(i1,a);extract(i2,b);extract(i3,c);

		d3op2(,b,-=,a,);
		d3op2(,c,-=,a,);
		d3kreuz(g,b,c);
	}

	mynum area();

	inline void createEqualDist(std::vector<mynum>&points,int numpts)
	{
		createEqualDensity(points,((mynum)numpts)/area());
	}
	void createEqualDensity(std::vector<mynum>&points,mynum density);

private:
	void addPointsForTri(std::vector<mynum>&points,mynum density,int i1,int i2,int i3);


};











#endif

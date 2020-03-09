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
#include<stdlib.h>
#include<vector>
#include<cmath>

#include "d3op.h"
#include "momentSet.h"
#include "equaldist.h"	

using namespace std;

mynum equaldist::area()
{
	mynum ret=0;

	for(int i=0;i<tris.size();i+=3)
	{
		ret+=area(tris[i],tris[i+1],tris[i+2]);
	}
	for(int i=0;i<quads.size();i+=4)
	{
		ret+=area(quads[i],quads[i+1],quads[i+2]);
		ret+=area(quads[i],quads[i+2],quads[i+3]);
	}
	for(int i=0;i<polygons.size();++i)
	{
		vector<int>&p=polygons[i];

		mynum vec[3]={0,0,0};mynum g[3];
		for(int j=2;j<p.size();++j)
		{
			areaNormal(g,p[0],p[j-1],p[j]);
			d3op2(,vec,+=,g,);	
		}

		ret += sqrt(d3prod(vec,vec))/2;
	}
	return ret;
}


void equaldist::addPointsForTri(std::vector<mynum>&points,mynum density,int i1,int i2,int i3)
{
		mynum num = density*area(i1,i2,i3);

		int n=(int)floor(num);
		//force expectance of n to be num
		if( rand() < RAND_MAX *( num-n) )
			++n;

		mynum a[3],b[3],c[3],ret[3];
		extract(i1,a);
		extract(i2,b);
		extract(i3,c);
		
		for(int i=0;i<n;++i)
		{
			mynum u = rand()*1.0/RAND_MAX;
			mynum v = rand()*1.0/RAND_MAX;
			mynum w=1-u-v;
			if(w<0){
				u=1-u;
				v=1-v;
				w=-w;
			}
			d3op4(,ret,=u*,a,+v*,b,+w*,c,);

			points.push_back(ret[0]);
			points.push_back(ret[1]);
			points.push_back(ret[2]);
		}
}
void equaldist::createEqualDensity(std::vector<mynum>&points,mynum density)
{
	mynum ret=0;
	points.clear();
	for(int i=0;i<tris.size();i+=3)
	{
		addPointsForTri(points,density,tris[i],tris[i+1],tris[i+2]);
	}
	for(int i=0;i<quads.size();i+=4)
	{
		addPointsForTri(points,density,quads[i],quads[i+1],quads[i+2]);
		addPointsForTri(points,density,quads[i],quads[i+2],quads[i+3]);
	}
	for(int i=0;i<polygons.size();++i)
	{
		//todo:triangulate and process.
	}
}

/*
    I/O for osg file format

	This file is part of the source code used for the calculation of the moment invariants
	as described in the dissertation of Max Langbein
	https://nbn-resolving.org/urn:nbn:de:hbz:386-kluedo-38558

	Copyright (C) 2020 TU Kaiserslautern,
	              Prof.Dr. Hans Hagen (AG Computergraphik und HCI) (hagen@cs.uni-kl.de)

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

#define _USE_MATH_DEFINES
#include <cmath>
#include<string>
#include "osgio.h"
#include "cluster.h"
using namespace std;

//create set of circles intersecting each other
void osgcircball(std::ostream&out,int num)
{
	out<<
		"Group {\n"
		" UniqueID ballid\n"
		" MatrixTransform {\n"
		"  Matrix { \n"
		"	  1 0 0 0 \n"
		"	  0 1 0 0 \n"
		"	  0 0 1 0 \n"
		"	  0 0 0 1 \n"
		"	 } \n"
		"  Geode { \n"
		"   UniqueID circleid\n"
		"   Geometry {\n"
		"    PrimitiveSets 1 {\n"
    "     DrawArrays LINE_LOOP 0 "<<num<<"\n"
    "    }\n"
		"    VertexArray Vec3Array "<<num<<" {\n"
		;
	for(int i=0;i<num;++i)
	{
		double w=i* M_PI *2/num;
		out<<"    "<<sin(w)<<' '<<cos(w)<<" 0 \n";	
	}
	out<<
		"    }\n"
		"   }\n"
		"  }\n"
		" }\n"
		;

	int num2=24;
	for(int i=0;i<num2;++i){
		double w=i*M_PI*2/num2;
		double s=sin(w),c=cos(w);
		out<<
		" MatrixTransform {\n"
		"  Matrix { \n"
		"	  "<<c<<" "<<s<<" 0 0 \n"
		"	  0 0 1 0 \n"
		"	  "<<-s<<" "<<c<<" 0 0 \n"
		"	  0 0 0 1 \n"
		"	 } \n"
		"	 Geode Use circleid \n"
		" }\n";
	}

	int num3=12;
	for(int i=0;i<num3;++i){
		double w=i*M_PI/num3;
		double s=sin(w),c=cos(w);
		out<<
		" MatrixTransform {\n"
		"  Matrix { \n"
		"	  "<<s<<" 0 0 0 \n"
		"	  0 "<<s<<" 0 0 \n"
		"	  0 0 1 0 \n"
		"	  0 0 "<<c<<" 1 \n"
		"	 } \n"
		"	 Geode Use circleid \n"
		" }\n";
	}
	
	out<<	"}\n"
		;

}


void osgpointcloud(ostream&out,const vector<mynum> pts)
{
	int n=pts.size()/3;
	out<<		
		"Geode {\n"
		" StateSet { Point { size 1 }  }\n"
		" Geometry{\n"
		"  PrimitiveSets 1 {\n"
		"   DrawArrays POINTS 0 "<<n<<"\n" 
		"  }\n" 
		"  VertexArray Vec3Array "<<n<<" {\n" 
		;
	vector<mynum>::const_iterator it=pts.begin();
	for(int i=0;i<n;++i,it+=3)
		out<<it[0]<<' '<<it[1]<<' '<<it[2]<<'\n';

	out<<
		"  }\n" 
		" }\n" 
		"}\n" 
		;
}


void osgplaceandscalemat(ostream&out,const mynum*center,
														 mynum s
														 )
{
	out<<
		"  Matrix { \n"
		"	  "<<s<<" 0 0 0 \n"
		"	  0 "<<s<<" 0 0 \n"
		"	  0 0 "<<s<<" 0 \n"
		"	  "<<center[0]<<" "<<center[1]<<" "<<center[2]<<" 1 \n"
		"	 } \n"
	;
}

void osgplaceandscale(ostream&out,const mynum*center,
														 mynum s,
														const char*id)
{
	out<<
		" MatrixTransform {\n"
		;
	osgplaceandscalemat(out,center,s);
	out<<
		"	 Node Use "<<id<<"\n"
		" }\n"
		;

}

void osgconnector(ostream&out,
												 const mynum*center,
												 const mynum*center2
												 )
{
	out<<
		" Geode {\n"
		"  Geometry{\n"
		"	  PrimitiveSets 1 { DrawArrays LINES 0 2 } \n"
		"   VertexArray Vec3Array 2 {\n"
		<<center[0]<<" "<<center[1]<<" "<<center[2]<<"  \n"
		<<center2[0]<<" "<<center2[1]<<" "<<center2[2]<<"  \n"
		"   }\n"
		"  }\n"
		" }\n"
		;
}

void osgreadpoints(vector<mynum>&pts,std::istream& in)
{

	string buf;

	in>>buf;
	while(buf != "VertexArray")
		in>>buf;

	in>>buf;
	if(buf!="Vec3Array")
		return;

	int n;

	in>>n;
	pts.reserve(n*3);
	in>>buf;
	if(buf!="{")
		return;

	for(;;){
		mynum x,y,z;
		in>>x>>y>>z;
		if(in.fail())
			break;	
		pts.push_back(x);pts.push_back(y);pts.push_back(z);
	}

}


void points_colorcoded_osg
(const std::vector<mynum>& pts,std::vector<int> code,ostream&out)
{
	int n=pts.size()/3;



	out<<		
		"Geode {\n"
		" StateSet { \n"
		"   Point { size 4 }\n"
		" }\n"
		" Geometry{\n"
		"  PrimitiveSets 1 {\n"
		"   DrawArrays POINTS 0 "<<n<<"\n" 
		"  }\n" 
		"  VertexArray Vec3Array "<<n<<" {\n" 
		;
	vector<mynum>::const_iterator it=pts.begin();
	out.precision(15);
	for(int i=0;i<n;++i,it+=3)
		out<<it[0]<<' '<<it[1]<<' '<<it[2]<<'\n';

	out<<
		"  }\n"
		"  ColorArray Vec3Array "<<n<<" {\n" 
		;
	for(int i=0;i<n;++i)
	{
		int c=code[i];
		out
			<< .25 + .25*(c&3)<<' '
			<< .25 + .25*((c>>2)&3)<<' '
			<< .25 + .25*((c>>4)&3)<<'\n';
	}
	out<<
		"  }\n"
		" }\n" 
		"}\n" 
		;
}


//points colored accordin to distance to cluster
void points_colored_momdist_osg
(const std::vector<mynum>& pts,const std::vector<mynum>& moms,
 const std::vector<mynum>& clustercenters,int momdim,ostream&out)
{
	int n=pts.size()/3;



	out<<		
		"Geode {\n"
		" StateSet { \n"
		"   Point { size 4 }\n"
		" }\n"
		" Geometry{\n"
		"  PrimitiveSets 1 {\n"
		"   DrawArrays POINTS 0 "<<n<<"\n" 
		"  }\n" 
		"  VertexArray Vec3Array "<<n<<" {\n" 
		;
	vector<mynum>::const_iterator it=pts.begin();
	out.precision(15);
	for(int i=0;i<n;++i,it+=3)
		out<<it[0]<<' '<<it[1]<<' '<<it[2]<<'\n';

	out<<
		"  }\n"
		"  ColorArray Vec3Array "<<n<<" {\n" 
		;
	it=moms.begin();
	int k=clustercenters.size()/momdim;
	cluster::eukliddist<mynum> momdist;
	vector<mynum> dists(k);mynum sumdist;
	for(int i=0;i<n;++i,it+=momdim)
	{
		vector<mynum>::const_iterator cit=clustercenters.begin();
		sumdist=0;
		for(int j=0;j<k;++j,cit+=momdim)
		{
			dists[j]=momdist(&*it,&*cit,momdim);
			sumdist+=dists[j];
		}

		//make weighting factors for colors out of dists
		for(int j=0;j<k;++j)
			dists[j]=1 - dists[j] / sumdist;		

		out	<< dists[0] <<' '<<dists[1]<<' '<<dists[2]<<'\n'<<endl;
	}
	out<<
		"  }\n"
		" }\n" 
		"}\n" 
		;
}


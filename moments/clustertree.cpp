/*
	some simple clustering algorithms.

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
#include "clustertree.h"

#include<iostream>
#include<fstream>
#include<string>
#include<cmath>
using namespace std;


clustertree::clustertree(float*distmat,int dimm)
	:dim(dimm)
{
	int dd=dim*dim;

	//the representative id for a certain leaf/tree id
	//in the matrix
	int *noderep=new int[dim];
	tree=new node[dim-1];

	//value to mark something invalid
	float big=1e36;

	for(int i=0;i<dim;++i)
		noderep[i]=i;

	for(int n=0;n<dim-1;++n)
	{
		node&nd = tree[n];
		float mind=big;
		int mini=0,minj=0;
		
		
		float*dm=distmat+dim;
		//get smallest distance between two nodes/leaves
		for(int i=1;i<dim;++i,dm+=dim)
		{
			if(dm[0]<big)
				for(int j=0;j<i;++j)
					if(dm[j]<mind)
						{ mind = dm[j]; mini=i; minj=j; }
		}
		
		nd.left= noderep[minj];
		nd.right= noderep[mini];
		nd.distance = mind;
		
		nd.number=0;
		if(nd.left<0)
			nd.number+=tree[-nd.left-1].number;
		else
			nd.number++;

		if(nd.right<0)
			nd.number+=tree[-nd.right-1].number;
		else
			nd.number++;


		//nodes are represented by negative ids
		noderep[minj]= - n - 1;
		//this id shouldn't be under consideration any more
		noderep[mini]=0xffffffff;
		
		//calculate distances to current cluster
		//and write it to row or col minj
		//and clean mini row,col
		float *dmk  = distmat,
					*dmmini = distmat + mini*dim, 
					*dmminj = distmat + minj*dim;		

		for(int k=0;k<dim;++k,dmk+=dim)
		{
			float * jel,*iel;

			if(k==minj)
			{	dmmini[k]=big; continue;}
			else if(k==mini)
			{ dmk[minj]=big;	continue;}

			
			//the bigger coordinate always denotes the row
			if(k<minj)
				jel=&dmminj[k];
			else
				jel=&dmk[minj];

			if(k<mini)
				iel =&dmmini[k];
			else
				iel =&dmk[mini];

			
			if(*iel<*jel)
				*jel=*iel;

			*iel=big;					
		}

	}

	delete[] noderep;
}

clustertree::~clustertree()
{
	{delete[]tree;}
}



//draw tree with dimensions clustered,
//on a 100x100 rectangle
void drawCorrTreePS(const float*cov,int dim,ostream&out,const string*dimnames)
{
	
	float*distmat=new float[dim*dim];

	for(int i=0;i<dim*dim;++i)
		distmat[i]=1.0-fabs(cov[i]);

	for(int i=0;i<dim;++i,cout<<endl)
		for(int j=0;j<dim;++j)
			cout<<distmat[i*dim+j]<<' ';

	clustertree ct(distmat,dim);

	for(int i=0;i<dim-1;++i)
	{
		clustertree::node&n=ct.tree[i];
		cout<<"distance:" <<n.distance
			<<" children:"<<n.left<<" "<<n.right
			<<" #leaves:"<<n.number<<endl;
	}



	delete[] distmat;

	struct rec{
		const clustertree &t;
		const string*nam;
		ostream&out;
		float w,rh;

		rec(const clustertree &tt,const string*dimnames,ostream&outt)
			:t(tt),nam(dimnames),out(outt)
		{		
			w=100,rh=100./t.dim;
			out<<" 0 0 moveto "<<w<<" 0 lineto stroke "<<endl;
			out<<" 0 "<<rh*t.dim<<" moveto "<<w<<" "<<rh*t.dim<<" lineto stroke"<<endl;
			recurse(t.dim-2,0);
		}


		void recurse( int nodeid ,float height)
		{
			const clustertree::node & n = t.tree[nodeid];
			
			float nx = w * (1-n.distance);

			float lx,ly,rx,ry;
			int lnumber=1;
			if(n.left<0){
				clustertree::node & nn = t.tree[-n.left-1];
				lx = w * (1-nn.distance);
				lnumber=nn.number;
				ly = height+rh*lnumber/2;
			}
			else{
				lx = w;
				ly=height+rh/2;						
			}

			if(n.right<0)
			{
				clustertree::node & nn = t.tree[-n.right-1];
				rx = w * (1-nn.distance);
				ry = height + rh*lnumber + rh*nn.number/2;							
			}
			else
			{
				rx = w;
				ry = height + rh*lnumber+rh/2;
			}

			out<<lx<<' '<<ly<<" moveto ";
			out<<nx<<' '<<ly<<" lineto ";
			out<<nx<<' '<<ry<<" lineto ";
			out<<rx<<' '<<ry<<" lineto stroke "<<endl;

			if(n.left<0)
				recurse(-n.left-1,height);
			else
			{ 
				out<<lx<<' '<<ly-rh/4<<" moveto ";
				out
					<<"(" << n.left;
				if(nam)
					out<<":"<<nam[n.left];
				out<<") show"<<endl;

			}

			if(n.right<0)				
				recurse(-n.right-1,height + rh * lnumber);			
			else
			{  
				out<<rx<<' '<<ry-rh/4<<" moveto ";
				out
					<<"(" << n.right;
				if(nam)
					out<<":"<<nam[n.right];
				out<<") show"<<endl;			
			}

		}	
	} recu(ct,dimnames,out);

}




inline void testdrawCorrTreePS()
{
	static const int dim=5;

	ofstream out("covtree.ps");
	float mat[dim*dim]={
		6,6,6,6,6,
		1,6,6,6,6,
		0,0,6,6,6,
		0,0,0.7,6,6,
		0,0,0,0.5,6
	};

	out<<
		"/Times-Roman findfont 10 scalefont setfont "
		"10 10 translate "
		<<endl;
	drawCorrTreePS(mat,dim,out);

	out.close();

}



void testclustertree()
{
	/*
	static const int dim=5;
	float mat[dim*dim]={
		6,6,6,6,6,
		0.1,6,6,6,6,
		1,0.5,6,6,6,
		1,1,0.2,6,6,
		0.4,1,1,0.3,6
	};


	clustertree ct(mat,dim);

	for(int i=0;i<dim-1;++i)
	{
		clustertree::node&n=ct.tree[i];
		cout<<"distance:" <<n.distance
			<<" children:"<<n.left<<" "<<n.right
			<<" #leaves:"<<n.number<<endl;
	}
	
	for(int i=0;i<dim;++i,cout<<endl)
		for(int j=0;j<dim;++j)
			cout<<mat[i*dim+j]<<' ';
	*/
	testdrawCorrTreePS();
}

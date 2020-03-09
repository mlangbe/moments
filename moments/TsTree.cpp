/*
	octree with tensor set in every node used to accelerate the summing up of moment Tensors
	of the points or field inside a  ball mask.

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

#include "TsTree_impl.h"



int TsTree::numbersummed=0;
int TsTree::ptssummed=0;
/**
get All tsets of the blocks in a certain height
Action is a functor which takes three integers as input
*/
void TsTree::getBallTsets(
													std::vector<TsTree::tset> &inner,
													std::vector<TsTree::tset> &border,
													int minheight,int radius
													) const
{
	
	//impl::doforcenters(*this,minheight>0?minheight:0, inner,border, minheight,radius);

	impl::treeWalker<impl::getBallsForCentersAction> tw(*this);
	tw.doit.hmincenters = minheight>0?minheight:0;
	tw.doit.hmincalc=minheight;
	tw.doit.radius = radius*getScale();
	tw.walk();
	inner.swap(tw.doit.inner);
	border.swap(tw.doit.border);
}


/**
helper class for constructor:
build the tree from the buckets array,
build the off array and replace the numbers in buckets by ids in the off array
*/
TsTree::TsTree(const std::vector<long double> & pts,long double boxwidth,int method)
{

	assert(pts.size()%3 ==0);

	//calculate bbox, variance

	vector<long double>::const_iterator it=pts.begin();

	for(int i=0;i<3;++i,++it)
		min[i]=max[i]=*it,
		avg[i]=var[i]=0;

	it=pts.begin();
	while(it<pts.end())
	{
		for(int i=0;i<3;++i,++it)
		{
			if(min[i]>*it)min[i]=*it;
			if(max[i]<*it)max[i]=*it;
			avg[i]+=*it;
			var[i]+=*it * *it;
		}		
	}

	for(int i=0;i<3;++i)
	{
		avg[i]/= (pts.size()/3);
		var[i]/= (pts.size()/3);
		var[i]-=avg[i]*avg[i];
	}

	if(method==0)
		impl::buildFromPoints_old(*this,pts,boxwidth);
	else
		impl::buildFromPoints(*this,pts,boxwidth);

}



void TsTree::sumInCircle(
												 tset&innerSum,tset&borderSum,
												 const int center[3], int radius, int minh
												 )const
{
	impl::sumCircle2 s(center,radius,*this,minh);
	numbersummed+=s.numbersummed;
	ptssummed+=s.ptssummed;
	innerSum=s.ts;
	borderSum=s.bordts; 
}


void TsTree::sumInCircle(
															tset&innerSum,tset&borderSum,
															const mynum center[3], mynum radius, int minh
															) const
{
	impl::sumFloatCircle2 s(center,radius,*this,minh);
	numbersummed+=s.numbersummed;
	ptssummed+=s.ptssummed;
	innerSum=s.ts;
	borderSum=s.bordts; 
}

	void TsTree::sumInCircle(
			tset*tsets,int n,
			const mynum center[3], 
			const mynum* radii
			) const
			
	{

    //for small numbers, the single-radius procedure is faster
		if(n<20)
		{
			tset buf;
			for(int i=0;i<n;++i)
				sumInCircle(tsets[i],buf,center,radii[i],-1);
			return;
		}
	
		TsTree::impl::treeWalker<TsTree::impl::sumInCircleMultiRadiiAction> w(*this);
		vector<mynum> rsq(n);
		w.doit.rsqbegin =&rsq[0];
		w.doit.rsqend = w.doit.rsqbegin+n;
				
		std::fill(tsets,tsets + n, tset().reset());

		for(int i=0;i<n;++i)
			w.doit.rsqbegin[i]= radii[i]*radii[i];
		
		w.doit.tsbegin=tsets;
		w.doit.c = &center[0];
		w.walk();			
	}
	



TsTree::~TsTree()
{
	if(memfornodes.size())return;
	
		
	class rec
	{

		int h;

		void recurse(TsTree::ND&nd)
		{
			--h;

			assert(nd.f <= 0xff);
			//cout<<'<'<<hex<<nd->f<<flush;
			if(h>0){
				for(int i=0,f=nd.f;f!=0;f>>=1,++i){
					if(f&1){
						recurse(*nd.ch.n[i]);
						delete nd.ch.n[i];
					}
				}
			}
			//cout<<'>'<<flush;
			++h;
			//cout<<nd<<endl;
		}

	public:

		rec(TsTree&t)
		{
			h=t.depth;
			recurse(t.root);			
		}			
	} ;

	rec k(*this);



}


long double TsTree::fractdim() const
{
	struct rec
	{

		//number of blocks in each height
		vector<int> num;

		int h;

		void recurse(const TsTree::ND*nd)
		{
			--h;

			//cout<<'<'<<hex<<nd->f<<flush;
			int f,i;
			for(i=0,f=nd->f;f!=0;f>>=1,++i)
				if(f&1)
				{
					if(h>0)
						recurse(nd->ch.n[i]);
					++num[h];
				}
				//cout<<'>'<<flush;
				++h;
		}

	public:

		rec(const TsTree&t)
			:num(t.depth,0)
		{
			h=t.depth;
			//calculate number of blocks in each height
			recurse(&t.root);			


		}			
	} ;

	//calculate number of blocks in each height
	rec nums(*this);	

	//make lin regression of log(num[x]) versus x
	long double sx=0,sy=0,sxx=0,sxy=0,
		ln2=log(2.0l),y;

	int minn=0,maxn=depth;
	for(int x=minn;x<maxn;++x)
	{
		y=log((long double)nums.num[x])/ln2;
		//cout<<y<<' '<<nums.num[x]<<endl;
		sx+=x;sy+=y;
		sxx+=x*x;
		sxy+=x*y;
	}


	return - (sxy-sy*sx/(maxn-minn)) 
		/ ( sxx - sx*sx/(maxn-minn));



}

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

#include<math.h>
#include<cassert>
#include<iostream>
#include<algorithm>
#include "TsTree.h"
#include "zindex.h"
using namespace std;




struct TsTree::impl{
	//build the tree from the information in buckets
	struct buildFromBuckets;


	//update the tset sums from vector tsets
	struct updatetsets;


//help structures for sumInCircle
	struct sumCircle;
	struct sumCircle2;
	struct sumFloatCircle;
	struct sumFloatCircle2;


	struct doforcenters;


	








struct sumCircle
{
	const TsTree&tree;
	//the tensor set where all others are added to
	tset ts;

	//the sum of the tensorsets intersected by the ball border
	tset bordts;

	//the x,y,z of the center of the current node
	//relative to the center of the ball
	int x,y,z;

	//the real values of the Ball
	TsTree::mynum center_[3];
	TsTree::mynum rsq_;

	//radius of ball
	const int r;

	//radius,squared
	const int rsq;

	//the current height (coords for current node are
	//(x>>h)&1
	int h;

	//minumum height of nodes
	const int hmin;

	//half of the edge length of the current box
	int rb;


	//number of tsets summed 
	int numbersummed;

	//number of single points summed
	int ptssummed;

	//if true: box inside ball, if false:
	//box intersects ball or is outside
	inline bool inside()
	{		
		int 
			ax= x>=0?x:-x, 
			ay= y>=0?y:-y, 
			az= z>=0?z:-z;				
		ax+=rb;ay+=rb;az+=rb;

		return (ax*ax+ay*ay+az*az <= rsq );
	}


	//if true: box outside ball, if false:
	//box intersects ball or is outside
	inline bool outside()
	{		
		int 
			ax= x>=0?x:-x, 
			ay= y>=0?y:-y, 
			az= z>=0?z:-z;

		ax-=rb;ay-=rb;az-=rb;

		//a is now the distance vector from the nearest part
		//of the box to the center
		int c= (ax>0) | ((ay>0)<<1) | ((az>0)<<2);
		switch(c)
		{
		case 0: return false;
		case 1: return ax>=r;
		case 2: return ay>=r;
		case 3: return ax*ax+ay*ay>=rsq;
		case 4: return az>=r;
		case 5: return ax*ax+az*az>=rsq;
		case 6: return ay*ay+az*az>=rsq;
		case 7:	return (ax*ax+ay*ay+az*az >= rsq );				
		default: return true;
		}

	}



	inline void addpointsinside(int leafid)
	{

		vector<TsTree::mynum>::const_iterator 
			it = tree.coords.begin() + tree.off[ leafid ] ,
			itend = tree.coords.begin() + tree.off[ leafid+1] ;	


		mynum buf,lsq;			

		ptssummed += (itend-it)/3;

		for(;it!=itend;it+=3)
		{
			lsq=rsq_;
			buf = 	it[0]-center_[0]; lsq-=buf*buf;
			buf = 	it[1]-center_[1]; lsq-=buf*buf;
			buf = 	it[2]-center_[2]; lsq-=buf*buf;
			if(lsq>0)
				ts.add3dcoords(&*it);
		}

	}


	/*calculate ts as the sum of all leaves intersecting			//the ball
	* given by |x - 0.5*center | < 0.5*radius
	*/
	sumCircle(const int *center, int radius, const TsTree & tr,
		int minh)
		:r(radius),rsq(radius*radius),hmin(minh),tree(tr)

	{

		rsq_ = tr.rb0 * radius * .5; rsq_*=rsq_;

		center_[0] = center[0] * tr.rb0  * .5  + tr.min[0];
		center_[1] = center[1] * tr.rb0  * .5  + tr.min[1];
		center_[2] = center[2] * tr.rb0  * .5  + tr.min[2];


		numbersummed=0;
		ptssummed=0;

		h=tr.depth;
		//we use doubled coordinates
		rb=1<<h;
		x= rb-center[0],
			y= rb-center[1],
			z= rb-center[2],
			bordts.reset();
		ts.reset();
		if(inside())
		{
			ts = tr.root.t;
			return;
		}

		if(outside())			
			return;

		if(h==hmin)
			bordts = tr.root.t;

		else
			recurse(&tr.root);		
	}

	void recurse(const TsTree::ND * nd)
	{
		--h;rb>>=1;
		int f=nd->f;
		//test if center is inside ball
		int d=x*x+y*y+z*z;
		if(d < rsq )
		{ //leaves are inside or intersect				
			for(int i=0;i<8;++i,f>>=1)
			{
				if(f&1)
				{
					x+= (i&1)? rb :-rb;
					y+= (i&2)? rb :-rb;
					z+= (i&4)? rb :-rb;

					if( h>0)
					{
						if(inside())						
							ts+=nd->ch.n[i]->t,++numbersummed;					
						else if(h>hmin)
							recurse(nd->ch.n[i]);							
						else
							bordts+=nd->ch.n[i]->t,++numbersummed;					
					}
					else
					{
						++numbersummed;
						if(inside())						
							ts+=tree.tsets[ nd->ch.i[i] ];					
						else
							if(hmin==0)
								bordts+=tree.tsets[ nd->ch.i[i] ];						
							else{
								addpointsinside(nd->ch.i[i]);
							}
					}
					x-= (i&1)? rb :-rb;
					y-= (i&2)? rb :-rb;
					z-= (i&4)? rb :-rb;
				}
			}					
		}
		else if( d > rsq )
		{ //leaves are outside or intersect
			for(int i=0;i<8;++i,f>>=1)
			{
				if(f&1)
				{
					x+= (i&1)? rb :-rb;
					y+= (i&2)? rb :-rb;
					z+= (i&4)? rb :-rb;
					if(h>0)
					{
						if(!outside()){
							if(h>hmin)
								recurse(nd->ch.n[i]);
							else
								bordts+=nd->ch.n[i]->t,++numbersummed;
						}
					}
					else{
						if(!outside()){
							if(hmin==0)
								bordts+=tree.tsets[ nd->ch.i[i] ],++numbersummed;
							else
								addpointsinside(nd->ch.i[i]);
						}
					}
					x-= (i&1)? rb :-rb;
					y-= (i&2)? rb :-rb;
					z-= (i&4)? rb :-rb;
				}
			}												
		}					
		else
		{
			for(int i=0;i<8;++i,f>>=1)
			{
				if(f&1)
				{
					x+= (i&1)? rb :-rb;
					y+= (i&2)? rb :-rb;
					z+= (i&4)? rb :-rb;
					if(h>0)
					{
						if(inside())
							ts+=nd->ch.n[i]->t,++numbersummed;
						else
							if(!outside()){
								if(h>hmin)
									recurse(nd->ch.n[i]);
								else
									bordts+=nd->ch.n[i]->t,++numbersummed;
							}
					}
					else{
						++numbersummed;
						if(inside())
							ts += tree.tsets[ nd->ch.i[i] ];
						else
							if(!outside())
							{
								if(hmin==0)
									bordts += tree.tsets[ nd->ch.i[i] ];
								else
									addpointsinside(nd->ch.i[i]);
							}
					}
					x-= (i&1)? rb :-rb;
					y-= (i&2)? rb :-rb;
					z-= (i&4)? rb :-rb;
				}
			}																	
		}
		++h;rb<<=1;
	}






};

/**optimized version of sumcircle*/
struct sumCircle2
{
	const TsTree&tree;
	//the tensor set where all others are added to
	tset ts;

	//the sum of the tensorsets intersected by the ball border
	tset bordts;

	//the x,y,z of the center of the current node
	//relative to the center of the ball
	int x,y,z;

	//the real values of the Ball
	TsTree::mynum center_[3];
	TsTree::mynum rsq_;

	//radius of ball
	const int r;

	//radius,squared
	const int rsq;

	//the current height (coords for current node are
	//(x>>h)&1
	int h;

	//minumum height of nodes
	const int hmin;

	//half of the edge length of the current box
	int rb;


	//number of tsets summed 
	int numbersummed;

	//number of single points summed
	int ptssummed;

	//if true: box inside ball, if false:
	//box intersects ball or is outside
	inline bool inside()
	{		
		int 
			ax= x>=0?x:-x, 
			ay= y>=0?y:-y, 
			az= z>=0?z:-z;				
		ax+=rb;ay+=rb;az+=rb;

		return (ax*ax+ay*ay+az*az <= rsq );
	}


	//if true: box outside ball, if false:
	//box intersects ball or is outside
	inline bool outside()
	{		
		int 
			ax= x>=0?x:-x, 
			ay= y>=0?y:-y, 
			az= z>=0?z:-z;

		ax-=rb;ay-=rb;az-=rb;

		//a is now the distance vector from the nearest part
		//of the box to the center
		int c= (ax>0) | ((ay>0)<<1) | ((az>0)<<2);
		switch(c)
		{
		case 0: return false;
		case 1: return ax>=r;
		case 2: return ay>=r;
		case 3: return ax*ax+ay*ay>=rsq;
		case 4: return az>=r;
		case 5: return ax*ax+az*az>=rsq;
		case 6: return ay*ay+az*az>=rsq;
		case 7:	return (ax*ax+ay*ay+az*az >= rsq );				
		default: return true;
		}

	}



	inline void addpointsinside(int leafid)
	{

		vector<TsTree::mynum>::const_iterator 
			it = tree.coords.begin() + tree.off[ leafid ] ,
			itend = tree.coords.begin() + tree.off[ leafid+1] ;	


		mynum buf,lsq;			

		ptssummed += (itend-it)/3;

		for(;it!=itend;it+=3)
		{
			lsq=rsq_;
			buf = 	it[0]-center_[0]; lsq-=buf*buf;
			buf = 	it[1]-center_[1]; lsq-=buf*buf;
			buf = 	it[2]-center_[2]; lsq-=buf*buf;
			if(lsq>0)
				ts.add3dcoords(&*it);
		}

	}


	/*calculate ts as the sum of all leaves intersecting			//the ball
	* given by |x - 0.5*center | < 0.5*radius
	*/
	sumCircle2(const int *center, int radius, const TsTree & tr,
		int minh)
		:r(radius),rsq(radius*radius),hmin(minh),tree(tr)

	{

		rsq_ = tr.rb0 * radius * .5; rsq_*=rsq_;

		center_[0] = center[0] * tr.rb0  * .5  + tr.min[0];
		center_[1] = center[1] * tr.rb0  * .5  + tr.min[1];
		center_[2] = center[2] * tr.rb0  * .5  + tr.min[2];


		numbersummed=0;
		ptssummed=0;

		h=tr.depth;
		//we use doubled coordinates
		rb=1<<h;
		x= rb-center[0],
			y= rb-center[1],
			z= rb-center[2],
			bordts.reset();
		ts.reset();
		if(inside())
		{
			ts = tr.root.t;
			return;
		}

		if(outside())			
			return;

		if(h==hmin)
			bordts = tr.root.t;

		else
			recurse(&tr.root);		
	}









	inline void dorecurse(bool ins,bool out,int i,const TsTree::ND*nd)
	{
		if(out && outside())
			return;

		if(h>0){

			if(ins && inside() )						
			{
				++numbersummed;					
				ts+=nd->ch.n[i]->t;
				return;
			}

			if(h>hmin)
				recurse(nd->ch.n[i]);
			else
				bordts+=nd->ch.n[i]->t;

		}
		else{

			if(ins && inside() )						
			{
				++numbersummed;					
				ts+=tree.tsets[nd->ch.i[i]];
				return;
			}

			if(hmin==0)
				bordts+=tree.tsets[nd->ch.i[i]];
			else
				addpointsinside(nd->ch.i[i]);
		}	
	}

	void recurse(const TsTree::ND * nd)
	{
	
		bool couldbeinside,couldbeoutside;


		--h;
		int rb2=rb;
		
		rb>>=1;
		int f=nd->f;
		//test if center is inside ball
		int d=x*x+y*y+z*z;
		couldbeinside  = d <= rsq;
		couldbeoutside = d >= rsq;

		
		x-=rb;y-=rb;z-=rb;

		if(f&1)
			dorecurse(couldbeinside,couldbeoutside,0,nd);

		x+=rb2;

		if(f&2)
			dorecurse(couldbeinside,couldbeoutside,1,nd);

		y+=rb2;

		if(f&8)
			dorecurse(couldbeinside,couldbeoutside,3,nd);

		x-=rb2;

		if(f&4)
			dorecurse(couldbeinside,couldbeoutside,2,nd);

		z+=rb2;

		if(f&0x40)
			dorecurse(couldbeinside,couldbeoutside,6,nd);

		x+=rb2;

		if(f&0x80)
			dorecurse(couldbeinside,couldbeoutside,7,nd);

		y-=rb2;

		if(f&0x20)
			dorecurse(couldbeinside,couldbeoutside,5,nd);

		x-=rb2;

		if(f&0x10)
			dorecurse(couldbeinside,couldbeoutside,4,nd);

		z-=rb;
		x+=rb;y+=rb;

		++h;rb<<=1;

	}




	




};


struct sumFloatCircle
{
	typedef TsTree::mynum mynum;
	const TsTree&tree;
	//the tensor set where all others are added to
	tset ts;

	//the sum of the tensorsets intersected by the ball border
	tset bordts;

	//the x,y,z of the center of the current node
	//relative to the center of the ball
	mynum x,y,z;

	//the center of the ball
	mynum center_[3];

	//radius of ball
	const mynum r;

	//radius,squared
	const mynum rsq;

	//the current height (coords for current node are
	//(x>>h)&1
	int h;

	//minumum height of nodes
	const int hmin;

	//half of the edge length of the current box
	mynum rb;


	//number of tsets summed 
	int numbersummed;

	//number of single points summed
	int ptssummed;

	//if true: box inside ball, if false:
	//box intersects ball or is outside
	inline bool inside()
	{		
		mynum 
			ax= x>=0?x:-x, 
			ay= y>=0?y:-y, 
			az= z>=0?z:-z;				
		ax+=rb;ay+=rb;az+=rb;

		return (ax*ax+ay*ay+az*az <= rsq );
	}


	//if true: box outside ball, if false:
	//box intersects ball or is outside
	inline bool outside()
	{		
		mynum
			ax= x>=0?x:-x, 
			ay= y>=0?y:-y, 
			az= z>=0?z:-z;

		ax-=rb;ay-=rb;az-=rb;

		//a is now the distance vector from the nearest part
		//of the box to the center
		int c= (ax>0) | ((ay>0)<<1) | ((az>0)<<2);
		switch(c)
		{
		case 0: return false;
		case 1: return ax>=r;
		case 2: return ay>=r;
		case 3: return ax*ax+ay*ay>=rsq;
		case 4: return az>=r;
		case 5: return ax*ax+az*az>=rsq;
		case 6: return ay*ay+az*az>=rsq;
		case 7:	return (ax*ax+ay*ay+az*az >= rsq );				
		default: return true;
		}

	}



	inline void addpointsinside(int leafid)
	{

		vector<TsTree::mynum>::const_iterator 
			it = tree.coords.begin() + tree.off[ leafid ] ,
			itend = tree.coords.begin() + tree.off[ leafid+1] ;	


		mynum buf,lsq;			

		ptssummed += (itend-it)/3;

		for(;it!=itend;it+=3)
		{
			lsq=rsq;
			buf = 	it[0]-center_[0]; lsq-=buf*buf;
			buf = 	it[1]-center_[1]; lsq-=buf*buf;
			buf = 	it[2]-center_[2]; lsq-=buf*buf;
			if(lsq>0)
				ts.add3dcoords(&*it);
		}

	}


	/*calculate ts as the sum of all leaves intersecting			//the ball
	* given by |x - center | < radius
	*/
	sumFloatCircle(const mynum center[3], mynum radius, const TsTree & tr,
		int minh)
		:r(radius),rsq(radius*radius),hmin(minh),tree(tr)

	{
		copy(center,center+3,center_);


		numbersummed=0;
		ptssummed=0;

		h=tr.depth;

		rb = ldexpl(tree.getScale(),h-1);
		//x,y,z are now the distancecs of the center 
		//of the circle to the center of the current box
		const mynum*min = tr.getOrig();
		x= min[0]+rb-center[0];
		y= min[1]+rb-center[1];
		z= min[2]+rb-center[2];

		bordts.reset();
		ts.reset();
		if(inside())
		{
			ts = tr.root.t;
			return;
		}

		if(outside())			
			return;

		if(h==hmin)
			bordts = tr.root.t;

		else
			recurse(&tr.root);		
	}

	void recurse(const TsTree::ND * nd)
	{
		--h;rb*=0.5;
		int f=nd->f;
		//test if center is inside ball
		mynum d=x*x+y*y+z*z;
		if(d < rsq )
		{ //leaves are inside or intersect				
			for(int i=0;i<8;++i,f>>=1)
			{
				if(f&1)
				{
					x+= (i&1)? rb :-rb;
					y+= (i&2)? rb :-rb;
					z+= (i&4)? rb :-rb;

					if( h>0)
					{
						if(inside())						
							ts+=nd->ch.n[i]->t,++numbersummed;					
						else if(h>hmin)
							recurse(nd->ch.n[i]);							
						else
							bordts+=nd->ch.n[i]->t,++numbersummed;					
					}
					else
					{
						++numbersummed;
						if(inside())						
							ts+=tree.tsets[ nd->ch.i[i] ];					
						else
							if(hmin==0)
								bordts+=tree.tsets[ nd->ch.i[i] ];						
							else{
								addpointsinside(nd->ch.i[i]);
							}
					}
					x-= (i&1)? rb :-rb;
					y-= (i&2)? rb :-rb;
					z-= (i&4)? rb :-rb;
				}
			}					
		}
		else if( d > rsq )
		{ //leaves are outside or intersect
			for(int i=0;i<8;++i,f>>=1)
			{
				if(f&1)
				{
					x+= (i&1)? rb :-rb;
					y+= (i&2)? rb :-rb;
					z+= (i&4)? rb :-rb;
					if(h>0)
					{
						if(!outside()){
							if(h>hmin)
								recurse(nd->ch.n[i]);
							else
								bordts+=nd->ch.n[i]->t,++numbersummed;
						}
					}
					else{
						if(!outside()){
							if(hmin==0)
								bordts+=tree.tsets[ nd->ch.i[i] ],++numbersummed;
							else
								addpointsinside(nd->ch.i[i]);
						}
					}
					x-= (i&1)? rb :-rb;
					y-= (i&2)? rb :-rb;
					z-= (i&4)? rb :-rb;
				}
			}												
		}					
		else
		{
			for(int i=0;i<8;++i,f>>=1)
			{
				if(f&1)
				{
					x+= (i&1)? rb :-rb;
					y+= (i&2)? rb :-rb;
					z+= (i&4)? rb :-rb;
					if(h>0)
					{
						if(inside())
							ts+=nd->ch.n[i]->t,++numbersummed;
						else
							if(!outside()){
								if(h>hmin)
									recurse(nd->ch.n[i]);
								else
									bordts+=nd->ch.n[i]->t,++numbersummed;
							}
					}
					else{
						++numbersummed;
						if(inside())
							ts += tree.tsets[ nd->ch.i[i] ];
						else
							if(!outside())
							{
								if(hmin==0)
									bordts += tree.tsets[ nd->ch.i[i] ];
								else
									addpointsinside(nd->ch.i[i]);
							}
					}
					x-= (i&1)? rb :-rb;
					y-= (i&2)? rb :-rb;
					z-= (i&4)? rb :-rb;
				}
			}																	
		}
		++h;rb*=2;
	}






};







/** optimized version of sumFloatCircle*/
struct sumFloatCircle2
{
	typedef TsTree::mynum mynum;
	const TsTree&tree;
	//the tensor set where all others are added to
	tset ts;

	//the sum of the tensorsets intersected by the ball border
	tset bordts;

	//the x,y,z of the center of the current node
	//relative to the center of the ball
	mynum x,y,z;

	//the center of the ball
	mynum center_[3];

	//radius of ball
	const mynum r;

	//radius,squared
	const mynum rsq;

	//the current height (coords for current node are
	//(x>>h)&1
	int h;

	//minumum height of nodes
	const int hmin;

	//half of the edge length of the current box
	mynum rb;


	//number of tsets summed 
	int numbersummed;

	//number of single points summed
	int ptssummed;

	//if true: box inside ball, if false:
	//box intersects ball or is outside
	inline bool inside()
	{		
		mynum 
			ax= x>=0?x:-x, 
			ay= y>=0?y:-y, 
			az= z>=0?z:-z;				
		ax+=rb;ay+=rb;az+=rb;

		return (ax*ax+ay*ay+az*az <= rsq );
	}


	//if true: box outside ball, if false:
	//box intersects ball or is outside
	inline bool outside()
	{		
		mynum
			ax= x>=0?x:-x, 
			ay= y>=0?y:-y, 
			az= z>=0?z:-z;

		ax-=rb;ay-=rb;az-=rb;

		//a is now the distance vector from the nearest part
		//of the box to the center
		int c= (ax>0) | ((ay>0)<<1) | ((az>0)<<2);
		switch(c)
		{
		case 0: return false;
		case 1: return ax>=r;
		case 2: return ay>=r;
		case 3: return ax*ax+ay*ay>=rsq;
		case 4: return az>=r;
		case 5: return ax*ax+az*az>=rsq;
		case 6: return ay*ay+az*az>=rsq;
		case 7:	return (ax*ax+ay*ay+az*az >= rsq );				
		default: return true;
		}

	}



	inline void addpointsinside(int leafid)
	{

		vector<TsTree::mynum>::const_iterator 
			it = tree.coords.begin() + tree.off[ leafid ] ,
			itend = tree.coords.begin() + tree.off[ leafid+1] ;	


		mynum buf,lsq;			

		ptssummed += (itend-it)/3;

		for(;it!=itend;it+=3)
		{
			lsq=rsq;
			buf = 	it[0]-center_[0]; lsq-=buf*buf;
			buf = 	it[1]-center_[1]; lsq-=buf*buf;
			buf = 	it[2]-center_[2]; lsq-=buf*buf;
			if(lsq>0)
				ts.add3dcoords(&*it);
		}

	}


	/*calculate ts as the sum of all leaves intersecting			//the ball
	* given by |x - center | < radius
	*/
	sumFloatCircle2(const mynum center[3], mynum radius, const TsTree & tr,
		int minh)
		:r(radius),rsq(radius*radius),hmin(minh),tree(tr)

	{
		copy(center,center+3,center_);


		numbersummed=0;
		ptssummed=0;

		h=tr.depth;

		rb = ldexpl(tree.getScale(),h-1);
		//x,y,z are now the distancecs of the center 
		//of the circle to the center of the current box
		const mynum*min = tr.getOrig();
		x= min[0]+rb-center[0];
		y= min[1]+rb-center[1];
		z= min[2]+rb-center[2];

		bordts.reset();
		ts.reset();
		if(inside())
		{
			ts = tr.root.t;
			return;
		}

		if(outside())			
			return;

		if(h==hmin)
			bordts = tr.root.t;

		else
			recurse(&tr.root);		
	}

	inline void dorecurse(bool ins,bool out,int i,const TsTree::ND*nd)
	{
		if(out && outside())
			return;

		if(h>0){

			if(ins && inside() )						
			{
				++numbersummed;					
				ts+=nd->ch.n[i]->t;
				return;
			}

			if(h>hmin)
				recurse(nd->ch.n[i]);
			else
				bordts+=nd->ch.n[i]->t;

		}
		else{

			if(ins && inside() )						
			{
				++numbersummed;					
				ts+=tree.tsets[nd->ch.i[i]];
				return;
			}

			if(hmin==0)
				bordts+=tree.tsets[nd->ch.i[i]];
			else
				addpointsinside(nd->ch.i[i]);
		}	
	}

	void recurse(const TsTree::ND * nd)
	{
	
		bool couldbeinside,couldbeoutside;


		--h;
		mynum rb2=rb;
		
		rb*=0.5;
		int f=nd->f;
		//test if center of box is inside ball
		mynum d=x*x+y*y+z*z;
		couldbeinside  = d <= rsq;
		couldbeoutside = d >= rsq;

		
		x-=rb;y-=rb;z-=rb;

		if(f&1)
			dorecurse(couldbeinside,couldbeoutside,0,nd);

		x+=rb2;

		if(f&2)
			dorecurse(couldbeinside,couldbeoutside,1,nd);

		y+=rb2;

		if(f&8)
			dorecurse(couldbeinside,couldbeoutside,3,nd);

		x-=rb2;

		if(f&4)
			dorecurse(couldbeinside,couldbeoutside,2,nd);

		z+=rb2;

		if(f&0x40)
			dorecurse(couldbeinside,couldbeoutside,6,nd);

		x+=rb2;

		if(f&0x80)
			dorecurse(couldbeinside,couldbeoutside,7,nd);

		y-=rb2;

		if(f&0x20)
			dorecurse(couldbeinside,couldbeoutside,5,nd);

		x-=rb2;

		if(f&0x10)
			dorecurse(couldbeinside,couldbeoutside,4,nd);

		z-=rb;
		x+=rb;y+=rb;

		++h;rb*=2;

	}




	
};



struct treeWalkerData
{
	treeWalkerData(const TsTree&tr)
		:tree(tr)
	{}

	typedef TsTree::mynum mynum;
	const TsTree&tree;

	//current x,y,z values of the center of the current node
	mynum center[3];

	//the current height (current node width= ldexp(tree.getscale(),h))
	//ldexp((x>>h)&1
	int h;

	//half of the edge length of the current box
	mynum rb;

};


	

template<class action>
struct treeWalker :public treeWalkerData
{
	//the action to do for node or leaves
	action doit;

	treeWalker(const TsTree & tr)
		:treeWalkerData(tr)
	{}

	void walk()
	{
		h=tree.depth;

		rb = ldexpl(tree.getScale(),h-1);
		//x,y,z are now the distancecs of the center 
		//of the circle to the center of the current box
		const mynum*min = tree.getOrig();
		center[0]= min[0]+rb;
		center[1]= min[1]+rb;
		center[2]= min[2]+rb;

		doit(*this,tree.root.t,-1,&tree.root);		
	}

	template<int i>
	inline void dorecurse(const TsTree::ND*nd)
	{
		if(nd->f & (1<<i))
		{
			if(h>0){
				doit(*this,nd->ch.n[i]->t,-1,nd->ch.n[i]);
			}
			else{
				doit( *this, tree.tsets[nd->ch.i[i]], nd->ch.i[i],0 );
			}
		}
	}

	void recurse(const TsTree::ND * nd) 
	{
		--h;
		mynum rb2=rb;
		
		rb*=0.5;	
		
		center[0]-=rb;center[1]-=rb;center[2]-=rb;

		dorecurse<0>(nd);
		center[0]+=rb2;
		dorecurse<1>(nd);
		center[1]+=rb2;
		dorecurse<3>(nd);
		center[0]-=rb2;
		dorecurse<2>(nd);
		center[2]+=rb2;
		dorecurse<6>(nd);
		center[0]+=rb2;
		dorecurse<7>(nd);
		center[1]-=rb2;
		dorecurse<5>(nd);
		center[0]-=rb2;
		dorecurse<4>(nd);
		
		center[2]-=rb;
		center[0]+=rb;center[1]+=rb;
		++h;rb*=2;

	}




	
};

struct getBallsForCentersAction
{
	//minimum height to use when calculationg the tsets
	//(can be negative)
	int hmincalc;

	//minimum height to go down to determine which centers
	//to use
	int hmincenters;

	vector<TsTree::tset> inner;
	vector<TsTree::tset> border;

	mynum radius;
		

	inline void operator ()(	
		treeWalker<getBallsForCentersAction> & w, 
		const TsTree::tset & t,
		int i,
		const TsTree::ND*n)
	{
					
		sumFloatCircle2 s(w.center,radius,w.tree,hmincalc);
		
		inner.push_back(s.ts);
		border.push_back(s.bordts);
		
		if ( w.h>hmincenters && n )
			w.recurse(n);	
				
	}

		
};

#ifndef BIG_ENDIAN
#define fabsify(x) ((unsigned char*)&x)[sizeof(x)-1] &= ~0x80;
#else
#define fabsify(x) ((unsigned char*)&x)[0] &= ~1;
#endif

struct sumInCircleMultiRadiiAction
{ 
	
	//iterators on tsets for different radii
	TsTree::tset *tsbegin;
	TsTree::mynum *rsqbegin,*rsqend; 
	
	//the center of the balls
	const TsTree::mynum *c;
	
	
	//buffer variables
	TsTree::mynum smallrsq,bigrsq,absc[3], rbp[3],rbm[3],buf,lsq;

	inline void addpointsinside(const TsTree&tree, int leafid)
	{

		vector<TsTree::mynum>::const_iterator 
			it = tree.coords.begin() + tree.off[ leafid ] ,
			itend = tree.coords.begin() + tree.off[ leafid+1] ;	
			

		for(;it!=itend;it+=3)
		{
			
			buf = 	it[0]-c[0]; lsq=buf*buf;
			buf = 	it[1]-c[1]; lsq+=buf*buf;
			buf = 	it[2]-c[2]; lsq+=buf*buf;
			TsTree::tset*ts=tsbegin+(rsqend-rsqbegin);
			TsTree::mynum*rsq=rsqend;
			for( --ts,--rsq; rsq>=rsqbegin && *rsq >= lsq ;--ts,--rsq)				
				ts->add3dcoords(&*it);
		}

	}


	inline void operator ()(	
		treeWalker<sumInCircleMultiRadiiAction> & w, 
		const TsTree::tset & t,
		int i,
		const TsTree::ND*n)
	{
			absc[0]= c[0] -w.center[0];
			absc[1]= c[1] -w.center[1];
			absc[2]= c[2] -w.center[2]; 
			fabsify(absc[0]);
			fabsify(absc[1]);
			fabsify(absc[2]);
			rbp[0]= absc[0]+w.rb;
			rbm[0]= (absc[0]-w.rb) * ( absc[0] > w.rb );
			rbp[1]= absc[1]+w.rb; 
			rbm[1]= (absc[1]-w.rb) * ( absc[1] > w.rb );
			rbp[2]= absc[2]+w.rb; 
			rbm[2]= (absc[2]-w.rb) * ( absc[2] > w.rb );

	
			smallrsq = rbm[0]*rbm[0]+rbm[1]*rbm[1]+rbm[2]*rbm[2];
			bigrsq = rbp[0]*rbp[0]+rbp[1]*rbp[1]+rbp[2]*rbp[2];

		
			TsTree::mynum *rs1,*rs2;

			
			rs1=rsqbegin;
			while( rs1!=rsqend && *rs1 <= smallrsq) ++rs1;
			rs2=rs1;
			while( rs2!=rsqend && *rs2 < bigrsq) ++rs2;			
		
			//for the big radii where the current box is totally inside
			for(TsTree::tset * ts=tsbegin+(rs2-rsqbegin);ts!=tsbegin+(rsqend-rsqbegin) ;++ts)
				*ts += t;


		
			//for the radii where the current box is half inside
			if(rs2>rs1)
			{
				

				TsTree::tset *otsb=tsbegin;
				TsTree::mynum *orsb=rsqbegin,*orse=rsqend;

			
				tsbegin+=rs1-rsqbegin;
				rsqbegin=rs1;
				rsqend=rs2; 

				if(n)
				 	w.recurse(n);
				else
					addpointsinside(w.tree,i);
			
				tsbegin=otsb;
				rsqbegin=orsb;rsqend=orse;

			}
				
	}
};







struct doforcenters
{

private:
	//width of current box
	int w;

	const TsTree&t;

	//minimum width
	int minw;
	//coordinates of lower-left corner of current box,
	int x,y,z;


	//data for sumInCircle
	vector<tset> ts,bordts;
	tset bufts,bufbordts;
	int rad;
	int cent[3];
	int minh;



	void recurse(const TsTree::ND &nd)
	{
		w>>=1;
		int f=nd.f;
		//test if center is inside ball
		for(int i=0;i<8;++i,f>>=1)
		{
			if(f&1)
			{
				x+= (i&1)? w :0;
				y+= (i&2)? w :0;
				z+= (i&4)? w :0;

				if( w>minw )
					recurse(*nd.ch.n[i]);							
				else
				{
					//we use doubled coordinates
					int xx[]={(x<<1)+w,(y<<1)+w,(z<<1)+w};
					t.sumInCircle(bufts,bufbordts,xx,rad<<1,minh);				
					ts.push_back(bufts);
					bordts.push_back(bufbordts);
				}

				x-= (i&1)? w : 0;
				y-= (i&2)? w : 0;
				z-= (i&4)? w : 0;
			}
		}			
		w<<=1;
	}
public:

	doforcenters(const TsTree&tt,unsigned minh_,
		std::vector<tset> &inner,
		std::vector<tset> &border,
		int minhCircle,
		int radius
		)
		:t(tt),minh(minhCircle),rad(radius)
	{
		inner.swap(ts);
		border.swap(bordts);

		w=(1<<t.depth);
		x=y=z=0;
		minw=1<<minh_;
		recurse(t.root);

		inner.swap(ts);
		border.swap(bordts);
	}
};



struct buildFromBuckets
{
private:
	std::vector<int> buckets;
	std::vector<size_t> off;
	//x and y width of buckets grid
	const int nx,ny,nz,nynx,nznynx;
	//width of current box,multiplied by the resp.offsets
	int w,wnx,wnynx;
	//coordinates of lower-left corner of current box,
	//multiplied by their offsets
	int x,ynx,znynx;


	int recurse(TsTree::ND::children &ch)
	{
		w>>=1;wnx>>=1;wnynx>>=1;

		TsTree::ND::children chbuf;

		int f=0;
		if(w>1)
		{
			for(int i=0;i<8;++i){

				x+= w & (-(i&1));
				ynx+= wnx & (-((i&2)>>1));
				znynx+= wnynx & (-((i&4)>>2));

				if((x<nx) &(ynx<nynx) &(znynx<nznynx))
				{
					int bf=recurse(chbuf);		
					if(bf)
					{				
						ND * chni=new ND;
						chni->ch = chbuf;
						chni->f = bf;
						ch.n[i]=chni;
						f|=1<<i;
					}
					else
						ch.n[i]=0;
				}
				else
					ch.n[i]=0;

				x-= w & (-(i&1));
				ynx-= wnx & (-((i&2)>>1));
				znynx-= wnynx & (-((i&4)>>2));

			}
		}
		else {
			for(int i=0;i<8;++i){

				x+= (i&1);
				ynx+= nx & (-((i&2)>>1));
				znynx+= nynx & (-((i&4)>>2));

				if((x<nx) &(ynx<nynx) &(znynx<nznynx))
				{
					int & bid =  buckets[ znynx+ynx+x ];
					if(bid)
					{
						off.push_back(off.back()+3*bid);
						ch.i[i]=bid=off.size()-3;
						f|=1<<i;
					}					
					else
						ch.i[i]=bid=-1;
				}
				else
					ch.i[i]=-1;

				x-= (i&1);
				ynx-= nx & (-((i&2)>>1));
				znynx-= nynx & (-((i&4)>>2));

			}

		}
		w<<=1;wnx<<=1;wnynx<<=1;
		return f;
	}	

public:

	//init t->off, replace numbers in buckets by ids in t->off
	buildFromBuckets( TsTree&t,const int *num,vector<int> &buck)
		:nx(num[0]),ny(num[1]),nz(num[2]),
		nynx(num[0]*num[1]),
		nznynx(num[0]*num[1]*num[2])
	{
		x=ynx=znynx=0;
		t.depth=0;
		int e;
		for(int i=0;i<3;++i)
		{
			frexpl(num[i]-1,&e);
			if(e>t.depth)
				t.depth=e;
		}

		w = 1<<t.depth;
		wnx = nx<<t.depth;
		wnynx = nynx<<t.depth;

		off.swap(t.off);
		buckets.swap(buck);

		t.root.f = recurse(t.root.ch);

		off.swap(t.off);			
		buckets.swap(buck);
	}

};


struct updatetsets
{
private:
	std::vector<TsTree::tset> tsets;

	//current height
	int h;

	void recurse(TsTree::ND & nd)
	{
		--h;

		int f,i;

		if(h>0)
		{
			//first, init tsets of children
			for(i=0,f=nd.f; f!=0 ; ++i,f>>=1)
				if(f&1)
					recurse(*nd.ch.n[i]);


			//then, add them to current tset

			//overjump empty nodes
			for(i=0,f=nd.f;
				(f&1)==0;
				++i,f>>=1);

			//init nd
			nd.t = nd.ch.n[i]->t;

			//add other nd's
			for(++i,f>>=1; f!=0; ++i, f>>=1)
				if(f&1)
					nd.t += nd.ch.n[i]->t;

		}
		else {

			//overjump empty nodes
			for(i=0,f=nd.f;
				(f&1)==0;
				++i,f>>=1);

			//init nd
			nd.t =  tsets[ nd.ch.i[i] ];
			//add other nd's
			for(++i,f>>=1; f!=0; ++i, f>>=1)
				if(f&1)
					nd.t += tsets[ nd.ch.i[i] ];

		}

		++h;
	}	

public:

	//init t->off, replace numbers in buckets by ids in t->off
	updatetsets( TsTree&t )
	{
			//initialize the tsets
		t.tsets.resize(t.off.size()-1);	

		//set tsets to 0
		memset(&t.tsets[0],0,t.tsets.size()*sizeof(t.tsets[0]));

		for(int id=0;id<t.off.size()-1;++id)
		{
			TsTree::tset &ts=t.tsets[id];
			vector<long double>::iterator 
				it    =  t.coords.begin() + t.off[id],
				itend	=  t.coords.begin() + t.off[id+1];
			for(;it!=itend;it+=3)
				ts.add3dcoords(&*it);
		}

		h=t.depth;
		tsets.swap(t.tsets);
		recurse(t.root);
		tsets.swap(t.tsets);
	}

};



struct buildFromPoints
{

	typedef unsigned long int luint;
  

	  
  typedef long double mynum;

  struct pt
  {
    luint index;	
    const mynum * p;
    bool operator <( const pt&x ) const
    {
      return index<x.index;
    }
  };

	/**
	*build tree from points in pts
	*\param boxwidth: minimum boxwidth
	*\param t: the tree to be initialized;
	* max,min should be already built.
	*\pre: max,min of TsTree are initialized
	*/
  buildFromPoints( TsTree&t,const vector<mynum> &pts,mynum boxwidth)
  {
    vector<pt> p(pts.size()/3);

		static const int maxdepth=sizeof(luint)*8/3;
		
		{
			//calculate leaf width and height of tree
			
			//calculate max edge length of bb.
			mynum wmax=t.max[0]-t.min[0];
			mynum v=t.max[1]-t.min[1];
			if(wmax<v)wmax=v;
			v=t.max[2]-t.min[2];
			if(wmax<v)wmax=v;


			boxwidth *= sqrt(t.var[0]+t.var[1]+t.var[2]);


			int pboxwidth,pwmax;
		  frexpl(boxwidth,&pboxwidth);
			frexpl(wmax,&pwmax);
			if(pwmax > maxdepth +pboxwidth|| boxwidth==0 )
				pboxwidth = pwmax-maxdepth;

			t.rb0=ldexpl(1.,pboxwidth);
			t.depth=pwmax-pboxwidth;
		}
    mynum invrb0=1.0/t.rb0;
    
    
    vector<pt>::iterator i;
		vector<mynum>::const_iterator ipts = pts.begin();

    for(i=p.begin();i!=p.end();++i,ipts+=3)
    {
    	i->p= &*ipts; 
	i->index=zord<luint>(
		(luint)( (ipts[0]-t.min[0])*invrb0), 
		(luint)( (ipts[1]-t.min[1])*invrb0), 
		(luint)( (ipts[2]-t.min[2])*invrb0) 
	);    	   
    }
    
    sort(p.begin(),p.end());
    vector<size_t> off(p.size()+1);     
    vector<size_t>::iterator offit=off.begin();
    *offit=0;++offit;
    
    //number of nodes that will be in the octree
    //in addition to the root node;
    size_t nodenumber= t.depth-1;
    if(p.size()>=1) 
		{
			i=p.begin();
			do {
				++i; 
				while( i!=p.end() && i->index == (i-1)->index )
					++i;

				if(i!=p.end()){
					luint delt= i->index ^ (i-1)->index;
					for(delt>>=3;delt;delt>>=3)
						++nodenumber;
				} 

				*offit=i-p.begin();++offit;       	
			} while( i!=p.end() );
		}
		off.resize(offit-off.begin());
    
    
    vector<mynum> coords(pts.size());
    vector<mynum>::iterator it=coords.begin();
    for(i=p.begin();i!=p.end();++i)
    {
    	*it=i->p[0];++it;
    	*it=i->p[1];++it;
    	*it=i->p[2];++it; 
    }
    
    t.root.f=0;
    
    //now build tree with z-index and offsets
    t.memfornodes.resize(nodenumber);
    vector<TsTree::ND>::iterator memptr=t.memfornodes.begin();
		vector<TsTree::ND*> st(t.depth);
		st.back() = &t.root;

		size_t currleafid;
		luint oindex=(luint)-1; 
    for(offit=off.begin();offit!=off.end()-1;++offit)
    {
    	luint currindex=p[*offit].index;
			
			int h=0;
			if(offit==off.begin())
				h=t.depth;
			else
			{
				luint deltindex = (oindex^currindex);
				while(deltindex)
					++h,deltindex>>=3 ;
			}

			int h3=h*3;
			for(--h,h3-=3 ; h>0; --h,h3-=3)				
			{
				int childid=(currindex>>h3)&7;
				TsTree::ND&nd = *st[h];
				( st[h-1] = nd.ch.n[childid] = &*(memptr++) )->f = 0 ;				
				nd.f |= 1<<childid;
			}
			int childid=currindex&7;
			TsTree::ND&nd = *st[0];
			nd.ch.i[childid] = offit-off.begin();
			nd.f|=1<<childid;

			oindex=currindex;
      
			//off should be an index in the coords array,
			//so multiply by three
			*offit *= 3;
    }	

		*offit*=3;
    
		//return arrays to tree
    t.coords.swap(coords);
    t.off=off;
    
    updatetsets updts(t);        
  }
  

};


static void buildFromPoints_old(TsTree&t, const std::vector<long double> & pts,long double boxwidth)
{

	t.rb0 = sqrt(t.var[0]+t.var[1]+t.var[2]) * boxwidth;

	int num[3]={2,2,2};

	if(t.rb0==0) t.rb0=1e-6;
	int expr;
	frexpl(t.rb0,&expr);

	t.rb0=ldexpl(1,expr);

	//maximum memory for the grid(512Mbytes).
	size_t maxmem=512<<20;	


	long double invrb0;

	for(;; t.rb0*=2 ){

		invrb0=1.0l/t.rb0;

		//number of box grid in every dim
		for(int i=0;i<3;++i){ 
			num[i]=int( floor( (t.max[i]-t.min[i])*invrb0 ) )+2;
			if(num[i]<2)num[i]=2; 
		}

		if(maxmem/sizeof(int)/num[0]/num[1]/num[2]>0)
			break;

	}

	//ids in the buckets array
	vector<int> ids(pts.size()/3);
	vector<int>::iterator iit;

	//a rectangular grid containing the buckets
	vector<int> buckets(num[0]*num[1]*num[2],0);

	//number of filled buckets
	int nfilled=0;
	
	vector<long double>::const_iterator it;
	for(it=pts.begin(),iit=ids.begin();
		iit!=ids.end();++iit,it+=3)
	{
		int id=(int) floor( (it[2] -t.min[2] )*invrb0 );
		id*=num[1] ;
		id+=(int)floor( (it[1] -t.min[1] )*invrb0 );
		id*=num[0];
		id+=(int)floor( (it[0] -t.min[0] )*invrb0 );

		if(!buckets[id])
			++nfilled;

		++buckets[id];
		*iit = id;
	}

	//cout<<nfilled<<" of "<<buckets.size()<<" buckets filled ";

	//the offsets in the coordinate array
	t.off.reserve( nfilled + 2 );
	t.off.push_back(0);//constant starting offset
	t.off.push_back(0);//offset of current point in 1st bucket


	//build the tree datastructure
	//(without initializin the tensorsets)
	//and build the off array
	impl::buildFromBuckets bfb(t,num,buckets);
	//cout<<"depth is "<<depth<<endl;


	t.coords.resize(pts.size());

	int iii=0;
	for(it=pts.begin(),iit=ids.begin();
		iit!=ids.end();
		++iit)
	{
		++iii;
		int id=buckets[*iit];
		if(id<0)continue;

		vector<long double>::iterator 
			cit= t.coords.begin()+t.off[id+1];

		//insert coords into list
		*cit=*it; ++cit,++it;
		*cit=*it; ++cit,++it;
		*cit=*it; ++cit,++it;

		//add 3 to off
		t.off[id+1]+=3;

	}



	//update the tensor sets in the nodes
	updatetsets upd(t);
};

};//end of TsTree::impl



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
#include "TsTree2.h"
#include "zindex.h"
using namespace std;




struct TsTree2::impl{


	//update the tset sums from vector tsets
	struct updatetsets;


//help structures for sumInCircle
	struct sumCircle2;
	struct sumFloatCircle2;


	struct doforcenters;


	




/**optimized version of sumcircle*/
struct sumCircle2
{
	const TsTree2&tree;
	//the tensor set where all others are added to
	tset ts;

	//the sum of the tensorsets intersected by the ball border
	tset bordts;

	//the x,y,z of the center of the current node
	//relative to the center of the ball
	int x,y,z;

	//the real values of the Ball
	TsTree2::mynum center_[3];
	TsTree2::mynum rsq_;

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



	inline void addpointsinside(TsTree2::Leaf*l)
	{

		const TsTree2::mynum 
		    *it = l->begin ,
				*itend = l->end ;	

		mynum buf,lsq;			

		ptssummed += (itend-it)/3;

		for(;it!=itend;it+=3)
		{
			lsq=rsq_;
			buf = 	it[0]-center_[0]; lsq-=buf*buf;
			buf = 	it[1]-center_[1]; lsq-=buf*buf;
			buf = 	it[2]-center_[2]; lsq-=buf*buf;
			if(lsq>0)
				ts.add3dcoords(it);
		}

	}


	/*calculate ts as the sum of all leaves intersecting			//the ball
	* given by |x - 0.5*center | < 0.5*radius
	*/
	sumCircle2(const int *center, int radius, const TsTree2 & tr,
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
			ts = tr.root->t;
			return;
		}

		if(outside())			
			return;

		if(h==hmin)
			bordts = tr.root->t;
		else if(tr.root->isleaf)
			addpointsinside((Leaf*)tr.root);
		else
			recurse((Node*)tr.root);		
	}









	inline void dorecurse(bool ins,bool out,int i,const TsTree2::Node*n)
	{
	        const TsTree2::ND * c = n->ch[i];
	        
	        if( !c )
                    return;
	        if(out && outside()) //totally outside
                    return;

                if(ins && inside() ) //totally inside						
                {
                    ++numbersummed;					
                    ts += c->t;	
                }
                else if(c->isleaf) 
									addpointsinside((TsTree2::Leaf*)c);
                else if(h>hmin) 
                  recurse((TsTree2::Node*)c);
                else 
                  bordts+=c->t;
                
	}

	void recurse(const TsTree2::Node * nd)
	{
	
		bool couldbeinside,couldbeoutside;


		--h;
		int rb2=rb;
		
		rb>>=1;
		//test if center is inside ball
		int d=x*x+y*y+z*z;
		couldbeinside  = d <= rsq;
		couldbeoutside = d >= rsq;

		
		x-=rb;y-=rb;z-=rb;

			dorecurse(couldbeinside,couldbeoutside,0,nd);

		x+=rb2;

			dorecurse(couldbeinside,couldbeoutside,1,nd);

		y+=rb2;

			dorecurse(couldbeinside,couldbeoutside,3,nd);

		x-=rb2;

			dorecurse(couldbeinside,couldbeoutside,2,nd);

		z+=rb2;

			dorecurse(couldbeinside,couldbeoutside,6,nd);

		x+=rb2;

			dorecurse(couldbeinside,couldbeoutside,7,nd);

		y-=rb2;

			dorecurse(couldbeinside,couldbeoutside,5,nd);

		x-=rb2;

			dorecurse(couldbeinside,couldbeoutside,4,nd);

		z-=rb;
		x+=rb;y+=rb;

		++h;rb<<=1;

	}




	




};


/** optimized version of sumFloatCircle*/
struct sumFloatCircle2
{
	typedef TsTree2::mynum mynum;
	const TsTree2&tree;
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



	inline void addpointsinside(const TsTree2::Leaf*l)
	{

		const mynum  
			*it = l->begin ,
			*itend = l->end ;	


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
	sumFloatCircle2(const mynum center[3], mynum radius, const TsTree2 & tr,
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
			ts = tr.root->t;
			return;
		}

		if(outside())			
			return;

		if(h==hmin)
			bordts = tr.root->t;
		else if(tr.root->isleaf)
			addpointsinside((Leaf*)tr.root);
		else
			recurse((Node*)tr.root);		
	}

	inline void dorecurse(bool ins,bool out,int i,const TsTree2::Node *n)
	{

		const TsTree2::ND* c=n->ch[i];
		if(!c)return;

		if(out && outside())return;

		if(ins && inside() ) //totally inside						
		{
			++numbersummed;					
			ts+=c->t;
		}
		else if(c->isleaf)
			addpointsinside((Leaf*)c);
		else if(h>hmin)
			recurse((Node*)c);
		else
			bordts+=c->t;
	}

	void recurse(const TsTree2::Node * nd)
	{
	
		bool couldbeinside,couldbeoutside;


		--h;
		mynum rb2=rb;
		
		rb*=0.5;
		//test if center of box is inside ball
		mynum d=x*x+y*y+z*z;
		couldbeinside  = d <= rsq;
		couldbeoutside = d >= rsq;

		
		x-=rb;y-=rb;z-=rb;

			dorecurse(couldbeinside,couldbeoutside,0,nd);

		x+=rb2;

			dorecurse(couldbeinside,couldbeoutside,1,nd);

		y+=rb2;

			dorecurse(couldbeinside,couldbeoutside,3,nd);

		x-=rb2;

			dorecurse(couldbeinside,couldbeoutside,2,nd);

		z+=rb2;

			dorecurse(couldbeinside,couldbeoutside,6,nd);

		x+=rb2;

			dorecurse(couldbeinside,couldbeoutside,7,nd);

		y-=rb2;

			dorecurse(couldbeinside,couldbeoutside,5,nd);

		x-=rb2;

			dorecurse(couldbeinside,couldbeoutside,4,nd);

		z-=rb;
		x+=rb;y+=rb;

		++h;rb*=2;

	}




	
};



struct treeWalkerData
{
	treeWalkerData(const TsTree2&tr)
		:tree(tr)
	{}

	typedef TsTree2::mynum mynum;
	const TsTree2&tree;

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

	treeWalker(const TsTree2 & tr)
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

		if(!tree.root->isleaf){
			doit( *this, tree.root->t, 0, (Node*)tree.root);
		}
		else{
			doit( *this, tree.root->t, (Leaf*)tree.root , 0 );
		}
	}

	template<int i>
	inline void dorecurse(const Node*nd)
	{
	  const ND * c=nd->ch[i];
		if(c)
		{
			if(!c->isleaf){
				doit( *this, c->t, 0, (Node*)c);
			}
			else{
				doit( *this, c->t, (Leaf*)c , 0 );
			}
		}
	}

	void recurse(const Node * nd) 
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

	vector<TsTree2::tset> inner;
	vector<TsTree2::tset> border;

	mynum radius;
		

	inline void operator ()(	
		treeWalker<getBallsForCentersAction> & w, 
		const TsTree2::tset & t,
		const TsTree2::Leaf*i,
		const TsTree2::Node*n)
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
	TsTree2::tset *tsbegin;
	TsTree2::mynum *rsqbegin,*rsqend; 
	
	//the center of the balls
	const TsTree2::mynum *c;
	
	
	//buffer variables
	TsTree2::mynum smallrsq,bigrsq,absc[3], rbp[3],rbm[3],buf,lsq;

	inline void addpointsinside(const TsTree2&tree, const TsTree2::Leaf*l)
	{

		const TsTree2::mynum
			*it = l->begin ,
			*itend = l->end ;	
			

		for(;it!=itend;it+=3)
		{
			
			buf = 	it[0]-c[0]; lsq=buf*buf;
			buf = 	it[1]-c[1]; lsq+=buf*buf;
			buf = 	it[2]-c[2]; lsq+=buf*buf;
			TsTree2::tset*ts=tsbegin+(rsqend-rsqbegin);
			TsTree2::mynum*rsq=rsqend;
			for( --ts,--rsq; rsq>=rsqbegin && *rsq >= lsq ;--ts,--rsq)				
				ts->add3dcoords(it);
		}

	}


	inline void operator ()(	
		treeWalker<sumInCircleMultiRadiiAction> & w, 
		const TsTree2::tset & t,
		const TsTree2::Leaf*l,
		const TsTree2::Node*n)
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

		
			TsTree2::mynum *rs1,*rs2;

			
			rs1=rsqbegin;
			while( rs1!=rsqend && *rs1 <= smallrsq) ++rs1;
			rs2=rs1;
			while( rs2!=rsqend && *rs2 < bigrsq) ++rs2;			
		
			//for the big radii where the current box is totally inside
			for(TsTree2::tset * ts=tsbegin+(rs2-rsqbegin);ts!=tsbegin+(rsqend-rsqbegin) ;++ts)
				*ts += t;


		
			//for the radii where the current box is half inside
			if(rs2>rs1)
			{
				

				TsTree2::tset *otsb=tsbegin;
				TsTree2::mynum *orsb=rsqbegin,*orse=rsqend;

			
				tsbegin+=rs1-rsqbegin;
				rsqbegin=rs1;
				rsqend=rs2; 

				if(n)
				 	w.recurse(n);
				else
					addpointsinside(w.tree,l);
			
				tsbegin=otsb;
				rsqbegin=orsb;rsqend=orse;

			}
				
	}
};

/*
struct doforcenters
{

private:
	//width of current box
	int w;

	const TsTree2&t;

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



	void recurse(const TsTree2::ND &nd)
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

	doforcenters(const TsTree2&tt,unsigned minh_,
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

*/

struct updatetsets
{
private:
	std::vector<Leaf> leaves;


	void recurse(Node * n)
	{

		int i;

		//first, init tsets of children
		for(i=0;i<8; ++i)
			if(n->ch[i] && ! n->ch[i]->isleaf )
				recurse( (Node*)n->ch[i]);


		//then, add them to current tset

		//overjump empty nodes
		i=0;
		while(i<8 && !n->ch[i] ) ++i;

		//init nd
		n->t = n->ch[i]->t;

		//add other nd's
		for(++i; i<8; ++i )
			if(n->ch[i])
				n->t += n->ch[i]->t;		

	}	

public:

	//init t->off, replace numbers in buckets by ids in t->off
	updatetsets( TsTree2&t )
	{
	  leaves.swap(t.leaves);
	        

		for(int id=0;id<leaves.size();++id)
		{
			TsTree2::tset &ts=leaves[id].t;
			ts.reset();
			mynum
				*it    =  leaves[id].begin,
				*itend	=  leaves[id].end;
			for(;it!=itend;it+=3)
				ts.add3dcoords(&*it);
		}

		if(!t.root->isleaf)
			recurse((Node*)t.root);
		
		leaves.swap(t.leaves);
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
	*\pre: max,min of TsTree2 are initialized
	*/
  buildFromPoints( TsTree2&t,const vector<mynum> &pts,mynum boxwidth,int ptperleaf)
  {
    vector<pt> p(pts.size()/3);

		static const int maxdepth=sizeof(luint)*8/3;

		//maximum number of points per leaf
		int maxptperleaf=ptperleaf;
		if(ptperleaf==0)maxptperleaf=32;
		
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
			if(pwmax-pboxwidth > maxdepth||boxwidth==0)
				pboxwidth = pwmax-maxdepth;
			if(pwmax<pboxwidth)pboxwidth=pwmax;

			t.rb0=ldexpl(1.,pboxwidth);
			t.depth=pwmax-pboxwidth;
		}
    mynum invrb0=1.0/t.rb0;


		if(t.depth==0)
		{
			t.leaves.resize(1);
			t.root=&t.leaves.front();
			t.coords=pts;
			t.root->isleaf=true;
			t.leaves.front().begin=&t.coords[0];
			t.leaves.front().end=&t.coords[0]+t.coords.size();
			updatetsets upd(t);
			return;
		}
    
    
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

    vector<mynum> coords(pts.size());
    vector<mynum>::iterator cit=coords.begin();
    for(i=p.begin();i!=p.end();++i)
    {
    	*cit=i->p[0];++cit;
    	*cit=i->p[1];++cit;
    	*cit=i->p[2];++cit; 
    }
    


		//calculate positions of msb of differences in zindex 
		//divided by 3 (equivalent to tree heights where differences apply)
		vector<unsigned char> heights(p.size()+1);		
		heights.front()=t.depth+1;
    if(p.size()>=1) 
		{
			vector<unsigned char>::iterator ht=heights.begin()+1;
			for(i=p.begin()+1;i!=p.end();++i,++ht){
				int h=0;
				for(luint delt= i->index ^ (i-1)->index; delt ; delt>>=3)	
					++h;

				*ht = (unsigned char)h;
			} 
		}
		heights.back() = t.depth+1;


		//number of nodes that will be in the octree
    size_t nodenumber= 0;    
		
		//the indices in the coords array where a leaf starts or ends.
		vector<size_t> off(p.size());     
    
		vector<size_t>::iterator offit=off.begin();		
    *offit=0;++offit;
		
		vector<unsigned char>::iterator ht=heights.begin();
		while(ht < heights.end()-1)
		{
			
			unsigned char maxht=0,currht=*ht;
			int nexti=0;

			/*
			by maxht<=currht it is also prevented that the index 
			runs over the border, because the last height is tree.depth.
			if maxht is still zero while i>maxptperleaf, there are more points 
			than the desired number that have to go into the same leaf. 
			It also means that the tree's depth has been calculated too low,
			so the zindex is not fine enough
			*/
			for(int i=1; (i<=maxptperleaf || maxht==0) && maxht<currht ;++i)
			{
				if(ht[i]>maxht)
					maxht=ht[i],nexti=i;		
			}
			
			ht+=nexti;
			
			*offit = ht-heights.begin(); ++offit;

			if(currht > maxht) nodenumber += currht-maxht;
		}

		off.resize(offit-off.begin());
        
    //now build tree with z-index and offsets
    vector<Node> nodes(nodenumber);
		if(nodes.size())
			memset(&nodes.front(),0,sizeof(nodes.front())*nodes.size());
		vector<Leaf> leaves(off.size()-1);

		//cout<<"no.leaves:"<<leaves.size()<<endl;
		//cout<<"no.nodes:"<<nodes.size()<<endl;

    vector<Node>::iterator nodeit=nodes.begin();
		
		//some type of stack
		vector<ND*> st(t.depth+1);

		size_t currleafid;
		
		//initialize leaves &nodes
		vector<Leaf>::iterator lit=leaves.begin();
    for(offit=off.begin();offit!=off.end()-1;++offit,++lit)
    {
			lit->isleaf=true;
			lit->begin  = &coords[0] + offit[0] *3;
			lit->end = &coords[0] + offit[1] *3; 

    	luint currindex=p[*offit].index;

			int 
				h0 = heights[offit[0]],
				h1 = heights[offit[1]];

			if(h1>h0)h1=h0;

			int h = h1;			
			currindex >>= 3*(h-1);
			st[h-1] = &*lit;
			while(h<h0){ //create nodes												
				st[h] = &*nodeit;					
				nodeit->ch[currindex&7] =st[h-1];
				++h;++nodeit;currindex>>=3;
			}
			//if last created node wasn't the root node, insert it 
			//into the next higher node
			if(h<st.size())
				((Node*)st[h])->ch[currindex&7] = st[h-1];			
    }	

		//return arrays to tree
    t.coords.swap(coords);
		t.nodes.swap(nodes);
		t.leaves.swap(leaves);
		t.root = st.back();
    
		//fill tensorsets with the correct values.
    updatetsets updts(t);        
  }
  

};



};//end of TsTree2::impl



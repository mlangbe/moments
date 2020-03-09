/*
	implementations of kdtree used to optimize the search in moment invraint sets

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

#include <math.h>
#include<time.h>
#include<algorithm>
#include<iostream>
#include<cassert>
#include<fstream>

#ifndef GNUCC
#pragma once
#endif


#ifndef OLD_VC
#define _TYPENAME_ typename
#else
#define _TYPENAME_
#endif

#include<vector>

template<class mynum>
class Momkdtree_templ
{
public:
	//the number type
	struct node{
		//the pivot value
		mynum p;
		//the coordinate id in which to split
		int i;

		node(){}
		node(const mynum& pp,int ii):i(ii),p(pp){}
	};

private:

	int dim; //dimension of the stored points
	int height; //height of the tree

	std::vector<mynum> bbox;

	//the original leaves in original order
	std::vector<mynum> leaves;	

	//the ids of the leaves in tree order
	//for every leaf one or two entries
	std::vector<int> leafids;

	//the pivot values
	std::vector<node> nodes;
	//the coordinate ids for the pivots
//	std::vector<int> pivotcids;

public:

	//namespace for implementation structs
	struct impl;

	int getDim() const {return dim;}
	int getHeight() const {return height;}

	const std::vector<mynum> &getleaves() const
	{return leaves;}

	//swap the leaves with the numbers in l
	void swapleaves(std::vector<mynum> & l)
	{leaves.swap(l);}

	Momkdtree_templ(void);

	//build kdtree from dim-dimensional points stored in mynum
	void init(std::vector<mynum> &moments,int dim);
	
	Momkdtree_templ(std::vector<mynum> &moments,int dim)
	{	init(moments,dim);}


	/**
	*get all leaves in the box.
	*\param moms: moments returned
	*\param box: an array with 2*dim entries
	*which are min,max-pairs for every dimension
	*/
	void momsInBox(std::vector<mynum> &moms,const mynum*box) const ;

	/**
	*get all leaves in the box.
	*\param momids: start index of moment set in moment set array
	*\param box: an array with 2*dim entries
	*which are min,max-pairs for every dimension
	*/
	void momIdsInBox(std::vector<int> &momids,const mynum*box) const;
	
	/**
	*get all leaves in the Ball.
	*\param momids: start index of moment set in moment set array
	*\param center: an array with dim entries containing 
	*the center of the hyperball
	*\param radius: the radius of the hyperball
	*/
	void momIdsInBall(std::vector<int> &momids,
		const mynum*center,mynum radius) const;

	~Momkdtree_templ(void);
};

using namespace std;

template<class mynum>
struct Momkdtree_templ<mynum>::impl{


	//recursion structure for building the kdtree
	struct build
	{
		struct vi
		{
			//the value in the current direction
			//of moms[i]
			long double v;
			//the index in moms
			int i;

			bool operator < (const vi&o) const
			{ return v<o.v; }

		};
		
		//the bbox of the current leaf	
		vector<mynum> bbox;
		
		//a buffer var for getmaxdimbb;
		vector<mynum> bbbuf;
		
		//the moments
		vector<mynum> moms;

		vector<node> nodes;
		
		typename vector<node>::iterator nodeit;
		vector<int> leafids;
		vector<int>::iterator leafit;
		
		int dim;//the dimension

		int h; //the current height
		
		//ids that will be re-sorted, along with their values
		vector< vi > ids;




		inline void initbb(typename vector<mynum>::iterator m,
												 typename vector<mynum>::iterator b)
		{
_TYPENAME_			vector<mynum>::iterator mend=m+dim;
			while(m<mend)
			{
				*b=*m; ++b; *b=*m; ++b,++m;
			}
		}
		
		inline void updatebb(typename vector<mynum>::iterator m,
													 typename vector<mynum>::iterator b)
		{
			#ifndef OLD_VC
			typename
			#endif
			vector<mynum>::iterator mend=m+dim;
			while(m<mend)
			{
				if(*b>*m)*b=*m;
				++b;
				if(*b<*m)*b=*m;
				++m,++b;
			}
		}

		//get the coordinate in which the bbox is greatest
		inline int getmaxdimbb(
			typename vector<vi>::iterator l,typename vector<vi>::iterator r
			)
		{		



			if(r-l<20)
			{
				//if only asmall number, ompute from all
				//points,
				initbb( moms.begin() + l->i,bbbuf.begin() ) ;				
				
				for(++l;l!=r;++l)
					updatebb(moms.begin() + l->i,bbbuf.begin() );			
			}
			else
			{
				//otherwise, compute from a random set of pts
				int s;
				s= int(float(r-l-1) * rand() / RAND_MAX);
				initbb( moms.begin() + l[s].i,bbbuf.begin() ) ;				
				
				for(int i=1;i<10;++i)
				{
					s= int(float(r-l-1) * rand() / RAND_MAX);
					updatebb(moms.begin() +l[s].i,
					bbbuf.begin() );			
				}
			}

			mynum maxd=bbbuf[1]-bbbuf[0];
			int maxi=0;
			for(int i=2;i<dim*2;i+=2)
			{
				mynum d=bbbuf[i+1]-bbbuf[i];
				if(d>maxd)
					maxd=d,maxi=i;
			}			
			return maxi/2;		

		}



		inline void initbbox()
		{
			//calculate bbox
			#ifndef OLD_VC
			typename
			#endif
			vector<mynum>::iterator it,bbit;

			//init bbox to the 1st laues
			initbb(moms.begin(),bbox.begin());

			for(it=moms.begin()+dim;it<moms.end();it+=dim)
			{
				updatebb(it,bbox.begin());
			}
		
		}

		build(Momkdtree_templ & t, vector<mynum> &moments )
			:dim(t.dim),bbox(t.dim*2),bbbuf(t.dim*2)
		{
			moms.swap(moments);		
			frexpl(moms.size()/dim -1 , &h);
			t.height=h;


			leafids.resize(1<<h);
			leafit=leafids.begin();

			nodes.resize((1<<h)-1);
			nodeit=nodes.begin();

			//calculate bbox from moms
			initbbox();

			//...............................
			//initialize ids
			ids.resize(moms.size()/dim);
			for(int i=0,j=0;i<ids.size();++i,j+=dim)
				ids[i].i=j;

			recurse(ids.begin(),ids.end());


			t.nodes.swap(nodes);
			t.leaves.swap(moms);
			t.leafids.swap(leafids);
			t.bbox.swap(bbox);
		
		}


		void recurse(typename vector<vi>::iterator l,typename vector<vi>::iterator r)
		{

			if(r-l<=1)
			{
				//fill pivots cids up so you still have a
				//complete binary tree in list repr.
				
				*(leafit++)=l->i;
	
				if(h>0)
				{
					*(leafit++)=l->i;
					*(nodeit++)=node(l->v,-1);
				}
						
				return;
			}

			//get coordinate id of the max.extension of the bb
			int ci=getmaxdimbb(l,r);
			
			//fill range with the values of the current coordinate:
			#ifndef OLD_VC
			typename
			#endif
			vector<vi>::iterator vit;
			for(vit=l;vit!=r;++vit)
				vit->v = moms[vit->i+ci];

			//split range at its median
			vit = l+(r-l)/2;
			nth_element(l, vit	,r);

			*(nodeit++) = node(vit->v,ci);



			mynum &bbmin=bbox[2*ci], &bbmax=bbox[2*ci+1];

			mynum obb;

			--h;			

			//do recursion

			obb=bbmax; 
			bbmax = vit->v;
			recurse(l,vit);
			bbmax = obb; 
			
			obb=bbmin;
			bbmin = vit->v;
			recurse(vit,r);
			bbmin = obb;

			++h;
		}
	

	};


	struct getbox{
		//statistic: for getBox: how many thrown out at leaf level
		static int numnodesvisited;
	
	
		int dim;
		//leaves in current subtree
		int lst;
 		vector<int>::const_iterator leafit;
		typename vector<node>::const_iterator nodeit;
	
		//returned leaf ids
		vector<int> ret;
		
		const Momkdtree_templ&t;
		const mynum*box;

		getbox(vector<int>&rett,const Momkdtree_templ&tt, const mynum boxx[])
			:dim(tt.dim),lst(1<<tt.height),t(tt),box(boxx)
		{
				rett.swap(ret);

				ret.clear();
				//test if box intersects t at all

				for(int i=0;i<2*dim;i+=2)
				{
					if( t.bbox[i] > box[i+1] ) return;
					if( t.bbox[i+1] < box[i] ) return;				
				}

				leafit = t.leafids.begin();
				nodeit = t.nodes.begin();
				
				recurse();			
				rett.swap(ret);
		}

		//test if the current leaf is inside the box
		//and add it if it is
		inline void addleaf()
		{
			#ifndef OLD_VC
			typename
			#endif

			vector<mynum>::const_iterator 
			it = t.leaves.begin() + *leafit,itend=it+dim;
						const mynum*b;
						
						for(b=box;it!=itend;++it,b+=2)
							if(*it < b[0] || *it > b[1] )
								break;

						if(it==itend)
							ret.push_back(*leafit);	
			
						++numnodesvisited;
		}

		void recurse()
		{
			++numnodesvisited;
				int cas=0;
				
				{
					int cid= nodeit->i;
					if(cid>=0)
					{
						mynum p=nodeit->p;
					
						cas = 
							/*does box intersect the left subtree?*/
							box[cid*2]<=p 
							/*does box intersect the right subtree?*/
							| (box[cid*2+1]>=p)<<1; 					
					}
					else
						//if this node has the same child twice:
						cas=2;
				}
				
				++nodeit;
			
				lst>>=1;
				
				if(lst>1)
				{
					if(cas&1)
						recurse();
					else
						leafit+=lst, nodeit+=lst-1;

					if(cas&2)
						recurse();
					else
						leafit+=lst, nodeit+=lst-1;
				}					
				else
				{
					if(cas&1)
						addleaf();

					++leafit;

					if(cas&2)
						addleaf();					

					++leafit;
				}

					
				lst<<=1;

		}		
	};

	struct getball{
		//statistic: for getBox: how many thrown out at leaf level
		static int numnodesvisited;
	
	
		int dim;
		//leaves in current subtree
		int lst;
		typename vector<int>::const_iterator leafit;
		typename vector<node>::const_iterator nodeit;
	
		//returned leaf ids
		vector<int> ret;
		
		const Momkdtree_templ&t;

		const mynum*ctr;
		
		//the current bounding box
		mynum *bb;		
		
		//the radius, the radius squared
		mynum r,rsq;


		mynum distinside,distoutside;

		mynum *disti, *disto;
		

		inline bool boxInsideBall() const
		{
			mynum delt=0,x,c,d;
			const mynum*b = bb;
			for(int i=0;i<dim;++i,b+=2)
			{
				c = ctr[i] - 0.5*(b[0]+b[1]);
				d = 0.5*(b[1]-b[0]);
				x = fabs( (mynum)c ) + d;
				delt+=x*x;				
			}			
			return delt<=rsq;
		}

		inline bool boxOutsideBall() const
		{
			mynum delt=0,x,c,d;
			const mynum*b = bb;
			for(int i=0;i<dim;++i,b+=2)
			{
				d = 0.5*(b[1]-b[0]);
				c = ctr[i] - (b[0]+d);
				x=fabs( (mynum)c);
				if(x>d)
				{
					x -= d;
					delt+=x*x;				
				}
			}			
			return delt>=rsq;
	
		}

		getball(vector<int>&rett,const Momkdtree_templ&tt, 
			const mynum*center,mynum radius)
			:
		dim(tt.dim),lst(1<<tt.height),t(tt),r(radius),
		rsq(radius*radius),ctr(center)
		{
				rett.swap(ret);

				bb=new mynum[dim*2];
				std::copy(tt.bbox.begin(),tt.bbox.end(),bb);

				disto=new mynum[dim];
				disti=new mynum[dim];
				
				
				distinside=0;distoutside=0;
				mynum x,c,d;

				const mynum*b = bb;
				for(int i=0;i<dim;++i,b+=2)
				{
					mynum d = .5 *(b[1] - b[0]);

					//the distance of the ball center to the box center
					mynum c = fabs((double)(center[i] - (b[0] + d)));

					x = c+d ;
					disti[i]=x*x;
					distinside+=disti[i];
	
					x = c;
					if(x>d)
					{
						x -= d;
						disto[i]=x*x;
						distoutside+=disto[i];				
					}	
					else
						disto[i]=0;
				}

				if(distoutside>=rsq)
					return;



				ret.clear();

				leafit = t.leafids.begin();
				nodeit = t.nodes.begin();

				if(distinside<rsq)
					addleaves();
				
				recurse();			
				rett.swap(ret);

				delete[]bb;
				delete[]disto;
				delete[]disti;
		}


		//if current leaf is inside,add it
		inline void addleaf()
		{
			#ifndef OLD_VC
			typename
			#endif
			vector<mynum>::const_iterator 
				it = t.leaves.begin()+*leafit;

			mynum sq=0,x;
			for(int i=0;i<dim;++i)
				x=it[i]-ctr[i],sq+=x*x;

			if(sq<=rsq)
				ret.push_back(*leafit);			

			++numnodesvisited;
		}

		//add all leaves inside the current box
		inline void addleaves()
		{
			vector<int>::const_iterator it=leafit,itend=leafit+lst;
			for(;it!=itend;++it){
				if(it+1 !=itend && it[0]==it[1] )++it;
				ret.push_back(*it);
			}
			
			++numnodesvisited;
		}

		inline void doforsubtree(int cid)
		{
			mynum odi=distinside,odo=distoutside;
			{
					const mynum*b= bb+(cid<<1);
					mynum d = .5 *(b[1] - b[0]);

					//the distance of the ball center to the box center
					mynum c = fabs((double) (ctr[cid] - (b[0] + d)));


					disti[cid] = c+d;
					disti[cid] *=disti[cid];
					
					if(c>d)
						disto[cid]=c-d,
						disto[cid] *=disto[cid];
					else
						disto[cid]=0;
			}
			distinside+=disti[cid];
			distoutside+=disto[cid];				

			assert( boxOutsideBall()==(distoutside>=rsq) );
			assert( boxInsideBall()==(distinside<=rsq) );

			if(distoutside<rsq){	//of box not totally outside ball			
				if(distinside<=rsq) //if box totally inside ball
				{
					addleaves();
					leafit+=lst, nodeit+=lst-1;
				}
				else //box intersects ball somehow
					recurse();
			}
			else //box totally outside ball
			{
				leafit+=lst, nodeit+=lst-1;
			}

			distinside=odi;distoutside=odo;
		}

		void recurse()
		{						
			++numnodesvisited;
				int cid = nodeit->i;
				mynum p = nodeit->p;
				
				++nodeit;
			
				lst>>=1;

				if(lst>1)
				{
					mynum buf,*bbb=bb+(cid<<1),
						odisti=disti[cid],odisto=disto[cid];

					distinside-=odisti;
					distoutside-=odisto;

					//do everything for first subtree
					buf = bbb[1];bbb[1]=p;
					doforsubtree(cid);
					bbb[1]=buf;
												
					//do everything for second subtree
					buf = bbb[0];bbb[0]=p;
					doforsubtree(cid);
					bbb[0]=buf;

					distinside+=odisti;
					distoutside+=odisto;

					disti[cid]=odisti;
					disto[cid]=odisto;

				}					
				else
				{
					if(cid<0)
					{
						addleaf();
						leafit+=2;
					}
					else
					{
						addleaf();
						++leafit;
						addleaf();
						++leafit;
					}							
				}

				lst<<=1;

		}		
	};






};

//for statistics: the number of leaves that 
//had to be checked in 
//last invocation of getbox
template<class mynum>
int Momkdtree_templ<mynum>::impl::getbox::numnodesvisited=0;

template<class mynum>
int Momkdtree_templ<mynum>::impl::getball::numnodesvisited=0;

//build kdtree from dim-dimensional points stored in mynum
template<class mynum>
void Momkdtree_templ<mynum>::init(vector<mynum> &moments,int dimm)
{
	dim=dimm;
	impl::build(*this,moments);
}










/**
*get all leaves in the box.
	*\param moms: moments returned
	*\param box: an array with 2*dim entries
	*which are min,max-pairs for every dimension
*/
template<class mynum>
void Momkdtree_templ<mynum>::momIdsInBox(vector<int> &moms,const mynum*box) const
{
	#ifndef OLD_VC
	typename
	#endif
	impl::getbox bx(moms,*this,box);

	//cout<<bx.ret.size()<< " , " <<bx.numnodesvisited<<'/'<< leafids.size()<<endl;
}


/**
*get all leaves in the box.
	*\param moms: moments returned
	*\param box: an array with 2*dim entries
	*which are min,max-pairs for every dimension
*/
template<class mynum>
void Momkdtree_templ<mynum>::momsInBox(vector<mynum> &moms,const mynum*box) const
{
	assert(!leaves.empty());

	if(leaves.empty())
		return;

	vector<int> retids;

	#ifndef OLD_VC
	typename
	#endif
	impl::getbox bx(retids,*this,box);

	moms.resize(retids.size()*dim);
	#ifndef OLD_VC
	typename
	#endif

	vector<mynum>::iterator momit;
	vector<int>::iterator bxit;
	for(bxit=retids.begin(),momit=moms.begin();
			bxit!=retids.end();++bxit,momit+=dim)
	{
		copy(leaves.begin() + *bxit,
				 leaves.begin() + *bxit + dim,
				 momit);	
	}
}

/**
	*get all leaves in the Ball.
	*\param momids: start index of moment set in moment set array
	*\param center: an array with dim entries containing 
	*the center of the hyperball
	*\param radius: the radius of the hyperball
	*/
template<class mynum>
void Momkdtree_templ<mynum>::momIdsInBall(std::vector<int> &momids,
		const mynum*center,mynum radius) const
{
			#ifndef OLD_VC
			typename
			#endif

	impl::getball(momids,*this,center,radius);
}

template<class mynum>
Momkdtree_templ<mynum>::Momkdtree_templ(void)
{
	dim=0;height=0;
}

template<class mynum>
Momkdtree_templ<mynum>::~Momkdtree_templ(void)
{

}



template<class mynum>
inline void testbuildspeed()
{
	cout<<"testin speed"<<endl;

	static const int dim=32;
	vector<mynum> v;
	for(int i=0;i<100000*dim;++i)
		v.push_back(rand()*1.0/RAND_MAX);

	clock_t tstart=clock();
	for(int j=1;j<=10;++j){
		int ni=100000/v.size() +1;
		for(int i=0;i<ni;++i)
		{
			vector<mynum> vv=v;
			Momkdtree_templ<mynum> mkd(vv,dim);
		}
		float t=1.0*( (clock()-tstart)/CLOCKS_PER_SEC ) /(j*ni);

		cout<<"time per build tree:"<<t<<"s"<<endl;
		cout<<"divided by vv.size() log(vv.size()) : "
			<< 1e9*t/v.size()/log((double)v.size()/dim)  <<"ns"<<endl;

		

	}
}

template<class mynum>
inline void testqueryspeed()
{

	static const int dim=32;

	vector<mynum> v;
	for(int i=0;i<100000*dim;++i)
		v.push_back(rand()*1.0/RAND_MAX);

	mynum b[2*dim];
	
	
	cout<<"building tree"<<endl;

	vector<mynum> vv=v;
	Momkdtree_templ<mynum> mkd(vv,dim);

	cout<<"testing query speed"<<endl;
	vector<int> mi;

	clock_t tstart;
		

	int ni=(int)( 100000/log((double)v.size()) +1 ),j;

	tstart=clock();

	float boxrad=0.01;
	
	cout<<"\n1st number: time per box query,half edge length"<<boxrad<<"\n"
		<<"2nd number: time per node visited "<<endl;
	j=0;
	do{
		++j;
		for(int i=0;i<ni;++i)
		{
			int s =	int( (1.0*rand()/RAND_MAX)*v.size()/dim)*dim;
			
			if(s<0)s=0;if(s>=v.size())s=v.size()-dim;
			#ifndef OLD_VC
			typename
			#endif
			vector<mynum>::iterator vit
				=v.begin() + s;
						
			for(int i=0;i<dim;i++)
				b[2*i]=vit[i]-boxrad,b[2*i+1]=vit[i]+boxrad;
			
			mkd.momIdsInBox(mi,b);
		}

	}while(clock()-tstart<2000);

	float t=1.0*( clock()-tstart)/CLOCKS_PER_SEC /(j*ni);
	float tn=1.0*( clock()-tstart)/CLOCKS_PER_SEC /
		Momkdtree_templ<mynum>::impl::getbox::numnodesvisited;

	cout<<t*1e3<<"ms "<<tn*1e9 <<"ns"<<endl	;	

	float ballrad=boxrad;

	cout<<"\n1st number: time per ball query,radius"<<ballrad<<" dim "<<dim<<"\n"
		<<"2nd number: the time per node visited"<<endl;
	tstart=clock();
	Momkdtree_templ<mynum>::impl::getball::numnodesvisited=0;
	j=0;
	do{
		++j;
		for(int i=0;i<ni;++i)
		{
			int s =	int( (1.0*rand()/RAND_MAX)*v.size()/dim)*dim;
			
			if(s<0)s=0;if(s>=v.size())s=v.size()-dim;
			
			
			mkd.momIdsInBall(mi,&mkd.getleaves()[s],ballrad);
		}

	}while(clock()-tstart<2000);
	{
	float t=1.0*( clock()-tstart)/CLOCKS_PER_SEC /(j*ni);
	float tn = 1.0*( clock()-tstart)/CLOCKS_PER_SEC /
	Momkdtree_templ<mynum>::impl::getball::numnodesvisited;

	cout<<t*1e3<<"ms "<<tn*1e9 <<"ns"<<endl	;	
	}
}


template<class mynum>
void testMomkdtree()
{
		

	static const int n=5,dim=3;
	//mynum testv[n*dim]={1,1,2,2,3,3,4,4,5,5};

	vector<mynum> v;//(&testv[0],&testv[0]+n*dim);


	for(int i=0;i<100000*dim;++i)
		v.push_back(rand()*1.0/RAND_MAX);

	mynum b[2*dim];
	
	for(int i=0;i<2*dim;i+=2)
		b[i]=v[i/2]-0.3,b[i+1]=v[i/2]+0.3;
	
	cout<<"building tree"<<endl;

	vector<mynum> vv=v;
	Momkdtree_templ<mynum> mkd(vv,dim);



	vector<int> momidsib,momidsibtest,intersect;


	for(int it=0;it<10;++it)
	{
		mynum c[dim];
		for(int i=0;i<dim;++i)
			c[i]= rand()*1.0/RAND_MAX;

		mynum r=rand()*1.0/RAND_MAX;

		cout<<"testing ball query"<<endl;
		momidsib.clear();
		mkd.momIdsInBall(momidsib,c,r);
		sort(momidsib.begin(),momidsib.end());

		cout<<"testing it with a simple test"<<endl;
		momidsibtest.clear();
		for(int j=0;j<mkd.getleaves().size();j+=mkd.getDim())
		{
			mynum rsq=0;
			const mynum * set= &mkd.getleaves()[j];
			for(int i=0;i<dim;++i)
			{
				mynum x= set[i]-c[i];
				rsq += x*x;
			}
			if(rsq<r*r)momidsibtest.push_back(j);
		}

		intersect.resize(momidsib.size());

		intersect.resize(set_difference(momidsibtest.begin(),momidsibtest.end(),
			momidsib.begin(),momidsib.end(),intersect.begin())-intersect.begin() );
		if(intersect.size()	)
		{
			ofstream ballpt("testmomkdtree.asc");
			for(int i=0;i<intersect.size();++i)
			{
				const mynum*st = &mkd.getleaves()[intersect[i]];
				for(int j=0;j<3;++j)
					ballpt<<st[j]<<' ';
				ballpt<<'\n';	
			}
			ballpt.close();

			cout<<"error: "
				<<intersect.size() 
				<<" of "<<momidsibtest.size()<<" not in real function"<<endl;
			return;
		}
		else
			cout<<"correct"<<endl;

		intersect.resize(momidsib.size());

		intersect.resize(set_difference(momidsib.begin(),momidsib.end(),
			momidsibtest.begin(),momidsibtest.end(),intersect.begin())-intersect.begin() );
		if(intersect.size()	)
		{
			ofstream ballpt("testmomkdtree.asc");
			for(int i=0;i<intersect.size();++i)
			{
				const mynum*st = &mkd.getleaves()[intersect[i]];
				for(int j=0;j<3;++j)
					ballpt<<st[j]<<' ';
				ballpt<<'\n';	
			}
			ballpt.close();

			cout<<"error: "
				<<intersect.size() 
				<<" of "<<momidsib.size()<<" not in tree function"<<endl;
			return;
		}
		else
			cout<<"2nd correct"<<endl;

	}





	vector<mynum> moms;
	//mkd.momsInBox(moms,b);

	/*
	for(int i=0;i<moms.size();i+=dim,cout<<endl)
		for(int j=0;j<dim;++j)
				cout<<moms[i+j]<<' ';
*/
	vector<int> ms;		
	//mkd.momIdsInBox(ms,b);



	testqueryspeed<mynum>();



/*
	for(int i=0;i<ms.size();++i)
		cout<<ms[i]<<endl;
*/


}

#undef _TYPENAME_

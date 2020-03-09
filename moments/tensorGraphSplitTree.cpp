/*
	create "tensor graph split tree" which
	creates a network of tensor products and contractions  which describe the calculation of a set of moment invariants
	in an optimized manner

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

#ifndef TENSOR_GRAPH_SPLIT_TREE_H
#define TENSOR_GRAPH_SPLIT_TREE_H
#include<map>
#include<cassert>
#include<cmath>
#include<sstream>
#include<algorithm>
//#define _USE_MATH_DEFINES
//#include<math.h>

#include "binaryTreesVisitor.h"
#include "tensorGraphSplitTree.h"
#include "optimizedFold.h"

using namespace std;

	//remove links to children and evtl. delete them
	void tensorGraphSplitTree::freechild(nodeset::iterator it)
	{
					//free old nodes
					if(it != nodes.end())
					{
						node *a = *it;
						a->nparents--;
						if(a->nparents==0)
						{
							freechild(a->a);
							freechild(a->b);
							a->a = a->b = nodes.end();
							freenodes.push_back(a);
							nodes.erase(it);
						}
					}			
	}


	tensorGraphSplitTree::nodeset::
		iterator 
		tensorGraphSplitTree::
		insertchild(node * it,bool&nw)
	{
		std::pair<nodeset::iterator,bool> ret	= nodes.insert(it);
		
		//if it was already there,
		//new node can be freed
		if(!ret.second){
			freenodes.push_back(it);
		}		
		//update number of parents
		(**ret.first).nparents ++ ;

		nw=ret.second;
		return ret.first;
	}


		int tensorGraphSplitTree::splittree::initRecurse(const node*root)
		{
			stnode buf;
			
			if(root->tg.nnodes()==1)
				return -1;
			if(root->tg.nnodes()==2)
			{
				buf.a=buf.b=-1;
			}
			else
			{
				buf.a=initRecurse(*root->a);
				buf.b=initRecurse(*root->b);
			}
			buf.split=root->split;

			nodes.push_back(buf);

			return nodes.size()-1;
		}

		int tensorGraphSplitTree::splittree::init(const binaryTreesVisitor::node*root)
		{
			stnode buf;

			if(root->a)
				buf.a = init(root->a);
			else
				buf.a = -1;

			if(root->b)
				buf.b = init(root->b);
			else
				buf.b = -1;

			nodes.push_back(buf);
			
			return nodes.size()-1;					
		}

		void tensorGraphSplitTree::splittree::recurseRecreate(int nodeid,node*root) const
		{
			if(nodeid<0)
				return;

			const stnode & buf = nodes[nodeid];	
			bool nwa,nwb;
			root->setSplit(buf.split,nwa,nwb);			

			//only set their splits if newly created.
			if( buf.a != -1 && nwa)
				recurseRecreate(buf.a,*root->a);			
			if( buf.b != -1 && nwb)
				recurseRecreate(buf.b,*root->b);			
		}



		//the effort for the current split
		//(a raw estimate)
		size_t tensorGraphSplitTree::node::getCurrentEffort()
		{
			int sumord =0; 
			for(int i=0; i < tg.nnodes(); ++i)
				sumord += tgs.orderpernodetype[tg.nodelabels[i]];

			sumord -= tg.edges.size()*2;

			//add effort for te edes summed over:
			sumord += s.edgessplit.size();

			int dim=tgs.dim;
			size_t ret = dim;
			for(int i=1;i<sumord;++i)
				ret*=dim;
			return ret;
		}


		void tensorGraphSplitTree::node::setSplit(unsigned newsplit,bool&nwa,bool&nwb)
		{
					tgs.freechild(a); 
					tgs.freechild(b);
					
					node *pa =tgs.getNode();
					node *pb =tgs.getNode();
					
					s.init(pa->tg,pb->tg,tg,newsplit);
					

					std::pair<nodeset::iterator,bool> ret;
					
					//normalize new children and store permutations
					pa->tg.normalize(&perma);
					pb->tg.normalize(&permb);

					//insert them into the tree
					a=tgs.insertchild(pa,nwa);	
					b=tgs.insertchild(pb,nwb);			
					split=newsplit;
		}


		


		size_t tensorGraphSplitTree::node::findOptimalSplit()
		{

			if(tg.nnodes() == 1 )
				return 0;
			
			size_t mineffort = 0xffffffff;
			unsigned nnodebits=(1<<tg.nnodes())-1;
			
			//iterate through the possible splits
			for(unsigned s=1; s < tg.nnodes();++s)
			{
				//use the numerically smaller instance of
				//the split and its inverse
				if( s <= (s ^ nnodebits ) )
				{
					bool nwa,nwb;
					
					cout<<oct<<s<<dec<<flush;

					//split current tree in two
					setSplit(s,nwa,nwb);
					
					//do't use splits where one component is unconnected afterwards.
					if( nwa &&  ! (**a).tg.connectedComps() 
						|| nwb && !(**b).tg.connectedComps() )
					{
						continue;
					}

					//the ffort for re-joinin the current split
					size_t effort = getCurrentEffort();
					cout<<'<'<<effort<<'>';
					//only count effort once
					//for creating the sub-splits
					if(nwa) {
						cout<<'(';
						effort += (**a).findOptimalSplit();
						cout<<')';
					}
					if(nwb) 
					{
						cout<<'(';
						effort += (**b).findOptimalSplit();
						cout<<')';
					}
					
					if(effort<mineffort)
					{
						mineffort=effort;
						optimalsplit.init(this);
					}
				}				
			}

			optimalsplit.recreate(this);
			
			return mineffort;

		}


	inline void checkSanity(const tensorGraph&tg,const vector<int>&orderpernodetype)
	{
		vector<int>::const_iterator it;
		for(it=tg.nodelabels.begin();it!=tg.nodelabels.end();++it)
			if(*it>=orderpernodetype.size() || *it<0)
				{
					cerr<<"Error: label"<<*it<<" not found in orderpernodetype"<<endl;
					throw;
				}
		if(!tg.connectedComps())
		{
			cerr<<"Error: tensor Graph has multiple connected components!"<<endl;
			throw;
		}
	}


	size_t tensorGraphSplitTree::findOptimalSplit( const tensorGraph& tg, 
		const std::vector<int> &orderpernodetype_, int dim_)
	{
		
		dim = dim_;
		orderpernodetype=orderpernodetype_;		
		//free all previous nodes used:
		freenodes.insert(freenodes.end(),nodes.begin(),nodes.end());
		//clear set of nodes
		nodes.clear();			
		roots.clear();

		checkSanity(tg,orderpernodetype);
		
		node * root = getNode();		
		root->tg=tg;

		root->tg.normalize();

		bool nw;		
		roots.push_back(insertchild(root,nw));
		return root->findOptimalSplit();		
	}





	size_t tensorGraphSplitTree::findOptimalSplits( const vector<tensorGraph>& tgs, 
		const std::vector<int> &orderpernodetype_, int dim_)
	{
		

		dim = dim_;
		orderpernodetype=orderpernodetype_;		
		//free all previous nodes used:
		freenodes.insert(freenodes.end(),nodes.begin(),nodes.end());
		//clear set of nodes
		nodes.clear();
			
		roots.resize(tgs.size());

		size_t effort=0;
		for(int i=0;i<tgs.size();++i)
		{
			checkSanity(tgs[i],orderpernodetype);

			node*root= getNode();
			root->tg=tgs[i];
			root->tg.normalize();
			bool nw;
			roots[i]=insertchild(root,nw);
			if(nw)
				effort += root->findOptimalSplit();
		}

		return effort;
	}



	void tensorGraphSplitTree::printPS(ostream&out,
		const std::vector<std::vector<int> >& groupSizPerType) const
	{
		for(nodeset::const_iterator it=nodes.begin(); it!=nodes.end();++it)
		{
			(**it).tg.printPS(out,groupSizPerType);
		}	
	}


/*****************************************************************************/
/*** now: tensorgraphsplitree2, splittin only one edge at atime.             */
/*****************************************************************************/

	//remove links to children and evtl. delete them
	void tensorGraphSplitTree2::freechild(nodeset::iterator it)
	{
					//free old nodes
					if(it != nodes.end())
					{
						node *a = *it;
						a->nparents--;
						if(a->nparents==0)
						{
							freechild(a->a);
							freechild(a->b);
							a->a = a->b = nodes.end();
							freenodes.push_back(a);
							nodes.erase(it);
						}
					}			
	}


	tensorGraphSplitTree2::nodeset::
		iterator 
		tensorGraphSplitTree2::
		insertchild(node * it,bool&nw)
	{
		std::pair<nodeset::iterator,bool> ret	= nodes.insert(it);
		
		//if it was already there,
		//new node can be freed
		if(!ret.second){
			freenodes.push_back(it);
		}		
		//update number of parents
		(**ret.first).nparents ++ ;

		nw=ret.second;
		return ret.first;
	}


		int tensorGraphSplitTree2::splittree::initRecurse(const node*root)
		{
			stnode buf;
			
			if(root->edgeid>=0)
			{
				if((**root->a).edgeid>=0)
					buf.a=initRecurse(*root->a);
				else
					buf.a=-1;
				
				if(root->compsb && (**root->b).edgeid>=0)
					buf.b=initRecurse(*root->b);
				else
					buf.b=-1;
			}
			else
				buf.a=buf.b=-1;
			
			buf.edgeid=root->edgeid;

			nodes.push_back(buf);

			return nodes.size()-1;
		}

		int tensorGraphSplitTree2::splittree::init(const binaryTreesVisitor::node*root)
		{
			stnode buf;

			if(root->a)
				buf.a = init(root->a);
			else
				buf.a = -1;

			if(root->b)
				buf.b = init(root->b);
			else
				buf.b = -1;

			nodes.push_back(buf);
			
			return nodes.size()-1;					
		}

		void tensorGraphSplitTree2::splittree::recurseRecreate(int nodeid,node*root) const
		{
			if(nodeid<0)
				return;

			const stnode & buf = nodes[nodeid];	
			bool nwa,nwb;
			root->setSplit(buf.edgeid,nwa,nwb);			

			//only set their splits if newly created.
			if( buf.a != -1 && nwa)
				recurseRecreate(buf.a,*root->a);			
			if( buf.b != -1 && nwb)
				recurseRecreate(buf.b,*root->b);			
		}



		//the effort for the current split
		//(a raw estimate, mostly an over-estimation)
		size_t tensorGraphSplitTree2::node::getCurrentEffort()
		{
			int freeedges =0; 
			for(int i=0; i < tg.nnodes(); ++i)
				freeedges += tgs.orderpernodetype[tg.nodelabels[i]];

			freeedges -= ( tg.edges.size() ) * 2;

			int dim=tgs.dim;
			size_t ret = dim;
			//effort for free indices + one for contraction of a single one
			for(int i=1;i<freeedges+1;++i)
				ret*=dim;

			return ret;
		}

		
		void tensorGraphSplitTree2::node::tryjoinwithchild()
		{
			//if there is not exactly one child
			//or that child has more than one parent,
			//don't do anything
			if( a==tgs.nodes.end() ||	
					b!=tgs.nodes.end() ||
					(**a).nparents!=1 )
					return;
			
			node  &n = **a;
			node* na=n.pa(),*nb=n.pb();

			//inverse of perma
			vector<int> iperma;
			//re-use mem of permb
			iperma.swap(permb);
			
			iperma.resize(perma.size());
			for(size_t i=0;i<perma.size();++i)
				iperma[perma[i]]=i;

			//insert edges split into this node with corrected 
			//nodeids
			for(size_t i=0;i<n.s.edgessplit.size();++i)
			{
				tensorGraph::edge::data d=n.s.edgessplit[i].d;
				
				tensorGraph::edge nw;

				nw.d.u = d.u;
				nw.d.v = d.v;

				nw.d.a = iperma[d.a];
				nw.d.b = iperma[d.b];
					
				s.edgessplit.push_back(nw);					
			}				


			if( !na )  //no sub-children
			{			
				//do nothing
			}else if (!nb) //one sub-child
			{
				
				//concatenate permutations of this node and it's subnode
				vector<int> operm;
				operm.swap(perma);
				tensorGraph::concatperm(perma,operm,n.perma);

				//now replace sub-child and erase node n.
				tgs.nodes.erase(a);
				a=n.a;				
				
				//put n into the free nodes list
				//but first set the children to non-existent.
				n.a = n.b = tgs.nodes.end();
				tgs.freenodes.push_back(&n);
				compsb=0;

				//because this node has still only one sub-child,
				//try join again
				tryjoinwithchild();

			}else //two sub-children
			{
				//apply operm to node id map:
				s.nodeidmap.resize(n.s.nodeidmap.size());
				for(size_t i=0;i<s.nodeidmap.size();++i)
					s.nodeidmap[i] = n.s.nodeidmap[perma[i]];
				
				compsb=0;
				for(size_t i=0;i<tg.nnodes();++i)
					compsb |= ((n.compsb>>perma[i])&1)<<i ;
			
				//now replace sub-children and erase node n.
				tgs.nodes.erase(a);
				a=n.a;
				b=n.b;
				

				perma.swap(n.perma);
				permb.swap(n.permb);

				//put n into the free nodes list,
				//but first set the children to non-existent.
				n.a = n.b = tgs.nodes.end();
				tgs.freenodes.push_back(&n);
			}

			

		
		}

		void tensorGraphSplitTree2::node::setSplit(int newsplit,bool&nwa,bool&nwb)
		{
					tgs.freechild(a); 
					tgs.freechild(b);
					edgeid=newsplit;
					
					compsb = tg.notConnectedWithFirst(newsplit);
					//if edge cut doesn't lead to split:
					if(!compsb )
					{
						node *pa =tgs.getNode();
						pa->tg=tg;
						s.edgessplit.clear();
						s.edgessplit.resize(1,tg.edges[newsplit]);
					

						pa->tg.edges.erase(pa->tg.edges.begin()+newsplit);
						pa->tg.normalize(&perma);
						a = tgs.insertchild(pa,nwa);		
						nwb = false;
						b = tgs.nodes.end();
						return;
					}
					else
					{	 //two connected comps: 				
						node *pa =tgs.getNode();
						node *pb =tgs.getNode();

						s.init(pb->tg, pa->tg, tg, compsb);

						//normalize new children and store permutations
						pa->tg.normalize(&perma);
						pb->tg.normalize(&permb);

						//insert children into tree
						a = tgs.insertchild(pa,nwa);		
						b = tgs.insertchild(pb,nwb);		
					}					
		}


		


		size_t tensorGraphSplitTree2::node::findOptimalSplit()
		{
			if(tg.edges.size()==0)
				return 0;

			size_t mineffort = 0xffffffff;
			unsigned nnodebits=(1<<tg.nnodes())-1;
			
			//iterate through the possible split edge
			for(int i=0; i < tg.edges.size();++i)
			{
					//overjump equal edges:
					while( i+1 < tg.edges.size() && tg.edges[i]==tg.edges[i+1] )++i;

					bool nwa,nwb;
					
					//cout<<oct<<i<<dec<<flush;

					++tgs.currentiter;
					if((tgs.currentiter%100)==0)
						cout<<tgs.currentiter<<" ( "<<100.*tgs.currentiter/tgs.maxiter<<"%)"<<endl;

					//split current tree in two
					setSplit(i,nwa,nwb);

					//the ffort for re-joinin the current split
					size_t effort=0;
					
				  effort+= getCurrentEffort();
					
					//cout<<'<'<<effort<<'>';
					//only count effort once
					//for creating the sub-splits
					if(nwa) {
						//cout<<'(';
						effort += (**a).findOptimalSplit();
						//cout<<')';
					}
					if(nwb) 
					{
						//cout<<'(';
						effort += (**b).findOptimalSplit();
						//cout<<')';
					}
					
					if(effort<mineffort)
					{
						mineffort=effort;
						optimalsplit.init(this);
					}
			}

			optimalsplit.recreate(this);
			
			return mineffort;

		}


	void tensorGraphSplitTree2::reset()
	{
		freenodes.insert(freenodes.end(),nodes.begin(),nodes.end());
		freenodes.insert(freenodes.end(),roots.begin(),roots.end());

		//clear set of nodes
		nodes.clear();			
		roots.clear();			
	}

	size_t tensorGraphSplitTree2::findOptimalSplit( const tensorGraph& tg)
	{
		
		checkSanity(tg,orderpernodetype);

		node * root = getNode();		
		root->tg=tg;
		root->tg.normalize();
		roots.push_back(root);

		currentiter=0;
		maxiter=0;
		for(int i=1;i<=tg.edges.size();++i)
		{maxiter+=1;maxiter*=i;}

		return root->findOptimalSplit();		
	}





	size_t tensorGraphSplitTree2::findOptimalSplits( const vector<tensorGraph>& tgs)
	{
		reset();
					
		roots.resize(tgs.size());

		size_t effort=0;
		for(int i=0;i<tgs.size();++i)
		{
			checkSanity(tgs[i],orderpernodetype);

			node*root= getNode();
			root->tg=tgs[i];
			root->tg.normalize();
			roots[i]=root;
			effort += root->findOptimalSplit();
		}

		return effort;
	}
	
	void tensorGraphSplitTree2::node::updateParentSetOf(const node* child ,tparentset &parents)const
	{
	
			set<const node*> & ps = parents[child];

			if(ps.empty())
				child->addToParentSet(parents);

			ps.insert(this);
	
	}

	void tensorGraphSplitTree2::node::addToParentSet(tparentset & parents) const
	{
		if( a!= tgs.nodes.end() )
		{
			updateParentSetOf(*a,parents);
		}

		if( b!= tgs.nodes.end() )
		{
			updateParentSetOf(*b,parents);
		}
	
	}


	//update the parent vector for each node
	void tensorGraphSplitTree2::getParentSet(
		tparentset & parents) const
	{
		for(int i=0;i<roots.size();++i)
			roots[i]->addToParentSet(parents);

		/*
		tparentset::iterator it;
		for(it=parents.begin();it!=parents.end();++it)
			if( 
				it->first->nparents 
				<
				it->second.size()
				)
				cout << it->first->nparents <<" !="<< it->second.size() ;
				*/

	}

		//for every level, a map with nodes ids
	void tensorGraphSplitTree2::createGraphLevelsBottomup
		(	vector< set< const node* > >& levels	) const
	{
		
		tparentset parents;		
		getParentSet(parents);		
		
		set<const node*> nchildren[3];		
		
		for(nodeset::const_iterator it=nodes.begin(); it!=nodes.end();++it)
		{
			nchildren[(**it).numDistinctChildren()].insert( *it ) ;
		}

		for(int i=0;i<roots.size();++i)
			nchildren[roots[i]->numDistinctChildren()].insert(roots[i]) ;

		while( !nchildren[0].empty() )
		{
			
			levels.push_back(set<const node*>());
			nchildren[0].swap(levels.back());
			set<const node*> & l = levels.back();						

			set<const node*>::iterator lit;
			for(lit=l.begin();lit!=l.end();++lit)
			{
				set<const node*> &s = parents[*lit];
				set<const node*>::iterator pit;
				
				//update num. children:
				for(pit=s.begin();pit!=s.end();++pit)
				{
					if(nchildren[1].erase(*pit))
						nchildren[0].insert(*pit);	
					else if(nchildren[2].erase(*pit))
						nchildren[1].insert(*pit);		
					else
						assert(false);
				}
			}

		}				
	}

	void tensorGraphSplitTree2::createGraphLevelsTopdown(std::vector< std::set< const node* > >& levels	) const
	{
		map<int,set<const node*> > nodesWithParentNum;
		map<const node*,int > parentsPerNode;
	
		for(nodeset::const_iterator it=nodes.begin(); it!=nodes.end();++it)
		{
			nodesWithParentNum[(**it).nparents].insert(*it);
			parentsPerNode[*it]=(**it).nparents;
		}

		set<const node*> &n0 =nodesWithParentNum[0];
		
		for(int i=0;i<roots.size();++i)
		{
			n0.insert(roots[i]);
			parentsPerNode[roots[i]] = roots[i]->nparents;
		}

		while( !n0.empty() )
		{
			levels.push_back(set<const node*>());
			set<const node*> & l = levels.back();						
			l.swap(n0);

			set<const node*>::iterator lit;
			for(lit=l.begin();lit!=l.end();++lit)
			{
				const node & n= (**lit) ;
				if( n.a!=nodes.end() )
				{
					int & ppn = parentsPerNode[*n.a];
					nodesWithParentNum[ppn].erase(*n.a);
					--ppn;
					nodesWithParentNum[ppn].insert(*n.a);
				}
				if( n.b!=nodes.end() )
				{
					int & ppn = parentsPerNode[*n.b];
					assert(ppn!=0);
					nodesWithParentNum[ppn].erase(*n.b);
					--ppn;
					nodesWithParentNum[ppn].insert(*n.b);
				}
			
			}			
		}
	
	}


	struct tensorGraphSplitTree2::npinfo
	{
		float x,y,w;	
		string name;
	};
	
	
	inline bool mycmp(const pair<float,const tensorGraphSplitTree2::node*>&a,
		const pair<float,const tensorGraphSplitTree2::node*>&b)

	{
		return a.first<b.first || 
			a.first==b.first && a.second->tg < b.second->tg;
	}



	void tensorGraphSplitTree2::printPS(ostream&out) const
	{

		//don't draw connections to single-node-tensors:
		bool nosingle=false;;

		//create graph bottom-up:
		bool bottomup=true;

		//scale the graphics
		float scale=1;

		out<<scale<< " dup scale\n";
		float pagewidth =  60 / 2.54 *72 /scale ;

		//a set with all nodes along with the child number

		vector< set<const node*> > levels;		

		if(bottomup)
			createGraphLevelsBottomup(levels);
		else
			createGraphLevelsTopdown(levels);

				

		//additional information for every node
		map< const node*, npinfo > nodeInfo;

		float y=0;
		//create procedure that draw the graphs:
		//out<<"\n/drawgraphs {\n";


		vector<pair< float, const node* > > desiredx;

		for(int i=0;i<levels.size();++i)
		{
			out<<"matrix currentmatrix\n";			
			out<<"0 "<<y<<" translate\n";
			float x=0;
			float maxw=0;
			set<const node*>::iterator sit;

			desiredx.clear();
			for(sit=levels[i].begin(); sit!=levels[i].end();++sit)
			{
				const node&n = **sit;
				float desx =  0;
				if(bottomup)
				{
					int nc=n.numDistinctChildren();
					if(nc==2)
						desx = 0.5*(nodeInfo[*n.a].x+nodeInfo[*n.b].x);
					else if(nc==1)
						desx = nodeInfo[*n.a].x;
				}

				desiredx.push_back(pair<float,const node*>(desx,*sit));
			}

			sort(desiredx.begin(),desiredx.end(),mycmp);
			
			vector<pair<float,const node*> >::iterator dit;

			int id=0;
			for(dit=desiredx.begin(); dit!=desiredx.end();++dit,++id)
			{
				npinfo& ni = nodeInfo[dit->second];
				ni.w = dit->second->tg.printPS(out,groupSizPerType);
				ni.x = x+ni.w/2;
				ni.y = y+ni.w/2;

				
				ostringstream nn; 
				nn<<dit->second->id<<char('A'+i)<<id;
				ni.name=nn.str();
				nn.clear();


				if(maxw<ni.w)maxw=ni.w;
				x+=ni.w;

				out<<"\n5 0 translate\n\n";
				x+=5;

				//line break
				if(x> pagewidth-maxw)
				{
					y+=maxw+10;
					x=0;maxw=0;
					out<<"setmatrix\n";
					out<<"matrix currentmatrix\n";
					out<<"0 "<<y<<" translate\n";
				}

			}
			y += maxw+100;
			out<<"setmatrix\n";
		}


		ostringstream os;
		//set moment names:
		for(int i=0;i<roots.size();++i)
		{
			os.str("");
			npinfo &np=nodeInfo[ roots[i] ];
			os<<"M"<<i;
			np.name =os.str();
			os.clear();
		}


		
		out<<
		"\n\n %draw graph names , draw connections between graphs\n"
		"/Times-Roman findfont 10 scalefont setfont\n"
		"0 setgray\n";

		set<const node*>::iterator it;
		for(int i=0;i<levels.size();++i)
		{
			set<const node*>::iterator it;
			for(it=levels[i].begin(); it!=levels[i].end();++it)
			{
				const node & n = **it;
				npinfo & p = nodeInfo[*it];
				
				int nc=n.numDistinctChildren();
				//draw name on top
				if(nc==0)
				{
					out<< p.x<<' '<<p.y+p.w/2<<" moveto\n";
					out<<"("<<p.name<<") show\n";
				}

				if(nc)
				{


					npinfo &a = nodeInfo[*n.a ];

					//draw name on top
					if(nc==1)
					{
						out<< p.x<<' '<<p.y+p.w/2<<" moveto\n";
						out<<"("<<p.name<<"\\("<<a.name<<"\\)) show\n";
					}
				  
					float dx=p.x-a.x, dy=p.y-a.y;	
					float invnorm = 1.0f / sqrt(dx*dx+dy*dy);			
					dx*=invnorm;dy*=invnorm;
					
					float asiz = a.w/2;
					float psiz = p.w/2;
					float 
							ax= a.x  + dx*asiz,
							ay= a.y + dy*asiz,
							px= p.x  - dx*psiz,
							py= p.y - dy*psiz;

					if(!nosingle || (**n.a).tg.nnodes()>1 )
					{
						if(n.a==n.b)
						{

							float mul=10;
							dx*=mul; dy*=mul;
							out<<px<<' '<<py<<" moveto\n";
							out<<px-dx+dy <<' '<<py-dy-dx <<"\n";
							out<<ax+dx+dy <<' '<<ay+dy-dx <<"\n";
							out<<ax<<' '<<ay<<" curveto\n";
							out<<ax+dx-dy <<' '<<ay+dy+dx <<"\n";
							out<<px-dx-dy <<' '<<py-dy+dx <<"\n";
							out<<px<<' '<<py<<" curveto\n";						
						}
						else
						{
							out<<px<<' '<<py<<" moveto\n";
							out<<ax<<' '<<ay<<" lineto\n";
						}
					}

					if(nc>1 && ( !nosingle || (**n.b).tg.nnodes()>1 ) )
					{

						npinfo &b = nodeInfo[* n.b ];
						
						//draw name on top right
						out<< p.x<<' '<<p.y+p.w/2<<" moveto ";
						out<<"("<<p.name<<"\\("<<a.name<<','<<b.name<<"\\))  show\n";
	
						dx=p.x-b.x, dy=p.y-b.y;	
						
						float invnorm = 1.0f / sqrt(dx*dx+dy*dy);						
						dx*=invnorm;dy*=invnorm;
					
						float bsiz = b.w/2;

						float 
							bx= b.x  + dx*bsiz,
							by= b.y + dy*bsiz,
							px= p.x  -dx*psiz,
							py= p.y -dy*psiz;
						out<<px<<' '<<py<<" moveto ";
						out<<bx<<' '<<by<<" lineto\n";
					}
				}
			}
		}
		out<<"1 setlinewidth .5 setgray stroke \n\n";
		//out<<"drawgraphs\n";

	}



/**convert a node permutation to an index permutation
*\param indperm: the new index permutation
*\param perm: the permutation of tensor nodes
*\param indsPerNode: number of indices per tensor node
* before the permutation
*/
	void tensorGraphSplitTree2::
		perm2indperm( vector<int>&indperm,
									const vector<int> &perm,
									const vector<int > &indsPerNode)
{
		assert(perm.size()==indsPerNode.size());
		vector<int> startinds(perm.size()+1);
		
		//invert perm
		vector<int> invperm(perm.size());
		for(int i=0;i<invperm.size();++i)
			invperm[perm[i]]=i;


		int numinds=0;
		startinds[0]=0;
		for(int i=0;i<perm.size();++i)
		{
			numinds+=indsPerNode[i];
			startinds[i+1] = numinds;
		}
	
			
		indperm.resize(numinds);		
		int numi=0;		
		for(int i=0;i<perm.size();++i)
		{
			int ip=invperm[i];
			int j = startinds[ip], jend = startinds[ip+1];
			while(j!=jend)
			{
				indperm[j] = numi;
				++numi;++j;
			}
		}			
}

//convert group sizes to offsets that point after the end 
//of one index group
inline void freeInds2Offsets(vector<vector<int> > &freeinds)
{
	int off=0;
	
	for(vvi::iterator i=freeinds.begin();
			i!=freeinds.end();
			++i)
	{
		for(vi::iterator j=i->begin();j!=i->end();++j)
		{
			off += *j;
			*j = off;
		}
	}

}

inline void permute(vector<vector<int> > &orig, const vector<int>&perm)
{
	vector< vector<int> > buf(orig.size());
	
	for(int i=0;i<buf.size();++i)
		buf[perm[i]].swap(orig[i]);

	buf.swap(orig);

}

/**inverse permutation*/
inline void invpermute(vector<vector<int> > &orig, const vector<int>&perm)
{
	vector< vector<int> > buf(orig.size());
	
	for(int i=0;i<buf.size();++i)
		buf[i].swap(orig[perm[i]]);

	buf.swap(orig);

}

inline bool isUnitPerm(const vector<int>&perm)
{
	for(int i=0;i<perm.size();++i)
		if(perm[i]!=i)
			return false;
	return true;
}


/**group the free index number by node
*/
void tensorGraphSplitTree2::groupFreeInds(vector<int>&indsPerNode,
													const vector<vector<int> >&fia)
{

	indsPerNode.resize(fia.size());
	for(int i=0;i<indsPerNode.size();++i)
	{
		const vector<int> &fiai=fia[i];
		int s=0;
		for(vi::const_iterator j=fiai.begin();j!=fiai.end();++j)
			s+=*j;
		indsPerNode[i]=s;
	}
}

repeatTensor * tensorGraphSplitTree2::node::createRepeatTensor
(const repeatTensor* ra,const repeatTensor* rb) const
{
	
	//permutation of indices
	vi indperm[2];
	
	//the number of free indices for every index group in a node
	vector<vector<int> >	freeindspg[2];

	//for debugging, store them after being converted to offsets:
	vector<vector<int> >	ofinds[2],ofi[2];

	//number of free indices of the tensors represented by the graphs
	vector<int> freeindspt[2];

	const tensorGraph * graphs[2]={ &((**a).tg) , &((**b).tg) };

	const vector<int> * perms[2]={ &perma, &permb };

	for(int i=0;i<2;++i){

		graphs[i]->getFreeIndsPerGroup
			(freeindspg[i],tgs.groupSizPerType);

		//restore original index groups
		invpermute(freeindspg[i], *perms[i] );	

		//calculate num free inds of every tensor node
		groupFreeInds(freeindspt[i],freeindspg[i] );

		//convert the permutation of nodes to
		//a permutation of indices
		perm2indperm(indperm[i], *perms[i] , freeindspt[i]);

		//convert the number of free indices per group
		//to the start index of the group plus the group size
		ofi[i]=freeindspg[i];

		freeInds2Offsets(freeindspg[i]);

		if(!isUnitPerm(*perms[i]) )
		{
			cout<<"<perm at Id= "<<(i==0?ra->id2name():rb->id2name())<<'>';
		}
		if( !isUnitPerm(indperm[i]) )
		{
			cout<<"<indperm at Id= "<<(i==0?ra->id2name():rb->id2name())<<'>';
		}
	}

	if(freeindspg[0].back().back()!=ra->ord())
	{
		assert(false);
		cerr<<"num inds wrong!"<<endl;
		throw;
	}

	//create "global" indices (the index offsets for rb
	//indices has ra.ord added)
	for(int i=0;i<freeindspg[1].size();++i){
		vector<int>& fi=freeindspg[1][i];
		for(int j=0;j<fi.size();++j)
			fi[j]+=ra->ord();
	}
	ofinds[0]=freeindspg[0];
	ofinds[1]=freeindspg[1];

	pr buf;
	vpr pairs;
	for(int i=0;i<s.edgessplit.size();++i)
	{
		tensorGraph::edge::data e = s.edgessplit[i].d;
		int
			ainb = (compsb>>e.a)&1,
			binb = (compsb>>e.b)&1;
		
		buf.a= --freeindspg[ainb][s.nodeidmap[e.a]][e.u];
		buf.b= --freeindspg[binb][s.nodeidmap[e.b]][e.v];

		if(buf.a>buf.b)
		{
			int xch=buf.a;buf.a=buf.b;buf.b=xch;
		}

		assert(buf.a!=buf.b && buf.a>=0 &&buf.b>=0);

		pairs.push_back(buf);
	}	

	repeatTensor *ret=new RTpartialFoldProd(ra,indperm[0],rb,indperm[1],pairs);

	return ret;
}

repeatTensor * tensorGraphSplitTree2::node::createRepeatTensor
	(const repeatTensor * ra) const
{
	vpr pairs;
	
	vi indperm;
	
	vector<vector<int> >
		freeindspg;

	(**a).tg.getFreeIndsPerGroup(freeindspg,tgs.groupSizPerType);
	
	invpermute(freeindspg,perma);	

	vector<int> freeindspt;
	groupFreeInds(freeindspt,freeindspg);
	freeInds2Offsets(freeindspg);

	pr buf;
	for(int i=0;i<s.edgessplit.size();++i)
	{
		tensorGraph::edge::data d = s.edgessplit[i].d;
		buf.a = --freeindspg[d.a][d.u];
		buf.b = --freeindspg[d.b][d.v];
		assert(buf.a!=buf.b && buf.a>=0 &&buf.b>=0);
		pairs.push_back(buf);
	}

	perm2indperm(indperm,perma,freeindspt);

	repeatTensor*ret=new RTpartialFold(ra,indperm,pairs);
	return ret;
}

void tensorGraphSplitTree2::updateNodeIds()
{
	nodeset::const_iterator ndit;
	
	//assign nodeids in sort order of the set.
	int id=0;
	for(ndit=nodes.begin();ndit!=nodes.end();++ndit,++id)
		(**ndit).id=id;
	
	//assign ids of roots
	for(int i=0;i<roots.size();++i,++id)
		roots[i]->id=id;

	for(int i=0;i<freenodes.size();++i,++id)
		freenodes[i]->id=id;
}




void tensorGraphSplitTree2::createRepeatTensors
( const vector< const repeatTensor* >& source, 
	vector< repeatTensor* > & intermediate,
	vector< repeatTensor* > & moments) const
{

		vector< set<const node*> > levels;	

		createGraphLevelsBottomup(levels);

		map<const node*,const repeatTensor*> rt ;

		map<const node*,int> rootid;

		for(int i=0;i<roots.size();++i)
			rootid[roots[i]]=i;

		moments.resize(roots.size());
		
		//insert node repeat tensors:
		set<const node*> &l = levels[0];
		set<const node*>::iterator it;
		for(it=l.begin();it!=l.end();++it)
		{
			rt[*it] = &*source[(**it).tg.nodelabels[0]];
		}

#define TESTPOLYS

#ifdef TESTPOLYS
		//for debugging:
		tensor<polyn> buf;
		vector<tensor<polyn> > tensPerNodeLabel(source.size());
		for(int i=0;i<source.size();++i)
		{
			if(source[i])
				tensPerNodeLabel[i] = *source[i];		
		}

#endif

		intermediate.clear();

		vpr pairs;

		for(int i=1;i<levels.size();++i)
		{
			set<const node*> &li=levels[i];
#ifdef TESTPOLYS			
			cout<<"level "<<char('A'+i)<<':'<<endl;
#endif
			for(it=li.begin();it!=li.end();++it)
			{

				const node &n = **it;				
				
				//nodes with only one parent should not be created.
				//(they will be joined with their parent)
				//if(n.nparents==1) continue;
				
				repeatTensor *r;				
				
				node* pa=n.pa(), *pb = n.pb();

#ifdef TESTPOLYS
				if(n.id==29||n.id==31||n.id==27)
				{
					cout<<'.';
				}
#endif

				if(pa&&pb)
					r=n.createRepeatTensor(rt[pa],rt[pb]);
				else
					r=n.createRepeatTensor(rt[pa]);

#ifdef TESTPOLYS

				
				cout<<"rt ids:";
				cout<<r->id<<'(';					
				cout<<rt[pa]->id;
				if(pb)cout<<','<<rt[pb]->id;
				cout<<")";

				cout<<" nd ids:";
				cout<<n.id<<'(';					
				cout<<pa->id;
				if(pb)cout<<','<<pb->id;
				cout<<")";

				if(n.nparents==0)
					cout<<"= M"<<rootid[*it];
				cout<<endl;

				n.tg.createPolynTensor(buf,
						tensPerNodeLabel, this->groupSizPerType );
				if(!( *r == buf ) )
				{
					cout<<"error! ";	
					cout<<r->id2name()<<endl;
					//cout<<"r=\n"<< *r <<'\n';
					//cout<<"buf=\n"<<buf<<'\n';

					cout<<" order is:"<<buf.ord();

					cout<<"tensorGraphs:\n";
					cout<<n.tg<<'\n'<<pa->tg;
					if(pb)cout<<'\n'<<pb->tg<<endl;




					cout<<"\n--------\n";
				}
#endif
				if(n.nparents==0)
					moments[rootid[*it]]=r;
				else
				{

					intermediate.push_back(r);				
					rt[*it]=r;
				}
			}
		}

}



void tensorGraphSplitTree2::joinSingleParentNodes()
{
	nodeset::iterator nit;
	for( nit = nodes.begin();nit!=nodes.end();++nit )
	{
		(**nit).tryjoinwithchild();		
	}

	for(size_t i=0;i<roots.size();++i)
		roots[i]->tryjoinwithchild();	
}

void tensorGraphSplitTree2::printForC
(ostream&out,const vector< const repeatTensor* >& source) const
{

	
	vector< repeatTensor* >  intermediate;
	vector< repeatTensor* >  moments;

	createRepeatTensors(source,intermediate,moments);


	for(size_t i=0;i<source.size();++i)
		if(source[i])
			source[i]->printForC(out);

	for(size_t i=0;i<intermediate.size();++i)
	{
		intermediate[i]->printForC(out);
		delete intermediate[i];
	}


	for(size_t i=0;i<moments.size();++i)
	{
		ostringstream oss;
		oss<<"M["<<i<<"]"<<ends;
		moments[i]->printForC(out,oss.str().c_str());
		delete moments[i];
	}

}


void tensorGraphSplitTree2::printGraphML(std::ostream&out) const
{
	out<<"<graph id=\"G\" edgedefault=\"directed\" >\n";
	out<<" <key id=\"M\" for=\"node\" attr.name=\"number of invariant\" "
		"attr.type=\"int\" />\n";


	//output nodes:		
	nodeset::const_iterator it;
	for(it= nodes.begin();it!=nodes.end();++it)
	{
		out<<"<node id=\""<<(**it).id<<"\">\n";
		(**it).tg.printGraphML(out,groupSizPerType);		
		out<<"</node>\n";
	}

	for(int i=0;i<roots.size();++i)
	{
		out<<"\n<node id=\""<<roots[i]->id<<"\">\n";
		out<<"<data key=\"M\">"<<i<<"</data>\n";
		roots[i]->tg.printGraphML(out,groupSizPerType);			
		out<<"</node>\n";
	}
	
	//output edges:
	for(it= nodes.begin();it!=nodes.end();++it)
	{
		const node&n=**it;
		node * a = n.pa(),* b = n.pb();
		if(a)
			out<<"<edge source=\""<<n.id<<"\" target=\""<<a->id<<"\"/>\n";
		if(b)
			out<<"<edge source=\""<<n.id<<"\" target=\""<<b->id<<"\"/>\n";
	}

	for(int i=0;i<roots.size();++i)
	{
		const node&n=*roots[i];
		node * a = n.pa(),* b = n.pb();
		if(a)
			out<<"<edge source=\""<<n.id<<"\" target=\""<<a->id<<"\"/>\n";
		if(b)
			out<<"<edge source=\""<<n.id<<"\" target=\""<<b->id<<"\"/>\n";
	}
	

	out<<"</graph>\n";
}



#endif

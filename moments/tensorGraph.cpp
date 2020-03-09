/*
	representation of a "tensor graph" from the symmetry information of a set of invarinat moments
	(of which each is a total contraction of tensor products)

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

#include<algorithm>
#define _USE_MATH_DEFINES
#include<math.h>
#include "tensWithSymm.h"
#include "binaryTreesVisitor.h"
#include "tensorGraph.h"
#include "tensorGraphSplitTree.h"
#include <time.h>
#include<fstream>
#include<string>

using namespace std;

//implementation-specific structs
struct tensorGraph::impl
{
	//for union_find
	struct unionFinder;

	//for ps output
	struct tensInfo;

	struct polynTensProd;

};



tensorGraph::tensorGraph(const tensSymmInfo&tws,const vpr& pairs)
{
	assert(sizeof(uint4) >= 4*sizeof(uchar));

	const vector<int>& tg=tws.tgroups();

	nodelabels.resize(tg.size());

	//first entry in pair:
	// node id per index in edges
	//=source tensor id for index
	//second entry: endpoint label
	//=number of symm. index group in tensor
	vector< pair<int,int> > buf;	
	vector<int>::const_iterator it = tws.groups().begin();
	
	for(int i=0;i<tg.size();++i)
	{
		//set nodelables to the number of indices in the tensor
		nodelabels[i] = 0;
		for(int j=0;j<tg[i];++j,++it)
		{
			nodelabels[i] += *it;
			pair<int,int> pr(i,j);
			buf.insert(buf.end(),*it,pr);				
		}
	}

	edges.resize(pairs.size());
	vector<edge>::iterator eit=edges.begin();
	vpr::const_iterator pit=pairs.begin();

	//fill edges
	for(;eit!=edges.end();++eit,++pit)
	{
		//for safety:
		eit->encompass=0;
		
		pair<int,int> & a = buf[pit->a];
		eit->d.a= a.first;
		eit->d.u= a.second;

		pair<int,int> & b = buf[pit->b];
		eit->d.b= b.first;
		eit->d.v= b.second;	
	}

}


void tensorGraph::permute(const vector<int>&perm)
{
	assert(sizeof(uint4) >= 4*sizeof(uchar));

	vector<int> buf(nodelabels.size());

	//permute the node labels
	for(int i=0;i<buf.size();++i)
		buf[perm[i]] = nodelabels[i];
	buf.swap(nodelabels);

	//correct the indices in the edges
	for(int i=0;i<edges.size();++i)
	{
		edge::data &d=edges[i].d;
		d.a=(uchar)perm[d.a];
		d.b=(uchar)perm[d.b];
		d.normalize();
	}

	std::sort(edges.begin(),edges.end());

}


void tensorGraph::allowedPermutationVisitor::recurse()
{
			if(depth+2<l.size())
			{
				//first, no perm in this level:
				++depth;
				recurse();
				--depth;
				//then: permutation with others having the same label
				for(int i=depth+1;i<l.size();++i)
					if(l[i]==l[depth])
					{
						xch(depth,i);
						++depth;
						recurse();
						--depth;
						//undo permutation
						xch(depth,i);				
					}			
			}
			else
			{
				accept();				
				for(int i=depth+1;i<l.size();++i)
					if(l[i]==l[depth])
					{
						xch(depth,i);
						accept();
						xch(depth,i);
					}
			}
}

/**
*create a permutation which corresponds to
*the permutations a and b executed one after another
*/
void tensorGraph::concatperm
(	vector<int>&ret,const vector<int>&a,const vector<int>&b )
{
	assert(a.size()==b.size());
	ret.resize(a.size());
	for(size_t i=0;i<a.size();++i)
		ret[i] = b[a[i]];
}



void tensorGraph::normalize( vector<int> * permut )
{

	int n=nnodes();
	
	//first, sort the nodes by their labels.
	vector<pair<int,int> >sorter(nnodes());
	for(int i=0;i<n;++i)
	{
		sorter[i].first = nodelabels[i];
		sorter[i].second = i;
	}

	sort(sorter.begin(),sorter.end());
	
	vector<int>perm(nnodes());
	for(int i=0;i<n;++i)
	{
		perm[ sorter[i].second ] = i;
	}


	permute(perm);
	
	//then, find the minimal graph under the equivalence class

	struct normalizeVisitor:public allowedPermutationVisitor
	{			
		tensorGraph mingraph;
		tensorGraph orig;
		tensorGraph tmp;
		vector<int> optimp;

		void accept()
		{
			tmp = orig;
			tmp.permute(p);
			if(tmp<mingraph)
			{
				mingraph.swap(tmp);
				optimp=p;
			}
		}

		normalizeVisitor(const tensorGraph & tg)
		:allowedPermutationVisitor(tg.nodelabels),
		orig(tg),mingraph(tg)
		{
			optimp=p;
		}
	};

	normalizeVisitor nv(*this);
	nv.visit();
	*this = nv.mingraph;
	
	if(permut)
	{
		concatperm(*permut,perm,nv.optimp);
	}

}

void tensorGraph::subgraph
(const tensorGraph& tg,const vector<int>& nis)
{

	int n=nis.size();
	
	nodelabels.resize(n);	

	//maps original node ids to new ones
	vector<int> nodeidmap(tg.nnodes(),-1);
	int nnodes=0;
	for(int i=0;i<n;++i)
	{
		nodeidmap[nis[i]]=i;
		nodelabels[i] = tg.nodelabels[nis[i]];
	}

	edges.clear();

	//collect edges in subgraph
	edge buf;
	for(int i=0;i<tg.edges.size();++i)
	{
		buf = tg.edges[i];
		int a = nodeidmap[buf.d.a];
		int b = nodeidmap[buf.d.b];
		if(a>=0 && b>=0)
		{
				buf.d.a=(uchar)a;
				buf.d.b=(uchar)b;
				buf.d.normalize();
				edges.push_back(buf);			
		}
	}
	//sort edges
	sort(edges.begin(),edges.end());
}


void tensorGraph::getConnMatrix(vector<unsigned> &conn) const
{
	conn.clear();
	conn.resize(nnodes(),0);
	
	for(int i=0;i<edges.size();++i)
	{
		edge buf = edges[i];
		conn[buf.d.a]|= 1<< buf.d.b;
		conn[buf.d.b]|= 1<< buf.d.a;	
	}
}

struct isnull
{
	inline bool operator()(unsigned x)const
	{return x==0;}
};



/**
implemented similar to
http://en.wikipedia.org/wiki/Union_find
still has te be debugged.
*/
struct tensorGraph::impl::unionFinder{
	
	struct nd{
		nd * parent;
		int rank;
	};

	nd nodes[tensorGraph::maxNodes];
	nd * buf[tensorGraph::maxNodes];
	
	int nnodes;
	
	//number of connected components
	int ncomps;

	/* find root element, update parent pointers on the way back*/
	nd* find(nd*x)
	{
		assert(x-nodes < nnodes );
		if(x->parent==x)
			return x;

		nd ** it = buf;
		
		do{
			//save pointer on node that has to be updated
			*it = x; ++it;
			//update x
			x = x->parent;
		} while(x->parent!=x);

		//update nodes		
		do{	
			--it;
			(**it).parent=x;
		}while(it!=buf);

		return x;
	}


	/**
	* init nodes as single connected components
	*/
	void init(int n)
	{

		assert(n<=tensorGraph::maxNodes);
		nnodes=n;
		ncomps=n;

		nd* x=nodes,*xend=nodes+n;
		for(;x!=xend;++x)
		{ 
			/** MakeSet in the wiki page*/
			x->parent=x;x->rank=0;
		}
	}	

	/* similar to union in the wiki page*/
	void join(tensorGraph::edge e)
	{
			//a, b are the roots of the nodes the edge connects
			nd * a = find(nodes+e.d.a);
			nd * b = find(nodes+e.d.b);
			if(a!=b)
			{
				--ncomps;
				if(a->rank > b->rank)
				{
					b->parent=a;
				}
				else if(a->rank < b->rank)
				{
					a->parent=b;
				}
				else
				{
					b->parent=a;
					a->rank++;
				}
			}	
	}


	/** join all nodes connected by an edge in e*/
	inline void join(const std::vector<tensorGraph::edge> &e)
	{
		std::vector<tensorGraph::edge>::const_iterator it;

		for(it=e.begin();it!=e.end();++it)
		{
			join(*it);
		}
	}

	inline void join(const std::vector<tensorGraph::edge> &e, int ignore)
	{
		std::vector<tensorGraph::edge>::const_iterator it,itign;
		itign=e.begin()+ignore;
		for(it=e.begin();it!=itign;++it)
			join(*it);

		for(++it;it!=e.end();++it)
			join(*it);	
	}


	inline unsigned notConnectedWithFirst(int n,
		const std::vector<tensorGraph::edge> &e,int ign)
	{
		unsigned ret=0;
		init(n);
		join(e,ign);
		nd * first = find(nodes+0);	
		for(int i=1;i<n;++i)
		{
			if( find( nodes+i ) !=first )
				ret|=1<<i;
		}
		return ret;
	}


};



void testUnionFinder()
{

	vector<tensorGraph::edge> e(6);

	e[0].d.a=0;
	e[0].d.b=1;
	e[1].d.a=0;
	e[1].d.b=2;
	e[2].d.a=4;
	e[2].d.b=5;
	e[3].d.a=5;
	e[3].d.b=6;
	e[4].d.a=6;
	e[4].d.b=7;
	e[5].d.a=0;
	e[5].d.b=5;
	
	/*
	*all in one group; except 3
	*/

	tensorGraph::impl::unionFinder f;

	f.init(8);

	cout<<"joining"<<endl;
	f.join(e);

	cout<<"groups:"<<endl;

	for(int i=0;i<f.nnodes;++i)
		cout<< (f.find(f.nodes+i)-f.nodes)<<endl;



}
	

struct 
	tensorGraph::impl::
	polynTensProd : public tens<polyn>
{

	const std::vector<tensor<polyn> > &t;
	const vector<int> & l;


	static int sumord(const std::vector<tensor<polyn> > &tensPerNodeLabel_,
						const vector<int> & nodelabels_)
	{
		int ret=0;
		for(int i=0;i<nodelabels_.size();++i)
			ret+=tensPerNodeLabel_[nodelabels_[i]].ord();
		return ret;		
	}


	polynTensProd(
		const std::vector<tensor<polyn> > &tensPerNodeLabel_,
		const vector<int> & nodelabels_
		)
		:
		t(tensPerNodeLabel_),
		l(nodelabels_),
		tens<polyn>(
			sumord(tensPerNodeLabel_,nodelabels_),
			tensPerNodeLabel_[nodelabels_[0]].dim()
			)
	{
	}


	virtual ~polynTensProd()
	{
	}

  polyn operator[](const int * inds) const
  {
		const tens<polyn> &t0=t[l[0]];
		polyn ret =  t0[inds];
		inds+=t0.ord();
		for(int i=1;i<l.size();++i)
		{
			const tens<polyn> &ti=t[l[i]];
			ret*=ti[inds];
			inds+=ti.ord();
		}
		return ret;
  }
};




/**create a polynomial corresponding to 
*the tensor Graph.
*\pre: the tensorGraph is complete(i.e. no open indices)
*/
void tensorGraph::createPolynTensor(
  tensor<polyn>& p,
	const std::vector<tensor<polyn> > &tensPerNodeLabel,
	const std::vector<std::vector<int> > & groupsPerNodeLabel) const
{
	vpr pairs;
	edges2indpairs(pairs,groupsPerNodeLabel);	
	impl::polynTensProd ptp(tensPerNodeLabel,nodelabels);
	
	//a basic check for completeness
	assert(ptp.ord() == pairs.size()*2);
	partialFold<polyn> pf(ptp,pairs);
	p=pf;
}

//convert the graph's edges to index pairs
void tensorGraph::edges2indpairs
( vpr&pairs,const std::vector<std::vector<int> > & groupsPerNodeLabel ) const
{
	//init offset vector
	vector<vector<int> > offsets( nnodes() );
	int off=0;
	for(int i=0;i<nnodes();++i){
		const vector<int>& gnl=groupsPerNodeLabel[nodelabels[i]];
		vector<int> &offs = offsets[i];
		
		offs.resize(gnl.size());
		for(int j=0;j<gnl.size();++j)
		{
			offs[j]=off;
			off+=gnl[j];
		}
	}

	//init pairs

	pairs.resize(edges.size());

	for(int i=0;i<edges.size();++i)
	{
		edge::data e=edges[i].d;
		pr & p = pairs[i];
		//note the post-fix ++, so the value before increment is assigned
		p.a = offsets[e.a][e.u]++;
		p.b = offsets[e.b][e.v]++;	
	}

}



/* 
 * gives all components not connected to first node
 * if cutedge is deleted
*/
unsigned tensorGraph::notConnectedWithFirst(int cutedge) const
{
	int gid[tensorGraph::maxNodes];

	int n=nnodes();
	for(int i=0;i<n;++i)
	{
		gid[i]=i;
	}

	//algo is of order n^2,
	//but very simple. because n is small, that's acceptable.
	for(int i=0;i<edges.size();++i)if(i!=cutedge)
	{
		edge buf = edges[i];
		int ga=gid[buf.d.a],gb=gid[buf.d.b];
		if(ga!=gb)
			for(int j=0;j<n;++j)
				if( gid[j]==gb )
					gid[j]=ga;
	}

	unsigned ret=0;
	for(int i=1;i<n;++i)
		ret|= (gid[i]!=gid[0])<<i;
	return ret;
}

bool tensorGraph::connectedComps(vector<unsigned> *comps) const
{

	int gid[tensorGraph::maxNodes];

	int n=nnodes();
	for(int i=0;i<n;++i)
	{
		gid[i]=i;
	}

	for(int i=0;i<edges.size();++i)
	{
		edge buf = edges[i];
		int ga=gid[buf.d.a],gb=gid[buf.d.b];
		if(ga!=gb)
			for(int j=0;j<n;++j)
				if( gid[j]==gb )
					gid[j]=ga;
	}

	if(!comps)
	{
		for(int i=1;i<n;++i)
			if(gid[i]!=gid[0])
				return false;
		return true;
	}
	else
	{

		comps->resize(n);
		for(int i=0;i<n;++i)
			(*comps)[gid[i]]|=1<<i;

		comps->resize(
			std::remove_if(comps->begin(),comps->end(),isnull())
			-comps->begin());
	
		return comps->size()==1;
	}
}

void tensorGraph::connectedComps(vector<tensorGraph> & ret)const
{
	vector<unsigned> groups;
	connectedComps(&groups);
	connectedComps(ret,groups);
}

void tensorGraph::connectedComps(vector<tensorGraph> & ret,const vector<unsigned>&groups)const
{
		//only one group
	if(groups.size()==1)
	{ ret.resize(1);ret[0]=*this; return;}

	//multiple groups:
	ret.resize(groups.size());
	vector<int> buf;
	for(int i=0;i<groups.size();++i)
	{
		buf.clear();
		unsigned grp=groups[i];
		for(int j=0;grp;grp>>=1,++j)
			buf.push_back(j);

		ret[i].subgraph(*this,buf);	
	}

}



bool tensorGraph::connectedComps(const vector<unsigned>&conn,
																 vector<unsigned> * ret)
{

	//1 bit for every node stating if it has been assigned already
	unsigned used=0;

	if(ret)
		ret->clear();

	int nodestack[tensorGraph::maxNodes];
	int * const stackStart = nodestack + tensorGraph::maxNodes;
	int * s = stackStart;

	//loop through all nodes
	for(int i=0; i<conn.size(); ++i)
	{
		//if node wasn't used already, start a new group
		if(((used>>i)&1)==0)
		{
			//start a new group
			unsigned newgrp=0;
			newgrp |= 1<<i;
			used |= 1<<i;

			//now collect all nodes in the same conn. comp.
			*(--s)=i;
			while(s!=stackStart)
			{
				int next = *s; ++s;

				//all nodes not already in use and connected to next:
				unsigned con_new = conn[next] & ~used;
				used |= con_new;
				newgrp |= con_new;
				
				for(int i=0;con_new;++i,con_new>>=1)
				{
					if(con_new&1)
					{
						*(--s)=i;
					}											
				}
			}		

			if(ret)
			{
				ret->push_back(newgrp);
			}
			else
			{
				//returns true if newgroup contains all nodes
				return newgrp == 	( 1 << conn.size() ) - 1;
			}
		}
	}

	return ret->size()==1;
}


void tensorGraph::split::init
(tensorGraph&tg1,
 const tensorGraph& orig, 
 int edgeid)
{
	tg1=orig;
	edgessplit.clear();
	edgessplit.resize(1,tg1.edges[edgeid]);					
	tg1.edges.erase(tg1.edges.begin()+edgeid);
	nodeidmap.clear();
}




void tensorGraph::split::init
(tensorGraph&tg1,
 tensorGraph&tg2,
 const tensorGraph& orig, 
 unsigned int nodein1)
{
	tg1.clear();
	tg2.clear();

	nodeidmap.resize(orig.nnodes());

	for(int i=0;i<orig.nnodes();++i)
	{
		if(nodeinfirst(i,nodein1))
		{
			nodeidmap[i] = tg1.nodelabels.size();
			tg1.nodelabels.push_back(orig.nodelabels[i]);
		}
		else
		{
			nodeidmap[i] = tg2.nodelabels.size();
			tg2.nodelabels.push_back(orig.nodelabels[i]);
		}
	}

	//split edges
	tensorGraph::edge buf;

	//the list of original edges
	edgessplit.clear();
	for(int i=0;i<orig.edges.size();++i)
	{
		buf=orig.edges[i];
		bool ba = nodeinfirst(buf.d.a,nodein1);
		bool bb = nodeinfirst(buf.d.b,nodein1);

		//if edge is not split:
		if( ba == bb) 
		{
			buf.d.a = (unsigned char) nodeidmap[buf.d.a];
			buf.d.b = (unsigned char) nodeidmap[buf.d.b];			
			buf.d.normalize();

			if(ba)
				tg1.edges.push_back(buf);
			else
				tg2.edges.push_back(buf);
		}
		else
			edgessplit.push_back(buf);		
	}		

	sort(tg1.edges.begin(),tg1.edges.end());
	sort(tg2.edges.begin(),tg2.edges.end());

}






struct tensorGraph::impl::tensInfo
{
	
	struct indgrpinfo
	{
		//endpoint coordinates=center of bubbles
		float ex,ey;

		//# internal pairs
		int npairs;

		//# indices not paired
		int freeinds;

		//total number of indices
		int totinds;
	};

	vector<indgrpinfo> ep;

	//sin,cos of angle
	float ca,sa;
	
	//total ord of tensor
	int ord;

};




//get the number of indices for every group in every node
void tensorGraph::getFreeIndsPerGroup
( vector<vector<int> >	&	freeIndsPerGroup,
			const vector<vector<int> > & groupsPerNodeLabel) const
{
	freeIndsPerGroup.resize(nodelabels.size());
	for(int i=0;i<nodelabels.size();++i)	
		freeIndsPerGroup[i]
			=groupsPerNodeLabel[nodelabels[i]];


	for(int i=0;i<edges.size();++i)
	{
		edge::data e = edges[i].d;
		freeIndsPerGroup[e.a][e.u]--;
		freeIndsPerGroup[e.b][e.v]--;		
	}
}


void tensorGraph::printGraphML(
		std::ostream&out,
		const std::vector< std::vector<int>  >& groupSizPerType,
		const char*name
	)
{
	std::vector< std::vector<int>  > freei;
	getFreeIndsPerGroup(freei,groupSizPerType);

	if(name)
		out<<"<node id=\""<<name<<"\"/>\n";
	else
	{
		out<<" <graph id=\"G\" edgedefault=\"undirected\">\n";
		out<<"  <key id=\"l\" for=\"node\" "
			"attr.name=\"tensor_id\" attr.type=\"int\" />\n";
	}

	string prefix;

	if(name)
		prefix=string(name)+"_";

	for(int i=0;i<nodelabels.size();++i)
	{
		int l=nodelabels[i];
		const vector<int> &g=groupSizPerType[l];
		out<<"  <node id=\""<<prefix<<i<<"\">";
		out<<"<data key=\"l\"> "<<l<<" </data>";				
		for(int j=0;j<g.size();++j)
		{
			out<<" <port name=\""<<j<<"\"/> ";
		}
		out<<"</node>\n";
		
		if(name)
		{
			//edge to supernode
			out<<
				"<edge "
				"source=\""<<prefix<<i<<"\" "
				"target=\""<<name<<"\""
				"/>"; 		
		}
		vector<int>&f=freei[i];
		for(int j=0;j<f.size();++j)
		{
			for(int k=0;k<f[j];++k)
			{
				out<<"  <node id=\""<<prefix<<i<<'_'<<j<<'_'<<k<<"\"/>";
				out<<" <edge "
					"source=\""<<prefix<<i<<"\"" 
					" sourceport=\""<<j<<"\"" 
					" target=\""<<prefix<<i<<'_'<<j<<'_'<<k<<"\""
					"/>\n";
			}
		}			
		out<<"\n";
	}


	for(int i=0;i<edges.size();++i)
	{
		edge::data d = edges[i].d;
		out<<
			"  <edge "
			" source=\""<<prefix<<int(d.a)<<"\""
			" sourceport=\""<<int(d.u)<<"\" "
			" target=\""<<prefix<<int(d.b)<<"\""
			" targetport=\""<<int(d.v)<<"\" "
			"/>\n";	
	}



	if(!name)
		out<<" </graph>\n";
}


float tensorGraph::printPS(ostream&out,const vector< vector<int>  >& groupSizPerType,float bd,bool epsheader) const
{
	

	vector< impl::tensInfo  > ti(nnodes());

	//collect stats about the tensors:
	int sumord=0;
	int maxord=0;
	for(int i=0;i<nodelabels.size();++i)
	{
		const vector<int>& grp=groupSizPerType[nodelabels[i]];
		impl::tensInfo & t = ti[i];
		t.ep.resize(grp.size());
		t.ord = 0;
		for(int j=0;j<grp.size();++j)
		{ 
			impl::tensInfo::indgrpinfo &g =t.ep[j];
			if(maxord<grp[j])maxord=grp[j];
			g.totinds = grp[j]; 
			g.freeinds = g.totinds;
			g.npairs=0;
			t.ord += g.totinds;
		}
		sumord +=t.ord;	
	}


	float bggray=1;
	if( sumord == edges.size()*2 )
		bggray=0;


	float r;
	float arcr;
	
	if(nnodes()>1)
	{
		r= bd * (sumord + nnodes())/2/M_PI ;
		if( maxord * bd > r)
			r= maxord*bd;
		arcr = r + maxord*bd;
	
	}
	else
	{
		r =  bd * sumord /2;
		arcr = r+bd/2;	
	}

	if(epsheader)
		out<<
		"%!PS-Adobe-2.0 EPSF-2.0\n"
		"%%BoundingBox: -1 -1 "<<int(arcr*2+2.5)<<' '<<int(arcr*2+2.5)<<'\n';

	if(nnodes()>1)
	{
		out<< arcr <<" dup translate\n" ;
		out<<" 0 0 "<< arcr<<" 0 360 arc "<<bggray<<" setgray fill\n";
		out<<" 0 0 "<< arcr<<" 0 360 arc 0 setgray stroke\n";

	}
	else
	{
		out<< arcr-r <<' '<< arcr <<" translate\n";
		out<< r <<" 0 "<< arcr <<" 0 360 arc "<<bggray<<" setgray fill\n";
		out<< r <<" 0 "<< arcr <<" 0 360 arc 0 setgray stroke\n";
	}


	float ang0 = 0;

	//calculate tensor-bubble centers
	for(int i=0;i<nodelabels.size();++i)
	{
		impl::tensInfo & t = ti[i];

		t.ca = cos(ang0); t.sa = sin(ang0);
		
		//center of a tensor-bubble- group
		float sx = t.ca*r, sy = t.sa*r;
		
		
		float grp2= 0.5*t.ord;
		float dx = -bd * t.sa, dy= bd * t.ca ;
		float ex= sx - dx * grp2;
		float ey= sy - dy * grp2;
	

		for(int j=0;j<t.ep.size();++j)
		{			
			impl::tensInfo::indgrpinfo &g= t.ep[j];
			//center of a tensor-bubble
			g.ex= ex +dx*g.totinds/2   ; g.ey= ey + dy*g.totinds/2;
			ex += dx*g.totinds;
			ey += dy*g.totinds;
		}


		if( i+1<nodelabels.size() )
			ang0 += (.5*t.ord + .5*ti[i+1].ord +1) *2*M_PI / (sumord+nnodes());
	}


	//draw arcs
	for(int i=0;i<edges.size();++i)
	{

		int mult=1;
		while( i+1 < edges.size() && edges[i] == edges[i+1]) 
			++i,++mult;

		edge::data e = edges[i].d;

		impl::tensInfo
			&a=ti[e.a],
			&b=ti[e.b];

		impl::tensInfo::indgrpinfo
			&au=a.ep[e.u],
			&bv=b.ep[e.v]	;

		au.freeinds-=mult;			
		bv.freeinds-=mult;			


		if(e.a==e.b&&e.u==e.v)
		{
			au.npairs += mult;			
		}		
		else{

			for(int i=0;i<mult;++i)
			{
				//the displacement perpendicular to the radius
				float p =  (i - (mult-1)*0.5 );					
				//the radial displacement of bezier control points
				float d =  0.2;
				if(e.a==e.b)
				{
					d /= 1+i ;
					p=0;
				}
					out<< au.ex <<' '<< au.ey <<" moveto\n";

					out
						<<au.ex*d - p* au.ey <<' '
						<<au.ey*d + p* au.ex <<' '
						<<bv.ex*d - p* au.ey <<' '
						<<bv.ey*d + p* au.ex <<' '
						;
					out<<bv.ex<<' '<<bv.ey<<" curveto\n";
			}
		}
	}
	out<< "1 setlinewidth "<<1-bggray<<" setcolor stroke \n";	

	//draw bubbles and self arcs and free arcs
	for(int i=0;i<ti.size();++i)
	{
		/*
		vector<pair<float,float> > &ep=endpoints[i];
		const vector<int> &fgrp = freepergroup[i];
		const vector<int> &grp = groupSizPerType[nodelabels[i]];
		vector<int> &pp = pairspergroup[i];
		*/
		impl::tensInfo &t = ti[i];

		float dx = -bd * t.sa, dy= bd * t.ca ;

		for(int j=0;j<t.ep.size();++j)
		{
			impl::tensInfo::indgrpinfo &e=t.ep[j];			

			//draw self arcs and free arcs			
			int n= e.freeinds + 2*e.npairs;

			float ex = e.ex - dx * (n-1)/2;
			float ey = e.ey - dy * (n-1)/2; 

			float fx=dx*e.totinds,fy=dy*e.totinds;

			//draw open arcs:
			int k=0;
			for(;k<e.freeinds;++k)
			{				
				out<<ex<<' '<<ey<<" moveto\n";
				out<<fy<<' '<<-fx<<" rlineto\n";				
				ex+=dx;ey+=dy;
			}
			
			fx*=1.2;fy*=1.2;

			//draw self-arcs
			for(;k<n;k+=2)
			{
				out<<ex<<' '<<ey<<" moveto\n";
				out<<ex+fy<<' '<<ey-fx<<"\n";				
				ex+=dx;ey+=dy;
				out<<ex+fy<<' '<<ey-fx<<"\n";				
				out<<ex<<' '<<ey<<" curveto\n";						
				ex+=dx;ey+=dy;
			}
			

			out<<1-bggray<<" setgray stroke\n";

			//draw bubbles
			out<<e.ex<<' '<<e.ey<<' '<<bd/2 * e.totinds <<" 0 360 arc\n";
			
			int mj=t.ep.size()-1;
			if(mj>0)
			{
				out<< 0.5 + 0.2 * (j*2.0/mj -1) <<" setgray fill\n";			
				//draw outline in color of brightest item

				if(j!=mj)
				{
					out<<e.ex<<' '<<e.ey<<' '<<bd/2 * e.totinds <<" 0 360 arc\n";
					out<< 0.5 + 0.2 <<" setgray 1 setlinewidth stroke\n";
				}
			}
			else
				out<< "0.5 setgray fill\n";			


		}	
	}

	if(nnodes()>1)
		out<<arcr<<" dup neg translate\n" ;
	else
		out<<arcr+r<<" "<<-arcr<<" translate\n";

	return arcr*2 ; 

}





ostream & operator <<(ostream&o,const tensorGraph &tg)
{
	o<<"nodelabels:";
	for(int i=0;i<tg.nodelabels.size();++i)
		o<<' '<<tg.nodelabels[i];
	o<<"\nedges:";

	const vector<tensorGraph::edge> & e = tg.edges;
	for(int i=0;i<e.size();++i)
	{
		tensorGraph::edge::data d=e[i].d;
		o<<'('<<int(d.a)<<' '<<int(d.u)
			<<','<<int(d.b)<<' '<<int(d.v)<<')';
	}	

	return o;

}

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

#ifndef TENSOR_GRAPH_H
#define TENSOR_GRAPH_H
#include<vector>
#include "vvpr.h"
#include "refp.h"
#include "binaryTreesVisitor.h"
#include "tensor.h"
#include "polynom.h"


struct tensSymmInfo;





struct tensorGraph:public refCountable
{
	
	static const int maxNodes = sizeof(unsigned)*8;

	inline tensorGraph(){}
	
	tensorGraph(const tensSymmInfo&tws,const vpr& edges);
	
	int nnodes()const {return (int)nodelabels.size();}

	typedef unsigned char uchar;
	typedef unsigned int uint4;

	//the edges

	union edge {
		struct data{

			//the node ids
			uchar a,b;
			
			//endpoint labels,
			//denotin to which index group in the tensor 
			//this edge connects
			uchar u,v;
			
			inline void normalize()
			{
				if(a>b)
				{
					uchar xch=a;a=b;b=xch;
					xch=u;u=v;v=xch;
				}
			}

		}d;

		uint4 encompass;

		inline bool operator == (const edge&x) const
		{
			return encompass==x.encompass;
		}
	
		inline bool operator < (const edge&x) const
		{
			return encompass<x.encompass;
		}

	};


	//invoke f for every permutation which is still allowed
	//(does not interchange labels)
	struct allowedPermutationVisitor{
		
		//labels
		const std::vector<int>& l;

		//permutation
		std::vector<int> p;
		
    //recursion depth
		int depth;

		inline void xch(int i,int j)
		{
			int xx=p[i];p[i]=p[j];p[j]=xx;
		}

		inline allowedPermutationVisitor(const std::vector<int> &l)
			:l(l),p(l.size())
		{
			//set to unit permutation
			for(int i=0;i<p.size();++i)
				p[i]=i;			
		}

		void visit()
		{
			depth=0;
			recurse();
		}

		void recurse();

		//accept the current permutation
		virtual void accept()=0;

	};


	//the edges. They will always be sorted.
	std::vector<edge> edges;
	
	//for every node, the tensor id
	std::vector<int> nodelabels;

  //permute the nodes
	void permute(const std::vector<int>& perm);

	
  //normalize this graph;
	//optionally return the permutation used for it
	void normalize(std::vector<int>*perm=0);

	inline bool operator == (const tensorGraph&x) const
	{
		return edges==x.edges && nodelabels==x.nodelabels;
	}

	inline bool operator < (const tensorGraph&x) const
	{
		return 
			nnodes() < x.nnodes() ||
			nnodes() == x.nnodes() &&
			(
				edges.size()< x.edges.size() ||
				edges.size()== x.edges.size() &&
				(
					lexicographical_compare
					( edges.begin(),edges.end(),
						x.edges.begin(),x.edges.end() )	|| 
					edges==x.edges && 
					lexicographical_compare
					( 
						nodelabels.begin(),nodelabels.end(),
						x.nodelabels.begin(),x.nodelabels.end() 				
					)
				)
			);	
	}


	friend struct tensorGraphSplit;

	void subgraph(const tensorGraph& tg,const std::vector<int>&nodesInSubgraph);
	

	inline void clear()
	{
		edges.clear();nodelabels.clear();
	}

	inline void swap(tensorGraph&tg)
	{
		tg.edges.swap(edges);
		tg.nodelabels.swap(nodelabels);
	}
	
	//every entry of conn is a bitvector,
	//with one bit for every node 
	//this could be conneted with.
	void getConnMatrix(std::vector<unsigned> &conn)const;

	static void 
		concatperm
	(	std::vector<int>&ret,const std::vector<int>&a,const std::vector<int>&b );

	/*
	*
	* \param conn: the matrix with the connections as returned by
	* getConnMatrix, the row index being the vector index
	* and the column index being the bit index
	* \param ret: 
	* if not null, a pointer on a 
	* vector with one unsigned for every component found
	*\param onlycheck
	*\return true if all nodes are in one connected component
	*/
	static bool connectedComps(const std::vector<unsigned>&conn,
		std::vector<unsigned> * ret=0);
	
	/**
	* \param ret: pointer on a vector which contains one vector with
	* node indices for every connected component
	*\param withoutedge: assumes edge number withoutedge wasn't present
	*\return true if all nodes are in one connected component
	*/
	bool connectedComps(std::vector<unsigned>*comps=0) const;

	void connectedComps(std::vector<tensorGraph> & ret,const std::vector<unsigned>&groups)const;
	void connectedComps(std::vector<tensorGraph> & ret)const;

	/* 
	* gives all components not connected to first node
	* if edge cutedge is ignored
	*\param cutedge: number of the edge that is ignored
	*\return:
	* bit i of the return value indicates if the i'th component is not connected
	*/
	unsigned notConnectedWithFirst(int cutedge) const;


	struct split
	{
		std::vector<tensorGraph::edge> edgessplit;
		std::vector<int> nodeidmap;	

		static inline bool nodeinfirst(int i,unsigned nodein1) 
		{ return (nodein1>>i)&1; }

		void init(tensorGraph&tg1,tensorGraph&tg2,
						const tensorGraph&orig,						
						unsigned int nodein1);

		//split only one edge
		void init(tensorGraph&tg1,const tensorGraph&orig,	int edgeid);
	};

	//for implementation structs;
	struct impl;

	/**
	print graph as a postscript picture.
	\param groupSizPerType has for every node type the sizes of the index groups.
	\param bd is the basic diameter for a 1-index group
	\return : the size of the image in native units
	*/
	float printPS(
		std::ostream&out,
		const std::vector< std::vector<int>  >& groupSizPerType,
		float bd=10,
		bool epsheader=false
		) const;

	/*
	 output tensor Graph in graphML
	 (witout headers)
		\param groupSizPerType has for every node type the sizes of the index groups.
		\param name:
		 base name for node ids and
		 node id of supernode
	*/
	void printGraphML(
		std::ostream&out,
		const std::vector< std::vector<int>  >& groupSizPerType,
		const char* name=0
	);

	//get the number of indices for every group in every node
	void getFreeIndsPerGroup
		( std::vector<std::vector<int> >	& freeIndsPerGroup,
			const std::vector<std::vector<int> > & groupsPerNodeLabel) const;


	//convert the graph's edges to index pairs
	void edges2indpairs(vpr&pairs,
		const std::vector<std::vector<int> > & groupsPerNodeLabel) const;

	/**create a tensor of polynomials corresponding to 
	*the tensor Graph.
	*/
	void createPolynTensor(
	  tensor<polyn> & p,
		const std::vector<tensor<polyn> > &tensPerNodeLabel,
		const std::vector<std::vector<int> > & groupsPerNodeLabel) const;


};


std::ostream & operator <<(std::ostream&o,const tensorGraph &tg);


#endif

/*
	create "tensor graph split tree" which
	creates a network of tensor products and contractions  which describe the calculation of a set of moment invariants
	in an optimized manner (main method: findOptimalSplits)

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

#include<set>
#include<map>

#include "tensorGraph.h"

class binaryTreesVisitor;
class repeatTensor;

struct tensorGraphSplitTree
{	

	struct node;

	//a stack with free node objects to be re-used
	std::vector< node* > freenodes;


	struct lt_pnode{
		inline bool operator ()(const node*a,const node*b) const
		{
			return a->tg < b->tg;
		}
	};

	//all nodes which are refernced in node objects
	typedef	std::set< node*, lt_pnode > nodeset;
	
	nodeset nodes;


	//remove links to children and evtl. delete them
	void freechild(nodeset::iterator it);
	
	//get a new free node
	node* getNode()
	{
		if(freenodes.empty())
			return new node(*this);
	
		node*x=freenodes.back();
		freenodes.pop_back();
		return x;
	}

	nodeset::iterator insertchild(node * it,bool&nw);
	//a stripped split info	
	//used to re-create the optimal split configuration
	struct splittree{

		struct stnode{
			//the split
			unsigned split;
			//the child indices =indices in the nodes array,
			//or -1 if leaf
			int a,b;

		};

		//the nodes in left-right-root-order
		std::vector<stnode> nodes;

		int initRecurse(const node*root);
		void init(const node*root)
		{
			nodes.clear();
			nodes.reserve(root->tg.nnodes()-1);
			initRecurse(root);			
		}

		int init(const binaryTreesVisitor::node * t);


		void recurseRecreate(int nodeid,node*root) const;
		/** recreate the split given by this tree,
			* but don't set the split for nodes used also by other subtrees
			*/
		void recreate(node*root) const
		{		
			recurseRecreate(nodes.size()-1,root);
		}




	};



	struct node{
		//the number of parent nodes
		int nparents;

		//the instance of tgs (used for getNode() etc)
		tensorGraphSplitTree & tgs;

		tensorGraph tg;
		
		//the permutations that led to the normal forms of
		//the children of this tensor graph;
		//needed for construction of foldings later
		std::vector<int> perma,permb;

		//the split metainfo
		tensorGraph::split s;

		//the current split
		unsigned split;
		

		//optimal split tree yet found:
		splittree optimalsplit;

		//the effort for the current split
		size_t getCurrentEffort();

		//the iterators of the subnodes in the set.
		nodeset::iterator a,b;

		node(tensorGraphSplitTree&tgs) : tgs (tgs) 
		{ 
			nparents=0;
			a=tgs.nodes.end();
			b=tgs.nodes.end();
		}



		void setSplit(unsigned newsplit,bool&nwa,bool&nwb);
		


		size_t findOptimalSplit();
	};




	tensorGraphSplitTree()
	{
	}


	~tensorGraphSplitTree()
	{

		for(int i=0;i<freenodes.size();++i)
			delete freenodes[i];

		nodeset::iterator it;
		for(it=nodes.begin();it!=nodes.end();++it)
			delete *it;

	}


	//the roots of the different tensor graphs
	std::vector<nodeset::iterator> roots;

	//the order of the tensor that corresponds to a node type
	std::vector<int> orderpernodetype;
	
	//the dimension of the tensors
	int dim;


	/**
	*the main routine.
	*find an optimal split tree of tensor Graph tg
	*where the nodes of tg represent
	*tensors of orders  as stored in oderpernode_
	*and dimension dim_
	*/
	size_t findOptimalSplit( 
		const tensorGraph& tg, 
		const std::vector<int> &orderpernodetype_,
		int dim);

	/**
	*the main routine.
	*find the split tree of tensor Graphs tgs
	*where the nodes of tg represent
	*tensors of orders  as stored in orderpernodetype_
	*and dimension dim_
	* and re-use subgraphs if they occur multiple times
	*/
	size_t findOptimalSplits( 
		const std::vector<tensorGraph> & tgs, 
		const std::vector<int> &orderpernodetype_,
		int dim);


	//print the entries of nodes
	void printPS(std::ostream&out,const std::vector<std::vector<int> > &groupSizPerType) const;
	
	
	
};



/***
* a split tree that splits only one edge at a time.
*/

struct tensorGraphSplitTree2
{	

	struct node;

	//the id to assign to the next node created.
	int nextnodeid;


	//update node Ids from the orderin in nodeset.
	void updateNodeIds();

	//a stack with free node objects to be re-used
	std::vector< node* > freenodes;


	struct lt_pnode{
		inline bool operator ()(const node*a,const node*b) const
		{
			return a->tg < b->tg;
		}
	};

	//all nodes which are refernced in node objects
	typedef	std::set< node*, lt_pnode > nodeset;
	
	nodeset nodes;


	//remove links to children and evtl. delete them
	void freechild(nodeset::iterator it);
	
	//get a new free node
	node* getNode()
	{
		if(freenodes.empty())
			return new node(*this);
	
		node*x=freenodes.back();
		freenodes.pop_back();
		x->edgeid=-1;
		return x;
	}

	nodeset::iterator insertchild(node * it,bool&nw);
	//a stripped split info	
	//used to re-create the optimal split configuration
	struct splittree{

		struct stnode{
			
			//the edgesplit
			int edgeid;
			//the child indices =indices in the nodes array,
			//or -1 if none there
			int a,b;

		};

		//the nodes in left-right-root-order
		std::vector<stnode> nodes;

		int initRecurse(const node*root);
		void init(const node*root)
		{
			nodes.clear();
			nodes.reserve(root->tg.nnodes()-1);
			initRecurse(root);			
		}

		int init(const binaryTreesVisitor::node * t);


		void recurseRecreate(int nodeid,node*root) const;
		/** recreate the split given by this tree,
			* but don't set the split for nodes used also by other subtrees
			*/
		void recreate(node*root) const
		{		
			recurseRecreate(nodes.size()-1,root);
		}




	};

	typedef std::map<const node*,std::set<const node*> > tparentset;


	struct node{

		int nparents;

		//a number to reference the node;
		//is assigned from tgs.nextnodeid in constructor.
		int id;


		//the instance of tgs (used for getNode() etc)
		tensorGraphSplitTree2 & tgs;

		tensorGraph tg;

		tensorGraph::split s;
		
		//which nodes are on the b side of the split ?
		unsigned compsb;

		//the permutations that led to the normal forms of
		//the children of this tensor graph;
		//needed for construction of foldings later
		std::vector<int> perma,permb;

		//the current edge to split
		int edgeid;		

		//optimal split tree yet found:
		splittree optimalsplit;

		//the effort for the current split
		size_t getCurrentEffort();

		//the iterators of the subnodes in the set.
		nodeset::iterator a,b;


		node * pa() const
		{  
			if( a!=tgs.nodes.end()) return *a;
			else return 0;
		}

		node * pb() const
		{  
			if( b!=tgs.nodes.end()) return *b;
			else return 0;
		}

		node(tensorGraphSplitTree2&tgs) : tgs (tgs) 
		{ 
			id=tgs.nextnodeid; 
			++tgs.nextnodeid;

			nparents=0;
			edgeid=-1;
			a=tgs.nodes.end();
			b=tgs.nodes.end();
		}

		/**
		*the number of distinct children
		*( if both children pointers represent the same node,
		*1 is returned)
		*/
		int numDistinctChildren()const{

			return 
				(a!=tgs.nodes.end())+
				( a!=b && (b!=tgs.nodes.end() ) );
		}

		void addToParentSet( tparentset & parents)const;

		void updateParentSetOf( const node* child ,tparentset &parents) const;


		void setSplit(int newsplit,bool&nwa,bool&nwb);

		/*
		*try to join this node with it's child,
		*if the child has only one parent	
		*/
		void tryjoinwithchild();


		size_t findOptimalSplit();

		
		//help functions for createRepeatTensors
		repeatTensor * createRepeatTensor
		(const repeatTensor*a,const repeatTensor*b) const;

		repeatTensor * createRepeatTensor
		(const repeatTensor*a) const;

	};


	void reset();

	
	static void groupFreeInds(
		std::vector<int>&indsPerNode,
		const std::vector<std::vector<int> >&fia
		);

	/*
	*convert a permutation of nodes
	*to a permutation of tensor indices
	*\param indperm:the index permuation
	*\param perm: the node permutation
	*\param indsPerNode: number of indices per node
	*/
	static void 
		perm2indperm( 
			std::vector<int>&indperm,
			const std::vector<int> &perm,
			const std::vector<int > &indsPerNode
			);

	tensorGraphSplitTree2(
		const std::vector< std::vector<int> > &groupSizPerType_,
		int dim_)
		:groupSizPerType(groupSizPerType_),dim(dim_)
	{
		nextnodeid=0;
		groupFreeInds(orderpernodetype,groupSizPerType);
	}


	~tensorGraphSplitTree2()
	{

		for(int i=0;i<freenodes.size();++i)
			delete freenodes[i];

		for(int i=0;i<roots.size();++i)
			delete roots[i];

		nodeset::iterator it;
		for(it=nodes.begin();it!=nodes.end();++it)
			delete *it;
	}


	//the roots of the different tensor graphs
	std::vector<node*> roots;

	//the order of the tensor that corresponds to a node type
	std::vector<int> orderpernodetype;

	//for every node type, the symmetric index group sizes
	std::vector< std::vector<int> > groupSizPerType;
	
	//the dimension of the tensors
	int dim;


	//to do progress-reports:
	int currentiter;
	int maxiter;

	/**
	*the main routine.
	*find an optimal split tree of tensor Graph tg
	*where the nodes of tg represent
	*tensors of orders  as stored in oderpernode_
	*and dimension dim_
	* uses also already present nodes.
	*/
	size_t findOptimalSplit( 	const tensorGraph& tg );

	/**
	*the main routine.
	*find the split tree of tensor Graphs tgs
	*where the nodes of tg represent
	*tensors of orders  as stored in orderpernodetype_
	*and dimension dim_
	* and re-use subgraphs if they occur multiple times
	*/
	size_t findOptimalSplits( const std::vector<tensorGraph> & tgs);
	
	
	
	//helper struct used in printps
	struct npinfo;

	//print the entries of nodes
	void printPS(std::ostream&out) const;
	
	//get the set of parents for every node
	void getParentSet(tparentset & parents) const;
	
	/**for every level, a map with nodes ids
	*1st level = nodes with no children,
	*2nd level = nodes with no children if 1st level is removed
	*etc.
	*/
	void createGraphLevelsBottomup(std::vector< std::set< const node* > >& levels	) const;

	/**for every level, a map with nodes ids
	*1st level = nodes with no parents =roots
	*2nd level = nodes with no parents if 1st level-nodes are removed
	*3rd level = nodes with no parents if 2nd level-nodes are removed
	*etc.
	*/
	void createGraphLevelsTopdown(std::vector< std::set< const node* > >& levels	) const;




	/**
	*creates a tree of mutually interconnected repeatTensors
	*\param nodes: for every node label, the corresponding repeatTensor
	*\param intermediate: repeat tensors in the middle
	*\param roots: for every root node, one repeat tensor
	*/
	void createRepeatTensors
		( 
		const std::vector<const repeatTensor* >& source, 
		std::vector<repeatTensor*> & intermediate,
		std::vector<repeatTensor*> & moments
		) const
		; 


	/**
	* join nodes which have only one parent
	* and whose parent has only one child 
	* with their parent.
	* ATTENTION:
	* this produces de-normalized nodes
	* in the nodes set
	*/
	void joinSingleParentNodes();


	/**
	*create a C-function out of the graph.
	*with A[...] as inputs and M[...] as outputs,
	* and buffer varaibles denoted woth lowercase names
	*/
	void printForC(std::ostream&out,const std::vector<const repeatTensor*> &source) const;


	void printGraphML(std::ostream&out) const;
	
};





/*
	header of binaryTreesVisitor visiting a rotated version of a binary tree.


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


#ifndef BINARY_TREES_VISITOR_H
#define BINARY_TREES_VISITOR_H
/**
* visits all possible binary trees
* with a fixed number of distinguishable
* leaves, leaving out the rotated instances.
*/
class binaryTreesVisitor
{

public:
  struct node
  {
    //bitvector with nodes contained
    unsigned nodescontained;
    
    //pointer on subnodes, if not null,
		node *a,*b;
		
		//else: these are the leaf ids of the children
    char cha,chb;
    
    //the level in the tree, level 0 = one above the leaves 
    char level;
  };

  void visit(int numleaves);
  
  virtual void accept(const node*root)=0;

private:  

	static const int maxNodes=32;
  node nodes[maxNodes];

	node * nodePerLevel[maxNodes];
    
  node * nodeIter;
  
  int numLeavesLeft;
  int level;
  int n;
    
  unsigned leavesUsed;
  
  void firstlevel();
  
  void upperlevel();
  
  void pairs1st(int minid,int nnodes);

  void pairsUpper(
    int minid,
    //#nodes to create in this level
    int nnodesleaf,
		int nnodesnoleaf,
    //unused subnodes
    unsigned freenodes
    );
  
};



#endif

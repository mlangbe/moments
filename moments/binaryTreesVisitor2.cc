/*
	Implementation part of binaryTreesVisitor2 visiting a rotated version of a binary tree.

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

#include "binaryTreesVisitor2.h"
#include<cassert>
  
  void binaryTreesVisitor2::visit(int numtypes_,const char*numPerType_)  {
		numtypes_=numtypes; numpertype = numPerType_;
		n=0;
		for(int i=0;i<numtypes;++i)
			n+=numpertype[i];
    level=0;nodeIter=nodes;
    leavesUsed=0;
    numLeavesLeft=n;
    firstlevel();
        
  }
  
	void binaryTreesVisitor2::firstlevel()
  {
    nodePerLevel[0]=nodeIter;
    for(int nnodes=1;nnodes <= n/2;++nnodes)
       pairs1st(0,nnodes);  
  }
  
  void binaryTreesVisitor2::upperlevel()
  {
    nodePerLevel[level]=nodeIter;
    int nperl= nodeIter - nodePerLevel[level-1];

    if(nperl==1 && numLeavesLeft==0)
      accept(nodeIter-1);
    else 
    {
			for(int nnodesnoleaf=0;nnodesnoleaf <= nperl/2;++nnodesnoleaf)
        pairsUpper(0,nperl-nnodesnoleaf*2,nnodesnoleaf,(1<<nperl)-1); 
    }
  }
  
  void binaryTreesVisitor2::pairs1st(int minid,int nnodes)
  {
    if(nnodes==0){
      ++level;
      
			if(level<n)
				upperlevel();
      
      --level;
      return;
    }
    for(int i=minid;i<n-1;++i)
      if(!((leavesUsed>>i)&1))
        for(int j=i+1;j<n;++j)
          if(!((leavesUsed>>j)&1))    
          {
            nodeIter->cha=i;
            nodeIter->chb=j;
						nodeIter->a=nodeIter->b=0;
            unsigned cnt = nodeIter->nodescontained = (1<<i)|(1<<j);
            ++nodeIter;
            
            numLeavesLeft-=2;
            leavesUsed |= cnt;
            pairs1st(  i+1, nnodes-1 );
            leavesUsed &= ~cnt;
            numLeavesLeft+=2;
            
            --nodeIter;          
          
          }
        
  }

  void binaryTreesVisitor2::pairsUpper(
    int minid,
    //#nodes to create in this level
    int nnodesleaf,
		int nnodesnoleaf,
    //unused subnodes
    unsigned freenodes
    )
  {
    
    if(nnodesleaf==0 && nnodesnoleaf==0){
      ++level;   
      if(level<n)upperlevel();
      --level;
      return;
    }
    
    //find afree subnode
    unsigned fn=freenodes>>minid;
    while( ((fn&1)==0) && fn!=0 ) ++ minid,fn>>=1 ;
		assert(fn);
    fn>>=1;
		
    //first, try to pair with other subnode:    
		if(nnodesnoleaf)
    for(int i=minid+1 ; fn!=0; ++i,fn>>=1) if(fn&1)
    {
    
            nodeIter->a = nodePerLevel[level-1] + minid;
            nodeIter->b = nodePerLevel[level-1] + i;
						nodeIter->cha=nodeIter->chb=0;
            nodeIter->nodescontained  = 
               (nodeIter->a)->nodescontained
               |
               (nodeIter->b)->nodescontained;
               
            ++nodeIter;
            pairsUpper(  minid+1, nnodesleaf ,nnodesnoleaf-1,freenodes & ~ ( (1<<minid)|(1<<i) ) );
            --nodeIter;
          
     }
     
     //then, try to pair whith leaves:
		 if(nnodesleaf)
     for(int i=0;i<n;++i) if(!((leavesUsed>>i)&1))    
          {
            unsigned cnt =  (1<<i);
            nodeIter->a=nodePerLevel[level-1] + minid;
            nodeIter->chb=i;
            nodeIter->cha=0;
            nodeIter->b=0;
						nodeIter->nodescontained = nodeIter->a->nodescontained | cnt;
            
            ++nodeIter;
            
            numLeavesLeft--;
            leavesUsed |= cnt;
            pairsUpper(  minid+1, nnodesleaf-1, nnodesnoleaf,freenodes & ~(1<<minid) );
            leavesUsed &= ~cnt;
            numLeavesLeft++;
            
            --nodeIter;          
          
          }
           
       
  }
  

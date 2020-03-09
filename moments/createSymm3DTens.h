/*
 	class to create polynimial tensor components for symmetric tensors
 	(having repetitions then)

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

#ifndef CREATE_SYMM_3DTENS_H
#define CREATE_SYMM_3DTENS_H

#include "pointCloud.h"
#include "tensWithSymm.h"
#include "polynom.h"

/*create tensor with entries represented the number from 
*a  tensorset  having maximum order of N
*(assigns indices created by calcIndex)
*/
struct createSymm3DTens
{

  tensor<polynom> x;

  static int nextid;


public:
  operator tensor<polynom>& (){
    return x;
  }

  createSymm3DTens ( int order,int N)
    :x(order,3)
  { 
    inditer ii(3,order);
    do{
			int nums[3]={0,0,0};

			const int* indices = *ii;
			for(int i=0;i<order;++i)
				nums[indices[i]]++;
							
      x[indices] = polynom(1,calcIndex(N,nums[0],nums[1],nums[2])) ;  

    }while(++ii!=0);
  }
  
};


/*create a tensor whose entries are polynomials with
*one tensor component id each.
*the tensor components ids are those given by calcindex
*and respect the symmetries of the tensor.
*/
struct createSymm3DTens2
{

  tensor<polynom> x;
	int fieldorder_;
	bool fieldsymm_;

public:
  operator tensor<polynom>& (){
    return x;
  }
	
	operator tensorWithSymm(){

		tensorWithSymm xx(x);
		if(x.ord()>fieldorder_)
			xx.grp.push_back(x.ord()-fieldorder_);
		
		if(fieldorder_)
		{
			if(fieldsymm_)
				xx.grp.push_back(fieldorder_);
			else for(int i=0;i<fieldorder_;++i)
			{
				xx.grp.push_back(1);
			}				
		}

		xx.tgrp.push_back(xx.grp.size());

    return xx;
  }

	/**
	* creates a tensor of a tensorset having two parts of indices:
	* one part consisting of "domain indices" which are
	* symmetric to each other (order), one part consisting of the  
	*indices of a tensorfield, with fieldsymm saying if it consists 
	* of symmetric tensors
	* \param N : the maximum domain order in the tensorset
	* \param order: the domain order of the tensor
	* \param fieldorder: the order of the field from which the indices are taken
	* \param fieldsymm: indicates if the field consists of symmetric tensors
	*(tensors which are symmetric to the interchangin of their indices)
	*/
  createSymm3DTens2 ( int order,int N,int fieldorder=0,bool fieldsymm=true,
		map<int,vector<int> > * decoder=0)
    :x(order+fieldorder,3),fieldorder_(fieldorder),fieldsymm_(fieldsymm)
  { 
    inditer ii(3,order+fieldorder);
		
		//size of a tensorset consisting of tot.symm.tensors up to order N
		int tssiz = getAsiz(N);
    do{
			int nums[3]={0,0,0};

			const int* indices = *ii;
			for(int i=0;i<order;++i)
				nums[indices[i]]++;

			int index = calcIndex(N,nums[0],nums[1],nums[2]) ;
							
			if(fieldorder==0){
				x[indices] = polynom(1,index);  
			}
			else if( fieldsymm )
			{
				nums[0]=nums[1]=nums[2]=0;
				for(int i=order;i<order+fieldorder;++i)
					nums[indices[i]]++;
				index += tssiz * calcIndex2(nums[0],nums[1],nums[2]);
			}
			else
			{
				int hiind=0;
				for(int i=order+fieldorder-1;i>=order;--i)
					hiind=hiind*3+indices[i];

				index += hiind*tssiz;
			}

			x[indices]=polynom(1,index);

			if(decoder){
				vector<int> & ind  =  (*decoder)[index];
				if(ind.size()==0)
					ind.insert(ind.begin(),indices,indices+order+fieldorder);
			}

    }while(++ii!=0);
  }
  
};


#endif //CREATE_SYMM_3DTENS_H

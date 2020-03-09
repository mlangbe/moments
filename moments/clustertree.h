/*
	tree of clusters
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


#ifndef GNUCC
#pragma once
#endif

#ifndef CLUSTERTREE_H
#define CLUSTERTREE_H

#include<string>
#include<ostream>

class clustertree
{
		

public:
	
	struct node{ 
		//the id of the left/right son if positive
		//or the -1-id of the left/right node if negative
		int left,right;
		
		//the distance between the left/right cluster
		float distance;

		//number of leaves under this node (useful for layout)
		int number;
	};

	//the number of leaves
	int dim;
	//the dim-1 nodes of the tree
	node * tree;	

	
	/**
	create a cluster tree of dim leaves
	using distmat as a dim x dim matrix
	containing the distances between the leaves.
	the entries of distmat with row>col will 
	be overwritten/destroyed	
	\param distmat: a dim x dim matrix
	\param dim : the number of leaves
	*/
	clustertree(float*distmat,int dim);

	~clustertree(void);

};

void drawCorrTreePS(const float*cov,int dim,std::ostream&out,const std::string*dimnames=0);


#endif

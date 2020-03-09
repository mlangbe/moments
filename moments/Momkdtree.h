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


#ifndef GNUCC
#pragma once
#endif

#ifndef MOMKDTREE_H
#define MOMKDTREE_H


#include<vector>
class Momkdtree
{
public:
	//the number type
	typedef long double mynum;

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

	Momkdtree(void);

	//build kdtree from dim-dimensional points stored in mynum
	void init(std::vector<mynum> &moments,int dim);
	
	Momkdtree(std::vector<mynum> &moments,int dim)
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
	get the box which contains coordinate or moment pos
	*/
	bool getBox(mynum*box,const mynum*pos) const;
	
	/**
	*get all leaves in the Ball.
	*\param momids: start index of moment set in moment set array
	*\param center: an array with dim entries containing 
	*the center of the hyperball
	*\param radius: the radius of the hyperball
	*/
	void momIdsInBall(std::vector<int> &momids,
		const mynum*center,mynum radius) const;




	~Momkdtree(void);
};

void testMomkdtree();
#endif

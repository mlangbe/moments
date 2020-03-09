/*
	octree with tensor set in every node used to accelerate the summing up of moment Tensors
	of the points or field inside a  ball mask.

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
#include "tensorSet.h"
#include<vector>

/**
 an octree storing a tensorset in every node
 every leaf has size 1x1x1,
 node at height i has size 2^i x 2^i x 2^i
*/
class TsTree2
{
public:
	typedef	tensorSet3D<4,long double> tset;

	typedef long double mynum;
	
private:

	struct ND
	{
		tset t;

		bool isleaf;
	};
	
	struct Node: public ND
	{	
		//child of type Node or Leaf 
		//the number is 4*z + 2*y + x
		//when x,y,z are the respective bits of the coordinates
		ND*ch[8];
	};
	
	struct Leaf:public ND
	{ 
	  //pointer on children point range
		mynum *begin, *end;
	};

	


  //pointer on the root node or root leaf
	ND * root;
	
	//the nodes
	std::vector<Node> nodes;
	
	//maximum depth of the tree
	int depth;

	//the leaves
	std::vector<Leaf> leaves;

	//the point coordinates
	std::vector<mynum> coords;

	//the min, max, average, variance of coords
	mynum min[3],max[3],avg[3],var[3];

	//edge length of one bucket(at minimum depth)
	mynum rb0;

	//name space for all implementation details.
	struct impl;

public:

	//for stats: number of tsets summed in suminCircle
	static int numbersummed;
	static int ptssummed;
	
	/**sum in Ball given in tree coordinates multiplied by 2.
	* to convert the integer coords to Float, you'd do the following:
	*	mynum s=2.0/rb0;
	*	int center[]={
	*			int( (c[0]-min[0]) *s ),
	*			int( (c[1]-min[1]) *s ),
	*			int( (c[2]-min[2]) *s )
	*		};
	*		int radius=int(r*s);
	*/
	void sumInCircle(
		tset&innerSum,tset&borderSum,
		const int center[3], int radius, int minh
		) const;
	
	//sum in ball given in real coordinates
	void sumInCircle(
														tset&innerSum,tset&borderSum,
														const mynum center[3], mynum radius, int minh
														) const;

	
	//sum in ball given in real coordinates with multiple radii
	//\pre: radii are sorted in ascending order
	void sumInCircle(
														tset*tsets,int n,
														const mynum center[3], 
														const mynum* radii
														) const;
	

	//the fractal dimension of the tree(max: 3)
	long double fractdim() const;

	//get the scale to be applied when converting fromm the tree
	//coordinates to the real coordinates
	long double getScale() const {return rb0;}

	//get the origin of the tree coordinate system
	const long double * getOrig() const {return min;}


	int getDepth()const {return depth;}

	void getBallTsets(
		std::vector<tset> &inner,
		std::vector<tset> &border,
		int minheight, int radius
		) const;

	/**
	*\param pts: points from which the tensorsets are computed
	*\param boxwidth: width of the leaves of the octree expressed 
	*as fraction of the std. deviation of the points
	*/
	TsTree2(const std::vector<long double>&pts , long double boxwidth=0,int ptperleaf=32);
	~TsTree2(void);
};

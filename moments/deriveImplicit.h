/*
	derive implicit representations from polynomial parametrics.

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

#pragma once
#include "polynom.h"

namespace di{

/**
* give the expected order for the implicit
* representation of the parametric manifold with nparam params 
*embedded in ndim dimensions.
*/
int neededOrder(int nparams, int order, int ndim);

/**
* give the expected order for the implicit
* representation of the parametric given by p.
* the dimension of the space is p.size(),
* the number od parameters is derived from the contents of p 
*/
int neededOrder(const vector<polynom> &p );


/**
* get a set of polynomial functions 
* which fullfill out(in)=0
* (i.e. the implicit representation of in)
*/
void getImplicits(vector<polynom> &out,const vector<polynom> &in);




};

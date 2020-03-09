/*
	tensor together with symmetry info

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

#ifndef TENS_WITH_SYMM_H
#define TENS_WITH_SYMM_H

#include "tensor.h"
#include "polynom.h"

/**
* information of the symmetric index groups
* and the tensors a tensor is composed of.
*/
struct tensSymmInfo
{
	//for every group of symmetric indices in this tensor
	//the number of indices belonging to it.
	virtual const std::vector<int>& groups() const =0;

	//for every source tensor of this tensor, the number of
	//index groups belonging to it
	virtual const std::vector<int>& tgroups() const =0;

	//get the orders per subtensor as indicated by tgroups
	void getTensorOrders(std::vector<int> &ord);


	void writeStructure(std::ostream&out) const;

};


struct basicTensSymmInfo:public tensSymmInfo
{
	vector<int> grp,tgrp;
	//for every group of symmetric indices in this tensor
	//the number of indices belonging to it.
	virtual const std::vector<int>& groups() const 
	{return grp;}

	//for every source tensor of this tensor, the number of
	//index groups belonging to it
	virtual const std::vector<int>& tgroups() const
	{return tgrp;}
};



/** a tensor along with information about symmetry groups inside */
struct tensWithSymm: public tensSymmInfo
{
	virtual operator const tens<polynom>& () const =0;
	//for convenience
	const tens<polynom>& t() const{return *this;}
};



struct tensorWithSymm:public tensWithSymm
{
	tensor<polynom> te;
	std::vector<int> grp,tgrp;

	tensorWithSymm(tensor<polynom> &ot):te(ot){}

	virtual operator const tens<polynom>& () const {return te;}
	virtual const std::vector<int>& groups() const {return grp;}
	virtual const std::vector<int>& tgroups() const {return tgrp;}
};



struct tensProdWithSymm:public tensWithSymm
{
	tensProd<polynom> tp;
	std::vector<int> grp,tgrp;
	
	inline tensProdWithSymm(const tensWithSymm&a,const tensWithSymm&b)
		:tp(a,b),grp(a.groups()),tgrp(a.tgroups())
	{
		grp.insert(grp.end(),b.groups().begin(),b.groups().end());
		tgrp.insert(tgrp.end(),b.tgroups().begin(),b.tgroups().end());
	}
	virtual operator const tens<polynom>& () const {return tp;}
	virtual const std::vector<int>& groups() const {return grp;}
	virtual const std::vector<int>& tgroups() const {return tgrp;}
};

#endif

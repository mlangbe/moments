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

#include "deriveImplicit.h"
#include "pointCloud.h"
#include "poly.h"


namespace di{








/**
* give the expected order for the implicit
* representation of the parametric manifold with nparam params 
*embedded in ndim dimensions.
*@return: the homogenous order maximally needed for the resulting polynomials,
*or -1 if more than 1000 
*   parameters are needed,i.e. the matrix is too big.
*/
int neededOrder(int nparams, int order, int ndim)
{

	for(int ord=order;ord<1000;++ord)
	{
		double equations = nuk(order*ord+nparams,nparams);

		double coefficients= nuk(ord+ndim,ndim);
	
		//cout<<' '<<ord<<' '<<equations<<' '<<coefficients<<' ';
		if(coefficients >= equations + ndim-nparams )
			return ord;
		
		if(coefficients>1e9)
			return -1;
	}
	return -1;
}

/**
* give the expected order for the implicit
* representation of the parametric given by p.
* the dimension of the space is p.size(),
* the number od parameters is derived from the contents of p 
*/
int neededOrder(const vector<polyn> &p )
{
	int ndim=p.size();
	int nparams=0;
	int order=0;
	for(int i=0;i<p.size();++i)
	{
		nparams = max(nparams,1+p[i].getMaxParamId());
		order=max(order,p[i].order()); 
	}

	return neededOrder(nparams,order,ndim);
}





void testImplicits()
{
	int p=2,d=3;
	for(int o=0;o<=10;++o)
	{
		int oord = neededOrder(p,o,d);
		cout<<"order"<<o<<"output order="<<oord
			<<" #coeffs:"<<nuk(oord+d,d)
			<<" #equations"<<nuk(oord*o + p,p)
			<<" #input coeffs:"<<nuk(o+p,p) <<endl;
	}
	

	char c;
	cin>>c;
}


};

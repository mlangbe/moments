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


#include "tensWithSymm.h"

//get the orders per subtensor as indicated by tgroups

void tensSymmInfo::getTensorOrders(std::vector<int> &ord)
{
	const vector<int> &tg = tgroups(), &grp =groups();
	
	ord.resize(tg.size(),0);

	vector<int>::const_iterator git=grp.begin();

	for(int i=0;i<tg.size();++i)
		for(int j=0;j<tg[i];++j,++git)
			ord[i] += *git;
}

void tensSymmInfo::writeStructure(std::ostream&out) const
	{
			int ordid=0;
			const vector<int>&tords=tgroups(),&ords=groups();

			for(int i=0;i<tords.size();++i){
				out<<'(';
				for(int j=0;j<tords[i];++j)
				{
					int n=ords[ordid++];
					for(int k=0;k<n;++k)
						out<<'*';
					out<<' ';
				}
				out<<')';
			}
			out<<endl;
	}

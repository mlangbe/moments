/*
     debug output helper

	This file is part of the source code used for the calculation of the moment invariants
	as described in the dissertation of Max Langbein
	https://nbn-resolving.org/urn:nbn:de:hbz:386-kluedo-38558

	Copyright (C) 2020 TU Kaiserslautern,
	              Prof.Dr. Hans Hagen (AG Computergraphik und HCI) (hagen@cs.uni-kl.de)

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

#include "outputflow.h"
using namespace std;

int outputflow::reclevel=0;
ostream * outputflow::dbgout;
int outputflow::tabwidth=2;



void outputflow::outputindent(ostream&out)	
{
	for(int i=0;i<reclevel*tabwidth&&i<20;++i)
		out<<' ';
}


std::ostream& outputflow::out(){
	ostream &o=myout();
	outputindent(o);
	return o;
}



outputflow::outputflow(const std::string& name_,const std::string&attstring)	
{
	ostream&out=myout();

	for(int i=0;i<name_.size();++i)switch(name_[i])
	{
		case':': name+='_';break;
		case'\\': name+="_sl_";break;
		case'.': name+="_dot_";break;
		default: name+=name_[i];
	}

	outputindent(out);

	if(attstring.size()==0)
		out<<'<'<<name<<'>'<<endl;
	else
	{
		out<<'<'<<name<<' '<<attstring<<'>'<<endl;
	}
	++reclevel;
}

outputflow::~outputflow()
{
	--reclevel;

	ostream&out=myout();
	outputindent(out);
	out<<"</"<<name<<'>'<<endl;
}

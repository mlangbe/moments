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

#pragma once

#ifndef OUTPUTFLOW_H
#define OUTPUTFLOW_H

#include<iostream>
#include<string>

class outputflow
{
	static int reclevel;

	std::string name;

	static std::ostream&myout(){return dbgout?*dbgout:std::cout;}

	static void outputindent(std::ostream&out);

public:


	static std::ostream& out();


	static int tabwidth;


	static std::ostream*dbgout;

	//output name and optionally some attributes
	outputflow(const std::string&name_,const std::string& attstring="");
	~outputflow();
};



#ifdef __FUNCTION__
#define OUTPUTFLOW() outputflow outputflow_f (__FUNCTION__)
#else
#define OUTPUTFLOW() outputflow outputflow_l(__FILE__)
#endif



#endif

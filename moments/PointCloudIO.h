/*
	I/O of pointclouds to/from different formats

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

#include<vector>
#include<iostream>
#include "momentSet.h"
typedef faceMomSet::valuetype mynum;



class PointCloudIO
{


protected:

	PointCloudIO(void);
	~PointCloudIO(void);

	virtual void read_(std::vector<mynum>&pts,const char*name)  {
		std::cerr<<"no reader for format "<<format()<<std::endl;
		throw;
	};

	virtual void write_(const std::vector<mynum>&pts,const char*name)  {
		std::cerr<<"no writer for format "<<format()<<std::endl;
		throw;
	};

	virtual const char * ext()  =0;
	virtual const char * format()  {return ext();}

	struct Factory{
		std::vector<PointCloudIO*> instances;		
		PointCloudIO*getIO(const char*name,const char*format=0) ;

	};
	
	static Factory* instance;

	static Factory*getInstance()
	{
		if(!instance)
			instance=new Factory();
		return instance;
	}

	void reg()
	{
		getInstance()->instances.push_back(this);
	}


public:

	bool accept(const char*name,const char*format=0);

	static void read(std::vector<mynum>&pts,const char*name,const char*format=0)
	{
		getInstance()->getIO(name,format)->read_(pts,name);
	}
	static void write(const std::vector<mynum>&pts,const char*name,const char*format=0)
	{
		getInstance()->getIO(name,format)->write_(pts,name);	
	}

	


};




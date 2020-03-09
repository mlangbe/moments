/*
	This file is part of example code using interval arithmetics for implicit surface triangulation

	Copyright (C) 2020 TU Kaiserslautern, Prof.Dr. Hans Hagen (AG Computergraphik und HCI)
	
	Author: Max Langbein	

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

#ifndef IMPLICIT_H
#define IMPLICIT_H

#include<myinterval.h> 




struct Implicit{

	typedef myinterval<float> tval;


	virtual tval operator()(const tval *x) const = 0;
	virtual void grad(tval*ret,const tval*x) const = 0;

	virtual void grad(tval*ret,const float *x) const
	{
		tval r[]={x[0],x[1],x[2]};
		grad(ret,r);
	}

	virtual tval operator()(const float *x) const
	{
		tval r[]={x[0],x[1],x[2]};
		return (*this)(r);
	}


};






#endif
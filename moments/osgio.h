/*
    optimizations operations for polynomials
    (common factor extraction e.g.)

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
#ifndef OSGOUT_H
#define OSGOUT_H

#include<iostream>
#include<vector>

#include "momentSet.h"
typedef faceMomSet::valuetype mynum;


//create set of circles intersecting each other
void osgcircball(std::ostream&out,int num=72);


void osgpointcloud(std::ostream&out,const std::vector<mynum> pts);


void osgplaceandscalemat(std::ostream&out,const mynum*center,
														 mynum s
														 );

void osgplaceandscale(std::ostream&out,const mynum*center,
														 mynum s,
														const char*id="ballid");
void osgconnector(std::ostream&out,const mynum*center,const mynum*center2);


void osgreadpoints(std::vector<mynum>&pts, std::istream& in);


void points_colored_momdist_osg
(const std::vector<mynum>& pts,const std::vector<mynum>& moms,
 const std::vector<mynum>& clustercenters,int momdim,std::ostream&out);

void points_colorcoded_osg
(const std::vector<mynum>& pts,std::vector<int> code,std::ostream&out);

#endif

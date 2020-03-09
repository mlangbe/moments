/*
	representation of a moment invariant
	as an undirected graph of the source tensors, every pair representing a contraction.
	has == and < to allow for sorting and thus normalization.



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

#ifndef VVPR_H
#define VVPR_H


#ifndef GNUCC
#pragma once
#endif

#include<vector>
#include<iostream>

struct pr{
  int a,b;
  inline bool operator < (const pr&x) const
  {
    return a!=x.a ? a < x.a : b<x.b;
  }

  inline bool operator != (const pr&x) const
  {
    return a!=x.a || b!=x.b;
  }

  inline bool operator == (const pr&x) const
  {
    return a==x.a && b==x.b;
  }
};

typedef std::vector< pr > vpr;
typedef std::vector< vpr > vvpr;

typedef std::vector<std::vector<int> > symm;


inline bool operator < (const vpr&x,const vpr&y)
{
  for(int i=0;i<(int)x.size();i++)
    if(x[i]!=y[i])
      return x[i]<y[i];

  return false;
}

inline bool operator == (const vpr&x,const vpr&y)
{
  for(int i=0;i<(int)x.size();i++)
    if(x[i].a!=y[i].a || x[i].b!=y[i].b)
      return false;

  return true;

}

inline std::ostream&operator<<(std::ostream&o,const vpr&x)
{
  for(int i=0;i<(int)x.size();i++)
     o<<"("<<x[i].a<<" \\; "<<x[i].b<<")";
	o<<"\\\\";
  return o;
}

#endif

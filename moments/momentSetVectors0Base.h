/*
	implementations of invariant computing methods from moment tensors.
	(generated code)

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

#ifndef OLD_VC
#define _TEMPLATE_ template<>
#else
#define _TEMPLATE_
#endif



#define MYMOMENTSET momentSet<MOMTYPE,27,0,0>

_TEMPLATE_
const int MYMOMENTSET::base::n;

_TEMPLATE_
const MYMOMENTSET::meta 
MYMOMENTSET::base::metainfo ={
30, 3,3,{
  {2,0,"1A 1A :(0 \; 1)\\"}
  ,
  {6,6,"1A 3A 1A 3A 1A 3A :(1 \; 2)(0 \; 11)(3 \; 9)(4 \; 7)(5 \; 8)(6 \; 10)\\"}
  ,
  {3,1,"1A 1A 2A :(0 \; 3)(1 \; 2)\\"}
  ,
  {1,1,"2A :(0 \; 1)\\"}
  ,
  {2,2,"2A 2A :(0 \; 2)(1 \; 3)\\"}
  ,
  {2,2,"2A 2A :(0 \; 3)(1 \; 2)\\"}
  ,
  {3,3,"2A 2A 2A :(0 \; 5)(1 \; 3)(2 \; 4)\\"}
  ,
  {3,3,"2A 2A 2A :(0 \; 5)(1 \; 2)(3 \; 4)\\"}
  ,
  {4,4,"2A 2A 2A 2A :(0 \; 7)(1 \; 5)(2 \; 6)(3 \; 4)\\"}
  ,
  {3,5,"2A 3A 3A :(2 \; 3)(5 \; 6)(0 \; 7)(1 \; 4)\\"}
  ,
  {3,5,"2A 3A 3A :(2 \; 3)(0 \; 5)(1 \; 4)(6 \; 7)\\"}
  ,
  {3,5,"2A 3A 3A :(5 \; 6)(0 \; 7)(1 \; 2)(3 \; 4)\\"}
  ,
  {2,4,"3A 3A :(0 \; 1)(3 \; 4)(2 \; 5)\\"}
  ,
  {2,4,"3A 3A :(0 \; 1)(2 \; 3)(4 \; 5)\\"}
  ,
  {2,4,"3A 3A :(0 \; 2)(1 \; 3)(4 \; 5)\\"}
  ,
  {2,4,"3A 3A :(0 \; 3)(1 \; 4)(2 \; 5)\\"}
  ,
  {2,4,"3A 3A :(0 \; 3)(1 \; 5)(2 \; 4)\\"}
  ,
  {4,8,"3A 3A 3A 3A :(0 \; 1)(3 \; 4)(6 \; 7)(2 \; 11)(5 \; 9)(8 \; 10)\\"}
  ,
  {4,8,"3A 3A 3A 3A :(0 \; 1)(3 \; 4)(2 \; 11)(5 \; 9)(6 \; 8)(7 \; 10)\\"}
  ,
  {4,8,"3A 3A 3A 3A :(0 \; 1)(3 \; 4)(2 \; 11)(5 \; 8)(6 \; 9)(7 \; 10)\\"}
  ,
  {4,8,"3A 3A 3A 3A :(0 \; 1)(3 \; 4)(2 \; 11)(5 \; 6)(7 \; 9)(8 \; 10)\\"}
  ,
  {4,8,"3A 3A 3A 3A :(0 \; 1)(2 \; 11)(3 \; 6)(4 \; 9)(5 \; 10)(7 \; 8)\\"}
  ,
  {4,8,"3A 3A 3A 3A :(0 \; 1)(2 \; 11)(3 \; 9)(4 \; 10)(5 \; 6)(7 \; 8)\\"}
  ,
  {4,8,"3A 3A 3A 3A :(0 \; 1)(2 \; 11)(3 \; 6)(4 \; 9)(5 \; 8)(7 \; 10)\\"}
  ,
  {4,8,"3A 3A 3A 3A :(0 \; 1)(2 \; 11)(3 \; 6)(4 \; 7)(5 \; 9)(8 \; 10)\\"}
  ,
  {4,8,"3A 3A 3A 3A :(0 \; 1)(2 \; 11)(3 \; 6)(4 \; 9)(5 \; 7)(8 \; 10)\\"}
  ,
  {4,8,"3A 3A 3A 3A :(0 \; 1)(2 \; 11)(3 \; 8)(4 \; 9)(5 \; 6)(7 \; 10)\\"},

 },
	"A moment set for tensor created from a vectorfield.\n"
"The tensor set A is given as 3 Tensorsets of the same maximum order, one after each other (one for every vector component)"
};

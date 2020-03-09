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

#define MYMOMENTSET momentSet<MOMTYPE,28,0>

#ifndef OLD_VC
template<>
#endif
const int MYMOMENTSET::base::n;

#ifndef OLD_VC
template<>
#endif
const MYMOMENTSET::meta 
MYMOMENTSET::base::metainfo ={
 35, 3,4, {
 {1,2,"2A :(0 \\; 1)\\ " } 
 ,
 {2,4,"2A 2A :(0 \\; 2)(1 \\; 3)\\ " } 
 ,
 {3,6,"2A 2A 2A :(0 \\; 2)(1 \\; 4)(3 \\; 5)\\ " } 
 ,
 {4,10,"2A 2A 2A 4A :(0 \\; 6)(1 \\; 7)(2 \\; 4)(3 \\; 8)(5 \\; 9)\\ " } 
 ,
 {4,10,"2A 2A 3A 3A :(4 \\; 5)(0 \\; 7)(1 \\; 8)(2 \\; 6)(3 \\; 9)\\ " } 
 ,
 {3,8,"2A 2A 4A :(4 \\; 5)(0 \\; 2)(1 \\; 6)(3 \\; 7)\\ " } 
 ,
 {3,8,"2A 2A 4A :(0 \\; 4)(1 \\; 5)(2 \\; 6)(3 \\; 7)\\ " } 
 ,
 {3,8,"2A 3A 3A :(2 \\; 3)(5 \\; 6)(0 \\; 4)(1 \\; 7)\\ " } 
 ,
 {3,8,"2A 3A 3A :(2 \\; 3)(0 \\; 5)(1 \\; 6)(4 \\; 7)\\ " } 
 ,
 {3,8,"2A 3A 3A :(0 \\; 2)(1 \\; 5)(3 \\; 6)(4 \\; 7)\\ " } 
 ,
 {2,6,"2A 4A :(2 \\; 3)(0 \\; 4)(1 \\; 5)\\ " } 
 ,
 {3,10,"2A 4A 4A :(2 \\; 3)(6 \\; 7)(0 \\; 4)(1 \\; 8)(5 \\; 9)\\ " } 
 ,
 {3,10,"2A 4A 4A :(2 \\; 3)(0 \\; 6)(1 \\; 7)(4 \\; 8)(5 \\; 9)\\ " } 
 ,
 {3,10,"2A 4A 4A :(0 \\; 2)(1 \\; 6)(3 \\; 7)(4 \\; 8)(5 \\; 9)\\ " } 
 ,
 {2,6,"3A 3A :(0 \\; 1)(3 \\; 4)(2 \\; 5)\\ " } 
 ,
 {2,6,"3A 3A :(0 \\; 3)(1 \\; 4)(2 \\; 5)\\ " } 
 ,
 {4,12,"3A 3A 3A 3A :(3 \\; 4)(6 \\; 7)(9 \\; 10)(0 \\; 5)(1 \\; 8)(2 \\; 11)\\ " } 
 ,
 {4,12,"3A 3A 3A 3A :(9 \\; 10)(0 \\; 3)(1 \\; 6)(2 \\; 11)(4 \\; 7)(5 \\; 8)\\ " } 
 ,
 {4,12,"3A 3A 3A 3A :(0 \\; 6)(1 \\; 9)(2 \\; 10)(3 \\; 7)(4 \\; 8)(5 \\; 11)\\ " } 
 ,
 {4,12,"3A 3A 3A 3A :(0 \\; 3)(1 \\; 6)(2 \\; 9)(4 \\; 7)(5 \\; 10)(8 \\; 11)\\ " } 
 ,
 {3,10,"3A 3A 4A :(0 \\; 1)(2 \\; 6)(3 \\; 7)(4 \\; 8)(5 \\; 9)\\ " } 
 ,
 {1,4,"4A :(0 \\; 1)(2 \\; 3)\\ " } 
 ,
 {2,8,"4A 4A :(0 \\; 1)(4 \\; 5)(2 \\; 6)(3 \\; 7)\\ " } 
 ,
 {2,8,"4A 4A :(0 \\; 4)(1 \\; 5)(2 \\; 6)(3 \\; 7)\\ " } 
 ,
 {3,12,"4A 4A 4A :(0 \\; 1)(4 \\; 5)(8 \\; 9)(2 \\; 6)(3 \\; 10)(7 \\; 11)\\ " } 
 ,
 {3,12,"4A 4A 4A :(4 \\; 5)(8 \\; 9)(0 \\; 6)(1 \\; 7)(2 \\; 10)(3 \\; 11)\\ " } 
 ,
 {3,12,"4A 4A 4A :(8 \\; 9)(0 \\; 4)(1 \\; 5)(2 \\; 6)(3 \\; 10)(7 \\; 11)\\ " } 
 ,
 {3,12,"4A 4A 4A :(0 \\; 4)(1 \\; 5)(2 \\; 8)(3 \\; 9)(6 \\; 10)(7 \\; 11)\\ " } 
	}
};

#undef MYMOMENTSET

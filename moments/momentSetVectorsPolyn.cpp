/*
	implementations of invariant computing methods from moment tensors.

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

#include "polynom.h"
#include "momentSet.h"

//we don't want to wait forever!
#pragma optimize("", off)

#define MOMTYPE polyn
#include "momentSetVectors0Base.h"
#include "momentSetVectors00.h"
#include "momentSetVectors01.h"
#include "momentSetVectors02.h"
//for computation time:
//#include "momentSetVectors03.h"
//#include "momentSetVectors04.h"
#undef MOMTYPE



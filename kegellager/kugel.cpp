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

#include<cmath>
#include<string.h>
#include<stdlib.h>
#include<iostream>

#include "Objects.h"
#include "extractsurface.h"

int main(int argc,const char*argv[])
{
 myinterval<float> inters[]={myinterval<float>(-2,2),myinterval<float>(-2,2),myinterval<float>(-2,2)};

 Kugel k(1);

 extractsurface().toStl(k,"kugel1_200.stl",1.0/200,inters); 

 return 0;
}

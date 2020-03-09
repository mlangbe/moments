/*
	implementations of operations on polynomial-valued tensors.

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

#ifndef MOMENTS_H
#define MOMENTS_H
#include <time.h>
#include <iostream>

#include "momentSet.h"

template<int numc,class momtype>
void testMomFunSpeed2()
{
	using namespace std;
	int num= 100*1000;
	int numm=momtype::n;

	cout<<"testing " << momtype::comment <<endl;
	cout<<"number components:"<<numc<<endl;
	mytype *A=new mytype[numc*num];
	mytype *M=new mytype[ momtype::n * num];

	for(int i=0;i<numc*num;++i)
	{
		A[i]=rand();
	}
	
	clock_t st=clock();
	for(int j=1;j<=100 && (clock()-st) < 2*CLOCKS_PER_SEC ;j++)
	{
		mytype* Ai=A,*Aend= A+numc*num,*Mi=M;
		for(;Ai!=Aend;Ai+=numc,Mi+=numm){
			momtype::compute(A,M);
		}

		double nspt=1e9*(clock()-st)/CLOCKS_PER_SEC/j/num;
		cout<<"time per tset :"<< nspt	
		<<" ns time per component:"<<nspt/numc<<" ns"<<endl;
	}

}

/*
*create a c++ file that is a specialization of moment set
*from the pairings 
*\param optimize: do quasi-horner optimization
*\param docommonsubopt: do common subexpression optimization
*/
void createVecMomSpecialization(std::ostream& of,
																int id=0,
																bool optimize=true,
																bool docommonsubopt=true,
																const char*comment=0);

/**
create an independent set of moment invariants of
scalar field tensor up to order 4 .
save the set in 
momScal_raw.txt
*/
void getMomentsScalarField(std::ostream&out,std::ostream& polysout);



#endif

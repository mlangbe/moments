/*
	implementations for operations to create sets of independent polynomials, i.e. there is no function with p_n(x)= f( p_1(x)... p_n-1(x))

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

#include "infbruch.h"
#include "polynom.h"


#include<vector>
#include<set>


//fill p with prime numbers;numthere:number thats  are already there
void getPrimes(vector<int>&p,int numthere=0);


class indepPolySet
{
public:
  //typedef bruch<long long int > mynum;
  typedef infbruch mynum;
  //typedef myinterval<long double> mynum;

	typedef std::set<polynom> tpolyset;
private:
	//the param ids in the order they occur in the matrix
  std::vector<int> paramIds;

	//the example values for the parameters(are sqrt's of prime numbers)
  std::vector<mynum> paramValues;
	
  tpolyset polys;

  // a vector of prime numbers
  std::vector<int> primes;

  //is parameter with a certain id present ?
  std::vector<bool> paramThere;


	//use the predefined values or use prime numbers ?
	bool usePresetParamValues;

  //matrix to determine if a polynom is independent of the others
  std::vector<std::vector<mynum> >mat;
 
public:

	indepPolySet()
	{usePresetParamValues=false;}

	indepPolySet(const std::vector<mynum>& paramValues_)
	{usePresetParamValues=true;paramValues=paramValues_;}

	void setParamValues(const std::vector<mynum>& paramValues_);

	void printMat(std::ostream&out);
  const tpolyset& p() const
  { return polys; }

  bool addpolynom(const polynom&x);
  
	int join(const indepPolySet&xx);
	
	static void save(ostream&o,const tpolyset&a);
	
	static void load(istream&in,tpolyset&a);
};







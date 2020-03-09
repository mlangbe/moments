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

#ifndef OPTIMIZED_POLY
#define OPTIMIZED_POLY

#include "polynom.h"

#ifndef GNUCC
#pragma once
#endif


class optimizedPoly
{
	//1st int in pair: param id if there are children in second
	//coefficient if there are no children in second

public:
	typedef std::pair<int, optimizedPoly > summand;

	typedef std::list<summand> summands;
	summands children;

	//default constructor so we can use the vector
	optimizedPoly()
	{}

	optimizedPoly(polyn p)
	{
		init(p.s);	
	}

	void init(std::list<polyn::summand>& s);

	//returns 0 if equal, -1 if less, 1 if greater
	static int vgl (
		pair<int,optimizedPoly> &a,
		pair<int,optimizedPoly> &b
		)
	{
		if( a.first < b.first ) return -1;
		if( b.first < a.first ) return 1;
		return a.second.vgl(b.second);
		return 0;
	}

	int vgl(optimizedPoly&b);


	bool operator < (optimizedPoly &b)
	{	return vgl(b)<0;	}

	bool operator == (optimizedPoly &b)
	{	return vgl(b)==0;	}

	template<class T>
	T eval(const std::vector<T> &params) const;

	//produce sth C++ can eat
	friend std::ostream& operator << (std::ostream&o, const optimizedPoly & p);

	friend class optimizedPolyCommonSub;
};

inline bool operator ==  (
		std::pair<int,optimizedPoly> &a,
		std::pair<int,optimizedPoly> &b
		)
{	return optimizedPoly::vgl(a,b)==0;}

inline bool operator < (
		std::pair<int,optimizedPoly> &a,
		std::pair<int,optimizedPoly> &b
		)
{	return optimizedPoly::vgl(a,b)<0; }







//optimize for common subexpressions,
//base for ... commonSub, ....Multi
class optimizedPolyCommonSubBase
{


	
protected:

	struct cmpProd
	{
		inline bool operator()(
			optimizedPoly::summands::iterator a,
			optimizedPoly::summands::iterator b
			) 
		{
			return *a < *b;	
		}
	};

	struct cmpSum
	{
		inline bool operator()(
			optimizedPoly* a,
			optimizedPoly* b
			) 
		{
			return *a < *b;	
		}
	};




	//the index will be the order of the polynomial in DNF
	typedef std::vector< optimizedPoly::summands::iterator > prods;
	std::vector< prods > prodOrders;
	typedef std::vector< optimizedPoly*> sums;
	std::vector< sums > sumOrders;
	
	int maxid;

	//the new subexpressions that get assigned a polynome id, the original 
	//coming first.
	//ids are starting with maxid, i.e. to got a poly with id k,
	//you have to access element k-maxid
	std::vector<optimizedPoly> np;


	//returns the maximum order of p,
	//fills prods,sums.
	int recinit(optimizedPoly & p);
		
	//used in the init functions of the subclasses
	//creates the np array and replaces
	//subexpressions by ids in np
	//maxorder=maximum order in parameters of all polynomes
	//pre: prods,sums have been filled
	void initBase(int maxorder);

public:
	int getMaxId(){return maxid;}

};

//optimize for common subexpressions
class optimizedPolyCommonSub : optimizedPolyCommonSubBase
{



	//the original with subexpressions replaced with ids in np
	optimizedPoly op;
		
	void init(const optimizedPoly& p);

public:
	optimizedPolyCommonSub(const optimizedPoly&p)
	{init(p);}

	friend 
	std::ostream&operator<<(std::ostream&o,optimizedPolyCommonSub &s);

	template<class T>
	T eval(const std::vector<T> &params) const;


};


//optimize for common subexpressions
class optimizedPolyCommonSubMulti : public optimizedPolyCommonSubBase
{


	/**the original polynomes with subexpressions replaced
	by ids in np
	*/
	std::vector<optimizedPoly> op;
		
	void init(const vector<optimizedPoly>& p);

public:
	optimizedPolyCommonSubMulti(const vector<optimizedPoly>&p)
	{init(p);}

	ostream& print(std::ostream &o,const std::vector<metainfo*>* mi=0) const;

	template<class T>
	void eval(
		const std::vector<T> & params,
		std::vector<T> & ret) const;

};

ostream&operator<<(ostream&o,optimizedPolyCommonSubMulti&s);


template<class T>
T optimizedPoly::eval(const std::vector<T> & params) const
{
	summands::const_iterator i;
	T ret=0.0,buf;
	for(i=children.begin();i!=children.end();++i)
	{
		if(i->second.children.empty())
			ret+=T(i->first);
		else
			ret+=params[i->first] * i->second.eval(params);
	}
	return ret;
}

template<class T>
T optimizedPolyCommonSub::eval(const std::vector<T> & params) const
{
	assert(params.size()>maxid);
	typename std::vector<T> params2(maxid+1+np.size());	
	std::copy(params.begin(),params.begin()+maxid+1,params2.begin());	
	
	for(int i=0;i<np.size();++i)
	{
		params2[maxid+1+i]=np[i].eval(params2);	
	}

	return op.eval(params2);
}


//not yet tested
template<class T>
void optimizedPolyCommonSubMulti::eval(const std::vector<T> & params,
																	 std::vector<T> & ret) const
{
	assert(params.size()>maxid);
	typename std::vector<T> params2(maxid+1+np.size());	
	std::copy(params.begin(),params.begin()+maxid+1,params2.begin());	
	
	for(int i=0;i<np.size();++i)
	{
		params2[maxid+1+i]=np[i].eval(params2);	
	}

	ret.resize(op.size());
	for(int i=0;i<op.size();++i)
		ret[i]=op[i].eval(params2);
}


#endif

/*
    contraction operations on tensors with polynomial entries and certain symmetries

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

#ifndef OPTIMIZED_FOLD_H
#define OPTIMIZED_FOLD_H

#include<set>
#include<sstream>
#include<algorithm>


#include "tensor.h"
#include "polynom.h"
#include "refp.h"


//storing of tensors with repetitions in the components
class repeatTensor: public tens<polyn>
{
public:
	static int nextid;
	
	int id;
	
	//tensor with ids into comps
	tensor<int> t;

	//neded to store number of comps if comps vector is empty.
	int numcomps;

	//the components the ids in t refer to
	vector<polyn> comps;

	//commment to be written into the c-code
	string comment;
	
	/** create the tensor indices for very component*/
	void createTindexPerComp(vector<int>&idxpert) const;
	
	/*decompose tensor index ii into sinle indices and write it to ret*/
	string& writeTindex(string&ret,int ii) const;

	virtual ~ repeatTensor();


protected:



	inline static string id2name(int i)
	{
	  string s="";
		do{
			s+= char( 'a' + (i%26) );
			i/=26;
		}while(i);

		return s;
	}


	inline repeatTensor(int ord,int dim,const string& nam="")
		:tens<polyn>(ord,dim),t(ord,dim),comment(nam)
	{
		numcomps=-1;
		id=nextid;
		++nextid;
	}
		
private:
		
	//create a map with unique polynomes
	//which are components of the tensor
	struct tvisitor{
	  int nextid;
		
		map<polyn,int> s;
		vector<int>::iterator tit;
		const tens<polyn> & t;
		
		inline void operator() (int * idx) 
		{
			polyn r=t[idx];
			
			pair< map<polyn,int>::iterator , bool > 
			found = s.insert(pair<polyn,int>(r,nextid) );
			if(found.second)
				++nextid;
			*tit = found.first->second;  				
			++tit;
		}
		
		tvisitor(const tens<polyn> &ot , tensor<int> &nt )
			:t(ot),tit(nt.t.begin())
		{
			nextid=0;	
		}
	};

	//help function for constructors
	void init(const tens<polyn>&tt);
		
public:
	repeatTensor(const tens<polyn>&tt,const string& nam="");
	repeatTensor(const tens<polynom>&tt,const string& nam="");

	inline string id2name() const	{	return id2name(id);	}


  inline const polyn & get(const int*inds) const
  {
		return comps[t.get(inds)];  
  }

	virtual polyn operator[](const int* inds) const;


	//print independent components as polynomes in 
	//the source tensor called A
	virtual void printForC(ostream&o,const char*varname=0) const;


};

/*
	a tensor where operator [] returns the id of the independent   component of t	, with addme added
*/
struct repeatTensorAdd: public tens<polyn>
{

	const repeatTensor & t;
	int addme;
	

	inline repeatTensorAdd(const repeatTensor& tt,int add)
		:tens<polyn>(tt),t(tt),addme(add)
	{	
	}

	virtual polyn operator[](const int* inds) const;

};



/**
*creates a repeat tensor from a tensor 
*whose components are unique if evaluated with the Belegung
*bel
*/
class repeatTensorUneval:public repeatTensor
{

public:

	//the unevaluated components of the tensor	
	vector<polyn> compsUneval;

protected:

	struct creationVisitor
	{
		int id;

		//the belegung of the parameters in ps
		const vector<polyn>& comps;

		//the set with unevaluated polynomes
		map< polyn,int > ps;

		//the map with evaluated polynomes
		//along with information about the unevalueted polyn
		//and the id
		typedef map<polyn, map<polyn,int>::iterator > tevalps;	
		tevalps evalps;

		const tens<polyn> & t;	

		//ids to be filled
		vector<int>::iterator indit;

		inline void operator()(int*inds)
		{
			polyn ret=t[inds];

			map<polyn,int>::iterator res = ps.lower_bound(ret);

			//if it is not in map: try to evaluate it and look if it's still new:
			if( res==ps.end() || res->first!=ret)
			{
				polyn evret = ret.eval(comps);			
				tevalps::iterator resEval = evalps.lower_bound(evret);

				//still new:
				if(resEval==evalps.end() || resEval->first!=evret )
				{
					map<polyn,int>::iterator
						insret = ps.insert(res , pair<polyn,int>(ret,id));				
					evalps.insert( resEval, tevalps::value_type(evret,insret));
					*indit = id;
					//use new id:
					++id;
				}
				else
				{
					*indit = resEval->second->second;
					//cerr<<"second check was useful!"<<endl;
				}
			}
			else //already there: use that id.
			{
				*indit = res->second;
			}

			++indit;
		}

		inline creationVisitor(const tens<polyn> & tt,vector<int>::iterator ii,const vector<polyn>&comps_)
			:t(tt),indit(ii),comps(comps_),id(0)
		{				
		}

	};


	inline repeatTensorUneval(int ord,int dim,const string&nam)
		:repeatTensor(ord,dim,nam)
	{
  }
  
  
  void init(const tens<polyn>&a , const vector<polyn>&bel);

	virtual ~repeatTensorUneval();


};




struct polyn2str
{
	ostringstream s;
	
	inline polyn2str(double x)
	{
		if(x!=0&&x!=1)
		{
			s<<x;
		}
	}

	inline polyn2str(const polyn2str&x)
	{
		s<<x.s.str();
	}
	
	inline polyn2str() {}
	
	inline polyn2str(const string&y)
	{
		s<<y;
	}
	
	inline void operator*= (const polyn2str & x)
	{
		if(s.str().size()!=0)
			s<<"*";
			
		s<<x.s.str();
	}
	
	inline void operator+= (const polyn2str &x)
	{
	if(s.str().size()!=0)
	  s<<"+";
	s<<x.s.str();	
	}


};


class RTpartialFold :public repeatTensorUneval
{
public:
	//a tensor whose components are 
	//polynomes in the independent components of the repeatTensor aa.t
	int aid;
	int acomps;

	vpr pairs;
	vi perm;


	RTpartialFold(const repeatTensor*a,
		const vi & perm,
		const vpr&inds,const string&nam="");

	
/**
* print this polynome as polynomial in the components of 
* the two tensors it was produced from.
* varname is inserted if result has only one variable
*/
  virtual void printForC(ostream&o,const char*varname=0) const;

	virtual ~RTpartialFold();

};

class RTpartialFoldProd: public repeatTensorUneval
{

public:

	int aid,bid;
	int acomps,bcomps;

	double estimateEvalPolyMem(const repeatTensor *a,const repeatTensor *b);

  vpr pairs;
	vi perma, permb;
  
	RTpartialFoldProd
		(	const repeatTensor *a,
			const vi&perma,
			const repeatTensor *b,			
			const vi&permb,
			const vpr &indpairs);

/**
* print this polynome as polynomial in the components of 
* the two tensors it was produced from.
* varname is inserted if result has only one variable
*/
	virtual void printForC(ostream&o,const char*varname=0) const;

	virtual ~RTpartialFoldProd();

};



#endif//OPTIMIZED_FOLD_H

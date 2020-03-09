/*
	operations on polynomials and polynomial-valued tensors.

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



#ifndef GNUCC
#pragma once
#endif

#ifndef POLYNOM
#define POLYNOM

#include<vector>
#include<map>
#include<set>
#include<list>
#include<iostream>
#include<algorithm>
#include<cassert>
#include<string.h>

#include "vvpr.h"

using namespace std;

struct polynom;
ostream&operator<<(ostream&o,const polynom&p);


typedef vector<int> vi;
typedef vector<vi> vvi;

inline bool operator < (const vi &a, const vi &b)
{
  return lexicographical_compare(a.begin(),a.end(),b.begin(),b.end());
}

template<class T>
inline bool veq (const vector<T>&a, const vector<T> &b)
{
  if(a.size()!=b.size())return false;
  for(int i=0;i<(int)a.size();++i)
    if(a[i]!=b[i])return false;
  return true;
}

inline bool operator == (const vi &a, const vi &b)
{ return veq(a,b);}


struct metainfo
{
	static vector<metainfo*> * types;

	virtual ~metainfo();
	virtual metainfo* clone() const =0;
	virtual const char* name() const =0;
	virtual void print(ostream&o) const=0;
	virtual void save_(ostream&o) const=0;
	virtual void load_(istream&in)=0;

	void save(ostream&o) const
	{
		o<<name()<<'\n';
		save_(o);
	}

	static metainfo * load(istream&in)
	{
		if(types==0)
			return 0;

		char x[80];
		in>>x;
		x[79]=0;
		for(unsigned i=0;i<types->size();++i)
			if(strcmp((*types)[i]->name(),x)==0)
			{
				metainfo*nw = (*types)[i]->clone();
				nw->load_(in);
				return nw;
			}
		return 0;
	}
};

template<class X>
struct tmireg
{
	X i;
	
	tmireg()
	{		
		if(!metainfo::types)
			metainfo::types=new vector<metainfo*>();
		
		metainfo::types->push_back(&i);
	}
};


inline ostream&operator<<(ostream&o,const metainfo&m)
{ m.print(o);return o;}



struct tensSymmInfo;

//meta information for one polynome
//created from a tensor
struct tensmetainfo :public metainfo
{
	/*a vector with the ids of the tensors
	the polynome was constructed from.	
	(e.g. a tensor of order 3 is represented by a 3)
	if you have a source tensor constructed as 3A 3A 4A
	vi will be 3 3 4 .
	*/
	vector<int> sourceTensors;

	//the variables from tensSymmInfo
	vector< int > grp,tgrp;

	
	/*the pairs for the paired summing*/
	vpr pairs;
	
	virtual ~tensmetainfo();
	
	metainfo*clone() const;

	void save_(ostream&out) const;
	void load_(istream&in);
	const char*name() const;
	void print(ostream&o) const;

	/*print it out in the "Einstein sum convention" format
	*/
	void printEsumForTex(ostream&o) const;

	/*print it out in the "pairwise summation of tensor product"
	*format
	*/
	void printPairsProdForTex(ostream&o) const;

	/*
	*print it out as a postcript pictore of the 
	*corresponding graph
	*\param grp: the symmetric index groups of every tensor type
	*/
	void printGraphPS(std::ostream&o) const;


	void init	( const tensSymmInfo&pi ,const vpr & pairs_);
};




//represents a polynom composed of a set of products,
//where the coefficients or variables are represented by their ids.
struct polynom
{

	metainfo * mi;

  vvi p;


  polynom(){mi=0;}


	void setmetainfo(const metainfo&m)
	{
		if(mi)delete mi;
		mi=m.clone();		
	}

	polynom(const polynom&o){ 
		p=o.p; mi=0;
		if(o.mi)
			mi =o.mi->clone();
	}
  
	~polynom(){delete mi;}

  polynom(const char* s)
  {
		mi=0;
    vi c;
    for(;*s!=0;++s)
    {
      if('a'<=*s && *s <= 'z')
        c.push_back(*s-'a');
      else
        sort(c.begin(),c.end()),
        p.push_back(c),
        c.clear();
    }    

    if(!c.empty())
      sort(c.begin(),c.end()),
      p.push_back(c);

    sort(p.begin(),p.end());
  }

	/*
  //a polynom consisting of one single variable
  explicit polynom(int id)
  {
		mi=0;
    p.push_back(vi(1,id));
  }
*/
	//a polynom consisting of a variable id number times
  explicit polynom(int number,int id)
		:p(number,vi(1,id)),mi(0)
  {}


  bool operator < (const polynom &x) const
  { return lexicographical_compare(p.begin(),p.end(),x.p.begin(),x.p.end()); }

  bool operator == (const polynom &x) const
  { return veq(p,x.p); }

  bool operator != (const polynom &x) const
  { return !veq(p,x.p); }


  void operator += (const polynom& x)
  {
    int osiz=p.size();
    p.resize(osiz + x.p.size());
    copy(x.p.begin(),x.p.end(),p.begin()+osiz);
    inplace_merge(p.begin(),p.begin()+osiz,p.end());
  }

  polynom operator * (const polynom &x) const
  {
    polynom ret;
    ret.p.resize(x.p.size()*p.size()) ; 

    vvi::const_iterator i,j;
    
    vvi::iterator r=ret.p.begin();

    for(i=p.begin();i!=p.end();++i)
      for(j=x.p.begin();j!=x.p.end();++j,++r)
      {
        r->resize(i->size()+j->size());
        merge(i->begin(),i->end(),j->begin(),j->end(),r->begin());
      }
    sort(ret.p.begin(),ret.p.end());

    //cout<<"="<<ret<<endl;
    return ret;
  }

  void operator *= (const polynom &x)
  {
    //cout<<"mul"<<*this<<' '<<x<<endl;
    (*this)=(*this)*x;
  }

  polynom operator + (const polynom &x) const
  {
    polynom ret=*this;
    ret+=x;
    return ret;
  }

  int getMaxParamID() const
  {
    int idmax=0;
    vvi::const_iterator i;
    vi::const_iterator j;
    for(i=p.begin();i!=p.end();++i)
      for(j=i->begin();j!=i->end();++j)
        if(*j>idmax)idmax=*j;
    return idmax;
  }


	template<class T>
  T eval(const vector<T> & values) const
  {
    T ret=(T)0.0;
    vvi::const_iterator i;
    vi::const_iterator j;
    for(i=p.begin();i!=p.end();++i){
      T x=(T)1.0;
      for(j=i->begin();j!=i->end();++j)
        x*=values[*j];
      ret+=x;
    }
    return ret;     
  }

	template<class T>
  //calculate value of derivation after parameter id
  T evalderiv(const vector<T> & values,int id) const
  {
    T ret=(T)0.0;
    vvi::const_iterator i;
    vi::const_iterator j;
		T x;
    for(i=p.begin();i!=p.end();++i){			

			//count the number this summand occurs
			
			int num=1;
			
			while((i+1)!=p.end() &&*i==*(i+1))
				++i,++num;

			//cout<<"num"<<num<<endl;
			
			//order of id in current summand
			int k=0;

			j=i->begin(); 
			while(j!=i->end() && *j!=id )++j;
			while(j!=i->end() && *j==id)     
            ++k,++j;

			//if order is zero, go to next summand
			if(!k) continue;

			x=k*num;
			for(j=i->begin();j!=i->end();){
        while(j!=i->end() && *j!=id)
          x*= values[*j],++j;

        if(j!=i->end()){          
            ++k,++j;
            while( j!=i->end() && *j==id) 
              x*= values[*j],++k,++j;
        }
      }

			ret+=x;

    }
    return ret;     
  }

	/*should we output the raw indices for the parameters, or
	* were the indices in the polynome created by encodeIndices,
	* so we should use thes in these output << ? */
	enum ioFormats{ RAW,PRINTFORC, PRINTINDEX };
	static int ioFormat;
	friend ostream&operator<<(ostream&o,const polynom&p);
	friend istream&operator>>(istream&in,polynom&p);
	
	/**print the polynome so it could be used in a C program
	* ( use A[id] for the parameters, don't use ^ for powers )
	*/
	void printForC(ostream&o) const;


	void printRaw(ostream&o) const;
	void readRaw(istream&in);


};


//create an index encoding order,dimension and the indices
int encodeIndices(const int*inds,int ord, int dim, int id);

//decode an index created by calcindex
void printindex(ostream&o,int i);


//derive all comparison ops for T from <, ==
#ifdef GNUCC
#define MYDERIV(T,op,negate,code) \
inline bool operator op (const T&b) const\
{ return negate (code);}
#else
#define MYDERIV(T,op,negate,code) \
inline bool operator op (const T&b) const\
{ return negate ## (code);}
#endif

#define DERIVECOMPARISONS(T) \
MYDERIV(T,>,,b<*this) \
MYDERIV(T,>=,!,*this<b) \
MYDERIV(T,<=,!,b<*this) \
MYDERIV(T,!=,!,b==*this)





/*a polynome which has powers for the factors and coefficients
*for the summands
*/
struct polyn
{
  struct factor {
    //id, power
    int i,p;
		
    bool operator < (const factor & x) const
    {
      return i < x.i || i==x.i && p<x.p;
    }

    bool operator == (const factor & x) const
    { return i == x.i && p==x.p; }

		DERIVECOMPARISONS(factor)
	};
	

	struct lti{
		bool operator()(const factor&x,const factor&y)
		{
			return x.i<y.i;
		}
	};


  struct summand {
    //coefficient
    double c;
    //ids with power
    vector<factor> f;

	/** homogenous order*/
	int order() const 
	{
		int ret=0;
		for(vector<factor>::const_iterator 
			cit=f.begin();cit!=f.cend();++cit)
				ret+=cit->p;
		return ret;
	}

	int getMaxParamId() const
	{
		int ret=0;
		for(vector<factor>::const_iterator 
		  cit=f.begin();cit!=f.cend();++cit)
			ret=max(ret,cit->i);
		return ret;
	}

	/**
	return 1 if this >b,
	0 if this=b
	-1 if this < b
	*/
	int cmpf(const polyn::summand &b) const
	{
		const polyn::summand &a=*this;
		//for every param id, the order
		typedef vector<polyn::factor>::const_iterator vpit;
		//remove common divisors from f1,f2

		vpit ita=a.f.begin(),itb=b.f.begin(), nita=ita, nitb=itb;

		while(ita!=a.f.end() && itb!=b.f.end() && *ita==*itb)
		{
			++ita;++itb;
		}

		//b cannot be greater a in this case: 
		if(itb==b.f.end()&&ita==a.f.end())
			return 0;

		//if there are still powers of b, but not of a, a<b.
		if(ita==a.f.end())
			return -1;
		
		if(itb==b.f.end())
			return 1;

		//next var. id with power!=0 in a is greater than that of b
		if(ita->i < itb->i)
			return 1;

		if(ita->i > itb->i)
			return -1;

		return ita->p-itb->p;
	
	}






	  /** only compare the variable orders, not coeffs*/
	//lexicograph order by variables:
	//compare power of var with lowest id first, then that with the second lowest, etc. 
	//cannot use std. lex.comp because var id compare is reversed to sorting of 
	//vars in factors.
	bool operator < (const polyn::summand &b) const
	{
		int i=cmpf(b);
		return i <0 || i==0 && c<b.c; 
	}



    bool operator == (const summand & x) const
    { return c==x.c && veq(f,x.f); }


		DERIVECOMPARISONS(summand)


    summand& operator *= (const summand&x)
    {
			//multiply coefficients
			c*=x.c;

			if(f.empty())
			{
				f=x.f;
				return*this;
			}
			if(x.f.empty())
			{
				return *this;
			}

			size_t osiz=f.size();
      f.resize(osiz + x.f.size());
      copy(x.f.begin(),x.f.end(),f.begin()+osiz);
			//merge the 2 ranges which are sorted by their param ids;
      inplace_merge(f.begin(),f.begin()+osiz,f.end(),lti());
      vector<factor>::iterator i,j;
      i=j=f.begin();
      for(++i;i<f.end();++i)
        if( i->i==j->i )
          j->p += i->p;
        else
          *(++j)=*i;
      //shorten to actual length
      f.resize(j-f.begin()+1);


      return *this;
    }

		friend ostream& operator <<(ostream&o,const summand&s);

		/*replace parameters by the strings given in names
		*and print it.
		*\param withcoeff: also print the coefficient.
		*\param expandPowers: 
		* should we expand powers (i.e x^3= x*x*x ?)
		*/
		void printWithParamNames(ostream&o,
			const std::vector<std::string>&names,
			bool expandPowers=true,
			bool withcoeff=true) const;

  };

  struct ltf
  {
	bool operator()(const summand&a,const summand&b)
	{
		return a.cmpf(b)<0;
	}
  };

 


	

	/*****************************************************/

  list<summand> s;

  /** homogenous order*/
  int order() const {
	  int ret=0;
	  for(list<summand>::const_iterator 
		  cit=s.begin();cit!=s.cend();++cit)
			ret=max(ret,cit->order());

	  return ret;
  }

  int getMaxParamId() const
  {
	  int ret=0;
	  for(list<summand>::const_iterator 
		  cit=s.begin();cit!=s.cend();++cit)
			ret=max(ret,cit->getMaxParamId());

	  return ret;
  }



	polyn(){};

	polyn(const polynom &p);


	/*
	explicit polyn(int id)
	{ 
		factor f;  f.p=1;f.i=id;
		summand ss;
		ss.c=1;
		ss.f.push_back(f);
		s.push_back(ss);  
	}
	*/
	polyn(double coeff)		
		:s(1)
	{
		if(coeff==0)
			s.clear();
		else
			s.front().c=coeff;
	}

	polyn(double coeff,int id)
	{ 
		if(coeff!=0){
			factor f;  f.p=1;f.i=id;
			summand ss;
			ss.c=coeff;
			ss.f.push_back(f);
			s.push_back(ss);  
		}
	}

  bool operator < (const polyn&b) const
  {
  	return 
  	  lexicographical_compare
  	  (
  	   s.begin() , s.end(),
  	   b.s.begin(),b.s.end()
  	  ); 
  }


	bool operator == (const polyn&b) const
	{
		return s==b.s;
		/*
		list<summand>::const_iterator ait,bit;

		for(ait=s.begin(),bit=b.s.begin();
				ait!=s.end() && bit!=b.s.end();
				++ait,++bit)
		{
			if(ait->c!=bit->c)return false;

			if(ait->f.size()!=bit->f.size())return false;

			for(int i=0;i<ait->f.size();++i)
				if(ait->f[i].p!=bit->f[i].p 
					|| ait->f[i].i!=bit->f[i].i )
						return false;						
		}
		
		return ait==s.end() && bit==b.s.end();
	*/
	}
  
	bool operator != (const polyn&b) const
	{ return !(*this==b);	}

	void removez(){
	//remove zeros:
		for(list<summand>::iterator i=s.begin();i!=s.end();){
        while(i!=s.end() && i->c == 0 )
          i=s.erase(i);
        while(i!=s.end() && i->c != 0 )
          ++i;
    }
	}

  polyn& operator += (const polyn&b)
  {
		assert(checkSorted());
		assert(b.checkSorted());
    if(b.s.empty())
      return *this;

    list<summand>::iterator 
      i = s.begin(),iend = s.end();
    list<summand>::const_iterator 
      J = b.s.begin(),Jend = b.s.end();

    while( i!=iend ){
      
			while(J!=Jend && ltf()(*J,*i) )
				s.insert(i,*J),++J;

      while(J!=Jend && veq(J->f,i->f) )
        i->c += J->c , ++J;

      if( J==Jend ) break;

      while(i!=iend && ltf()(*i,*J))
        ++i;    
    }

    s.insert(i,J,Jend);


		removez();

		
		assert(checkSorted());
	
		return *this;

  }


  //negate polynome
  polyn& neg()
  {
    for(list<summand>::iterator i=s.begin();
        i!=s.end();
        ++i)
        i->c=-i->c;         
    return *this;
  }

  inline polyn operator - () const 
  { polyn o=*this;o.neg();return o; }

  polyn& operator -= (const polyn&b)
  {

    neg();
    (*this)+=b;
		neg();
    return *this;

  }



	bool checkSorted() const
	{
		if(s.empty())return true;
		list<summand>::const_iterator i,ii;
		i=s.begin();ii=i;++ii;
		while(ii!=s.end())
		{
			bool lt=*i<*ii;
			assert(lt);
			if(!lt)return false;
			i=ii;		
			++ii;
		}
		return true;
	}



  polyn& operator *= (const polyn&b)
  {
    if(b.s.empty() || s.empty())
    { s.clear(); return *this; }

    list<summand>::iterator ii,i;
    list<summand>::const_iterator J;


    list<summand> ret;

    for(i=s.begin();i!=s.end();++i)
      for(J=b.s.begin();J!=b.s.end();++J)
        ret.push_back(*i),
        ret.back() *= *J;

    ret.sort(ltf());
    ret.swap(s);

    ii=i=s.begin();
    for(++ii; ii!=s.end(); ++ii){
      while(ii!=s.end() && veq(i->f,ii->f))
        i->c += ii->c, ++ii;      
      ++i;
      i = s.erase(i,ii);
      if(i==s.end())
        break;
    }


		removez();


		assert(checkSorted());
    return *this;
  }

	polyn operator * (const polyn&b) const
	{
		polyn ret=*this;
		ret*=b;
		return ret;
	}

	polyn operator + (const polyn&b) const
	{
		polyn ret=*this;
		ret+=b;
		return ret;
	}

	polyn& operator *= (double i)
	{
		list<summand>::iterator sit;
		for(sit=s.begin();sit!=s.end();++sit)
			sit->c *= i;
		return *this;
	}



	template<class T>
  T eval(const vector<T> & values) const
  {
    T ret=0;
    list<summand>::const_iterator i;
    vector<factor>::const_iterator j;
    for(i=s.begin();i!=s.end();++i){
      T x=i->c;
      for(j=i->f.begin();j!=i->f.end();++j)
        for(int k=0;k<j->p;++k)
          x*=values[j->i];
      ret+=x;
    }
    return ret;     
  }

  //calculate value of derivative respective to parameter id
	template<class T>
  T evalderiv(const vector<T> & values,int id) const
  {
    list<summand>::const_iterator i;
    vector<factor>::const_iterator j;
    T ret=(T)0.0;
    for(i=s.begin();i!=s.end();++i){
      T x=(T)i->c;
      for(j=i->f.begin(); j!=i->f.end() && j->i !=id ; ++j)
        for(int k=j->p;k;--k) 
          x*= values[j->i];

      if(j!=i->f.end()){  //if there is an entry of id:

        x*= j->p;

        for(int k=j->p-1;k;--k) 
          x*= values[id];

        for(++j; j!=i->f.end(); ++j)
          for(int k=j->p;k;--k) 
            x*= values[j->i];

        ret+= x;
      }

    }
    return ret;     
  }


	polyn operator /(const polyn &x) const
	{
		assert( x.s.end() == (++s.begin()) );
		assert( x.s.front().f.empty() );
		polyn ret=*this;
		ret *= 1.0/x.s.front().c;
		return ret;
	}
		


  //evaluate integral over variable id. the lower bound is given as values[id].
  template<class T>
  T evalint(const vector<T> &values,int id,const T&upper) const 
  {
   	T ret=(T)0.0;
   	T ret2=ret;
   	const T&lower=values[id];
  	list<summand>::const_iterator i;
  	for(i=s.begin();i!=s.end();++i)
  	{
  		T x=i->c;
		vector<factor>::const_iterator j;
		bool factorinside=false;
		for(j=i->f.begin();j!=i->f.end();++j){
		    	
		    if( j->i==id )
		    {
		    	factorinside=true;
		    	T up=upper;
		    	T lo=lower;
		    	for(int k=0;k<j->p;++k)
			    up*=upper,lo*=lower;
			up-=lo;
		    	x*=up/((T)(j->p+1));	
		    	
		    }
		    else 
		    {
			const T & vji=values[j->i];			
		    	for(int k=0;k<j->p;++k)
		    		x*=vji;
		    }		    			 
		}
		if(factorinside)
		   ret+=x;
		else
		   ret2+=x;
		   				 
  	}
  
  	ret2*=(upper-lower);
  	ret+=ret2;
  	return ret;
  
  }

	//should the ostream operator expand the powers ?
  static bool expandPowers;

	friend ostream&operator<<(ostream&o,const polyn&p);
	
	/*
	print polynomial 
	with common coefficients factored out
	and parameters replaced with paramnames;
	used for repeatTensor output.
	*/
	void printOptCommonCoeff
		( ostream&o,const std::vector<string> &paramnames ) const;

};


inline polyn operator *(double i ,polyn b)
{
	b*=i;	return b;
}

inline polyn operator -(polyn a ,const polyn& b)
{
	a-=b;	return a;
}









#endif

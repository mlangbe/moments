/*
	operations on tensors ignoring symmetry.

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

#ifndef TENSOR_H
#define TENSOR_H

#ifndef GNUCC
#pragma once
#endif
#include<iostream>
#include<vector>
#include<map>
#include<algorithm>
#include<cassert>
#include "vvpr.h"
#include "refp.h"

//abstract base class

//a tensor with only ord, dim defined
class odtens: public refCountable
{
protected:
 int o,d;

public:
  odtens(int order,int dim)
    :o(order),d(dim)
  {}
  
  odtens(const odtens&odt)
    :o(odt.o),d(odt.d)
  {}

  int ord() const {return o;}
  int dim() const {return d;}
};

/**
*iterates through the tensor indices 
*in an order as given in tensor<T>
*/
template<class visitor>
class visitTensorIndices :public odtens
{
	visitor& f;

	int * ind;
	int depth;

public:

	/** invokes f for every possible index
	* of a tensor with order ord and dimension d
	*/
	visitTensorIndices(const odtens&odt,visitor&ff)
		:odtens(odt),f(ff)
	{
		if(o>0)
		{
			ind=new int[o];
			std::fill(ind,ind+o,0);
			depth=0;
			recurse();
		}
		else
		{
			ind=0;
			f(0);
		}
	}

	~visitTensorIndices(){if(ind)delete[]ind;}

private:

	void recurse()
	{		
		if(depth<o-1)
		{
			for(ind[depth]=0;ind[depth]<d;++ind[depth])
			{
				++depth; 
				recurse(); 
				--depth;
			}       
		}
		else
		{
			for(ind[depth]=0;ind[depth]<d;++ind[depth])
				f(ind);		
		}
	}
};




/** a constant tensor */
template<class T>
class tens:public odtens{

public:
  tens(int order,int dim)
    :odtens(order,dim)
  {}
  
  tens(const odtens&odt)
    :odtens(odt)
  {}

  virtual T operator[](const int * inds) const=0;

	virtual ~tens(){}
  
	virtual int getNumDistinct() const
	{ assert(false/*not implemented here*/);return -1;}
private:
	
	struct equalityVisitor
	{
		const tens<T> &a, &b;
		bool isequal;

		equalityVisitor(const tens<T>&aa,const tens<T>&bb)
			:a(aa),b(bb)
		{ isequal=true; }
	

		void operator() (const int*inds)
		{
			if(isequal){
				isequal &= a[inds]==b[inds];
				/*
				if(!isequal)
				{
					cout<<a[inds]<<"\n\n!=\n\n"<<b[inds]<<endl;
				}
				*/
			}

		}
	
	};

	struct tensPrintVisitor
	{
		const tens<T>&a;
		std::ostream&out;
		tensPrintVisitor(std::ostream&o,const tens<T> &aa)
			:out(o),a(aa)
		{
		}

		void operator() (const int*inds)
		{
			for(int i=0;i<a.ord();++i)
				out<<inds[i]<<"=";
			out<<a[inds]<<'\n';
		}

	};

public:



	bool operator == (const tens<T>&x) const
	{
		if( x.o!=o | x.d!=d )return false;
		equalityVisitor v(*this,x);
		visitTensorIndices<equalityVisitor> vt(x,v);
		return v.isequal;
	}


	
	template<class S>
	friend std::ostream&operator<<(std::ostream&o,const tens<S> & x);


};
 

template<class T>
std::ostream&operator<<(std::ostream&o,const tens<T> & x)
{
	typename tens<T>::tensPrintVisitor v(o,x);
	visitTensorIndices< typename tens<T>::tensPrintVisitor> vt(x,v);
	return o;
}



/** a writable tensor */
template<class T>
class wtens:public tens<T>{
public:
  wtens(int order,int dim)
    :tens<T>(order,dim)
  {}

  wtens(const odtens&odt)
    :tens<T>(odt)
  {}
	
	virtual T & operator[](const int* inds)=0;

	virtual ~wtens(){}

	virtual T operator[](const int * inds) const
  {
    return const_cast<wtens&>(*this)[inds];    
  }  
  
};

template<class T>
class tensProd: public tens<T>
{
  const tens<T> &x,&y ;
public:
  tensProd(const tens<T>&xx,const tens<T>&yy)
    :x(xx),y(yy),tens<T>(xx.ord()+yy.ord(),xx.dim())
  {}

	virtual ~tensProd(){}

  T operator[](const int * inds) const
  {
    return x[inds] * y[inds+x.ord()];
  }
};

	

//a concrete implementation of tensors  
template<class T>
class tensor : public wtens<T>{
public:
#ifndef OLD_VC
  using odtens::o;
  using odtens::ord;
  using odtens::d;
 #endif

  typename std::vector<T> t;

  tensor(int order,int dim)
    :wtens<T>(order,dim)
  {

    int s=1;
    for(;order>0;--order)
      s*=dim;

    t.resize(s);    
  }

  tensor(const tensor& to)
    :wtens<T>(to)
  {
    t=to.t;  
  }

    

	
	struct copyVisitor
	{
      typename std::vector<T>::iterator a;
      const tens<T> &xx;						
			
			inline void operator()(int*indices)
			{ *(a++) = xx[indices];  }

			copyVisitor(const tens<T>&xxx, tensor<T>&x0)
				:xx(xxx),a(x0.t.begin())
			{}
	};


	tensor()
		:wtens<T>(0,0)
	{		
	}


	void operator=(const tens<T> & to)
	{
		copyFromTens(to);
	}

	void operator=(const tensor<T> & to)
	{
		o=to.ord();d=to.dim();
		t=to.t;
	}


  tensor(const tens<T>& to)
    :wtens<T>(to.ord(),to.dim())
  {
		copyFromTens(to);
	}

	void copyFromTens(const tens<T>& to)
	{
		o=to.ord();d=to.dim();
    
		int s=1;
    for(int i=ord();i>0;--i)
      s*=d;

    t.resize(s);    

		copyVisitor v(to,*this);
		visitTensorIndices<copyVisitor> vti(to,v);  		
	}

	virtual ~tensor(){}

  inline int getidx(const int* inds) const
  {
    int id=0;
    for(int i=0;i<o;i++)
      id*=d, id+=inds[i];  
		return id;
  }

  inline T & get(const int* inds)
  {
    return t[getidx(inds)];    
  }

  inline const T & get(const int* inds) const
  {
    return t[getidx(inds)];    
  }

  T & operator[](const int* inds)
  { return get(inds); }

  tensor& operator*=(const tensor& x)
  {
    typename std::vector<T> ret(t.size()*x.t.size());
    typename std::vector<T>::const_iterator i,j;
    typename std::vector<T>::iterator  retit=ret.begin();

    for(i=t.begin();i!=t.end();++i)
      for(j=x.t.begin();j!=x.t.end();++j)
        *retit= *i * *j , ++retit;
    ret.swap(t);
    o+=x.ord();
    return *this;
  }

  inline tensor operator *(const tensor&o)
  {
    tensor x=*this;
    x*=o;
    return x;
  }


  int getNumDistinct() const
  {
    typename std::vector<T> s=t;
    std::sort(s.begin(),s.end());
    int ret=0;
    for(int i=0;i<(int)s.size();++i,++ret)
      while( i+1<(int)s.size() && s[i]==s[i+1])
        ++i;
    
    return ret;
  }
  
  


};



class inditer
{
  int*ind;
  int d,o;

public:

  inditer(int dim, int ord)
    :d(dim),o(ord)
  { 
    ind=new int[o]; 
    std::fill(ind,ind+o,0); 
  }

  inditer(const inditer&x)
    :d(x.d),o(x.o)
  {
    ind=new int[o];
    std::copy(x.ind,x.ind+o,ind);
  }

  ~inditer()
  {delete[]ind;}

  const int* operator * ()
  {return ind;}


  inditer& operator ++ ()
  {
    int i=0;
    while( i<o && ind[i]==d-1)
      ind[i]=0,++i;

    if(i<o)
      ++ind[i];

    return *this;
  }

  inditer& operator -- ()
  {
    int i=0;
    while( i<o && ind[i]==0)
      ind[i]=d-1,++i;

    if(i<o)
      --ind[i];

    return *this;
  }

  bool operator == (const inditer&x)
  {
    assert(o==x.o && d==x.d );
    for(int i=0;i<o;++i)
      if(ind[i] != x.ind[i])
          return false;

    return true;
  }

  bool operator == (int x)
  {
    for(int i=0;i<o;++i)
      if(ind[i]!=x)
        return false;

    return true;
  }

  bool operator != (int x)
  {
    return ! (*this == x);
  }

  bool operator != (const inditer& x)
  {
    return ! (*this == x);
  }

};



template<class T>
class foldBase
{

protected:
  mutable int depth;

  const tens<T>&a;
  const vpr & indpairs ;

public:
  static const int maxOrd=32;
protected:
  mutable int ind[maxOrd];
  mutable T ret;

  foldBase(const tens<T>&a,
		const vpr & indpairs)
		:a(a),indpairs(indpairs)
  {		
    assert(a.ord()<=maxOrd);
  }
public:

  operator T () const
  { return ret; }

protected:
  void recurse() const     
  {                
     for(int i=0;i<a.dim();++i){
        ind[ indpairs[depth].a ] = i;
        ind[ indpairs[depth].b ] = i;
        if(depth){                   
           --depth; recurse(); ++depth;
        }        
        else{
          /*
          for(int i=0;i<a.ord();++i)
            std::cout<<ind[i]<<' '<<std::endl;*/
          ret += a[ind];
          //std::cout<<"bla"<<endl;
        }
     }
  }
  
};

template<class T>
const int foldBase<T>::maxOrd;

template<class T>
class fold :public foldBase<T>
{
#ifndef OLD_VC
	using foldBase<T>::depth;
	using foldBase<T>::ind;
	using foldBase<T>::indpairs;
	using foldBase<T>::maxOrd;
	using foldBase<T>::a;
	using foldBase<T>::recurse;
#endif
public:
  fold(const tens<T>&a,
		const vpr & indpairs)    
		:	foldBase<T>(a,indpairs)
  {		

		assert(sanitycheck());

		depth=int(indpairs.size()-1);
    recurse();
  }
private:
	bool sanitycheck()
	{
		/** first:some sanity checks:*/
    std::fill(ind,ind+maxOrd,-1);
		vpr::const_iterator it;

		for(it=indpairs.begin();it!=indpairs.end();++it)
		{
			assert( it->a < a.ord() );
			assert( it->b < a.ord() );
			//occupy indices
			ind[it->a]=-2;
			ind[it->b]=-2;
		}

		assert(indpairs.size()*2 == a.ord());

		//check if all indices are occupied
		for(int i=0;i<a.ord();++i)
			assert(ind[i]==-2);
	
		return true;
	}

};


template<class T>
class partiallyFolder:public foldBase<T>
{
  //the index positions of the free indices
protected:
  const int* freeInds;
	int myord;

public:
#ifndef OLD_VC
	using foldBase<T>::depth;
	using foldBase<T>::recurse;
	using foldBase<T>::indpairs;
	using foldBase<T>::ind;
	using foldBase<T>::ret;
#endif



  T operator [] (const int *index) const
  { 	
		ret = T(0);
		//set the free indices
		for(int i=0;i<myord;++i)
			ind[freeInds[i]]=index[i];

		//and sum over the paired ones.
    depth=int(indpairs.size()-1);
    recurse();			
		return ret; 
	}

  partiallyFolder(const int*freeInds,const tens<T>&a,const vpr & indpairs)
  :foldBase<T>(a,indpairs),freeInds(freeInds)
  {	
		myord = a.ord()-indpairs.size()*2;
  }	

};



template<class T>
class partialFold : public tens<T>,public partiallyFolder<T>
{

	//gives the index indices of the free indices
	int freeInds_[foldBase<T>::maxOrd];

public:
#ifndef OLD_VC
	using foldBase<T>::depth;
	using foldBase<T>::indpairs;
	using foldBase<T>::a;
#endif


  partialFold(const tens<T>&a_,
		const vpr & indpairs_)
		:tens<T>(a_.ord()-indpairs_.size()*2,a_.dim()),
		partiallyFolder<T>(freeInds_,a_,indpairs_)
		//a(a_),indpairs(indpairs_)
  {	
    assert(a.ord() < foldBase<T>::maxOrd );
		
		int ind[foldBase<T>::maxOrd];

		std::fill(ind,ind+a.ord(),-1);
		vpr::const_iterator it;

		for(it=indpairs.begin();it!=indpairs.end();++it)
		{
			assert( it->a < a.ord() );
			assert( it->b < a.ord() );
			//occupy indices
			ind[it->a]=-2;
			ind[it->b]=-2;
		}

		int freeIndId=0;
		for(int i=0;i<a.ord();++i)
			if(ind[i] != -2)
				freeInds_[freeIndId++]=i;
  }


  T operator [] (const int *index) const
  {
		return partiallyFolder<T>::operator[](index);
  }

	virtual ~partialFold(){}
  
};


/**
 a tensor wth indices permuted.
*/
template<class T>
class indexPerm : public tens<T>
{
	mutable std::vector<int> buf;

public:
	const std::vector<int> & perm;
	const tens<T>& a;

  T operator [] (const int *index) const
  {
		for(int i=0;i<a.ord();++i)
			buf[perm[i]] = index[i];
    return a[&buf[0]];  
  }

	indexPerm(const tens<T>&a_,const std::vector<int> & perm_)
		:tens<T>(a_),perm(perm_),a(a_),buf(a_.ord())
	{}

	virtual ~indexPerm(){}
	
};



template<class T>
class basisTransformed : public tens<T>
{

	struct visitor{
		
		const tens<T>&t;
		const T * const M;

		T sum;
		const int * index;

		
		void operator() (const int* inds)
		{
			T s= t[inds];
			for(int i=0;i<t.ord();++i)
				s *= M[ inds[i] * t.dim() + index[i]  ];
			sum+=s;
		}

		visitor(const tens<T> & tt,const T * MM)
			:t(tt),M(MM)
		{}
	};


	const tens<T>&t;
	const T * const M;


public:

  T operator [] (const int *index) const
  {
		visitor v(t,M);
		v.sum=0;
		v.index=index;
		visitTensorIndices<visitor> vv(t,v);
		return v.sum;
	}

	basisTransformed(const tens<T>&a,const T *mat )
		:tens<T>(a),t(a),M(mat)
	{}	

	~basisTransformed(){}
	
};




#endif

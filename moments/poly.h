/*
	polynomial representation with coefficients and exponents

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

#ifndef POLY_H
#define POLY_H

#include "polynom.h"

#ifndef OLD_VC
#define _NOTYPENAME_
#else
#define _NOTYPENAME_ typename
#endif

//a memory-saving and also more efficient way of storing the polynomes
//(if not too many variables are inolved)
template<int maxVar,class Tcoeffs=float>
struct poly{
	
	static const int n=maxVar;
	typedef Tcoeffs tcoeff;
	static const int maxVarInt = (n+sizeof(int)-1)/sizeof(int);
	struct summand{	
		Tcoeffs c;
		
		union powers{
			unsigned char p[n];
			unsigned int i[maxVarInt];

			bool operator==(const powers&s) const
			{
				for(int j=0;j<maxVarInt;++j)
					if(i[j]!=s.i[j])
						return false;
				return true;
			}
		}p;

		bool operator<(const summand&s) const
		{
			int i=0;
			while(i<n && p.p[i]==s.p.p[i])++i;			
			return i<n && p.p[i]<s.p.p[i]  || i==n && c<s.c ;
		}

		bool operator==(const summand&s) const
		{			
			return p==s.p && c==s.c;			
		}

	};



	_NOTYPENAME_ vector<summand> s;


	poly(){}


	polyn topolyn() const
	{
		polyn ret;
		typename polyn::summand ps;
		typename polyn::factor pf;
		typename vector<summand>::const_iterator i; 
		for(i=s.begin();i!=s.end();++i)
		{
			ps.c=(double)i->c;ps.f.clear();
			for(int j=0;j<maxVar;++j)
			{
				if(i->p.p[j])
				{
					pf.i=j;
					pf.p=i->p.p[j];
					ps.f.push_back(pf);
				}
			}
			ret.s.push_back(ps);
		}
		return ret;
	}



	void getVarIdsPresent(bool ids[maxVar]) const
	{
		summand buf;
		typename vector<summand>::const_iterator i;
		fill(buf.p.i,buf.p.i+maxVarInt,0);
		
		for(i=s.begin();i!=s.end();++i)
			for(int j=0;j<maxVarInt;++j)
				buf.p.i[j] |= i->p.i[j];

		for(int j=0;j<maxVar;++j)
			ids[j]= buf.p.p[j]!=0;
	}


	explicit poly(const polyn&p)
	{
		s.resize(p.s.size());
		typename std::vector<summand>::iterator sit;
		typename std::list<polyn::summand>::const_iterator psit=p.s.begin();

		for(sit=s.begin();sit!=s.end();++sit,++psit)
		{
			vector<polyn::factor>::const_iterator vit;
			fill(sit->p.i,sit->p.i+maxVarInt,0);
			for(vit=psit->f.begin();vit!=psit->f.end();++vit)
			{
				if( vit->i >= maxVar || vit->p > 0xff )
					throw;

				sit->p.p[vit->i] = vit->p;		
			}
			sit->c=(tcoeff)psit->c;
		}	
		sort(s.begin(),s.end());
	}


	template< class TT >
	explicit poly(const poly<maxVar,TT> &x)
	{
		s.resize(x.s.size());
		for(int i=0;i<x.s.size();++i){
			s[i].c=x.s[i].c;
			std::copy(x.s[i].p.i,x.s[i].p.i+maxVarInt,
									s[i].p.i);
								
		}
	
	
	}


	poly(const Tcoeffs& coeff)		
	{
		if(coeff!=0){
			s.resize(1);
			_NOTYPENAME_ summand&su = s[0];
			su.c=coeff;
			std::fill(su.p.i,su.p.i+maxVarInt,0);
		}

	}

	poly(const Tcoeffs& coeff,int id)
	{ 
		if(coeff!=0){
			s.resize(1);
			_NOTYPENAME_ summand&su = s[0];
			su.c=coeff;
			std::fill(su.p.i,su.p.i+maxVarInt,0);
			su.p.p[id]=1;
		}
	}

	bool operator == (const poly&a) const
	{
		return a.s==s;
	}
	inline bool operator != (const poly&a) const
	{
		return!(a.s==s);
	}

	//terme mit gleichen variablen-potenzen zusammenfassen
	//precondition: s is sorted
	void compact()
	{
		if(s.empty())return;
		typename vector<summand>::iterator i,ii;

		checksorted();

		ii=i=s.begin();
		++ii;		
		while(ii!=s.end() && i->p == ii->p)
			i->c += ii->c, ++ii;      
		
		while(ii!=s.end())
		{
			if(i->c!=Tcoeffs(0))++i;
			*i=*ii;
			++ii;
			while(ii!=s.end() && i->p == ii->p)
		     i->c += ii->c, ++ii;      
		}
		if(i->c!=Tcoeffs(0))++i;
		s.resize(i-s.begin());
	}


	void mul(const poly&b,poly&ret) const
	{
		if(s.size()==0||b.s.size()==0) {ret.s.clear();return;}

		ret.s.resize(b.s.size()*s.size());
		typename vector<summand>::iterator r;
		typename vector<summand>::const_iterator i,j;

		r = ret.s.begin();
		for(i=s.begin();i!=s.end();++i)
			for(j=b.s.begin();j!=b.s.end();++j,++r)
			{
				//multiply coeffs
				r->c=i->c;
				r->c*=j->c;
				//add powers 
				//(for efficiency, the integers
				// and not the char powers are added here
				//(so 4 char powers are added at once)
				for(int k=0;k<maxVarInt;++k)
					r->p.i[k] = i->p.i[k]+j->p.i[k];						
			}				
		sort(ret.s.begin(),ret.s.end());
		ret.compact();
  }



  
	poly operator * (const poly&b) const
  {
		poly ret;
		mul(b,ret);
		return ret;
	}

	poly& operator *= (const poly&b)
  {
		poly ret;
		mul(b,ret);
		s.swap(ret.s);
		return *this;
	}

	poly& operator *= (Tcoeffs v)
	{
		typename vector<summand>::iterator r;
		for(r=s.begin();r!=s.end();++r)
			r->c *=v;			
		return *this;
	}

	poly& operator /= (Tcoeffs v)
	{
		typename vector<summand>::iterator r;
		for(r=s.begin();r!=s.end();++r)
			r->c /= v;			
		return *this;
	}


	void add (const poly&b,poly&ret) const
  {		
		ret.s.resize(s.size()+b.s.size());
		merge(s.begin(),s.end(),b.s.begin(),b.s.end(),ret.s.begin());
		ret.compact();
  }


	poly& operator += (const poly&b)
	{
		size_t osiz=s.size();
		s.resize(s.size()+b.s.size());
		copy(b.s.begin(),b.s.end(),s.begin()+osiz);
		inplace_merge(s.begin(),s.begin()+osiz,s.end());
		compact();
		return *this;
	}


	void checksorted() const
	{
		for(int i=1;i<s.size();++i)
			assert( !(s[i]<s[i-1]) );
	}


	poly& operator -= (const poly&b)
	{
		checksorted();b.checksorted();
		size_t osiz=s.size();
		s.resize(s.size()+b.s.size());
		copy(b.s.begin(),b.s.end(),s.begin()+osiz);
		
		//negate the b part
		typename vector<summand>::iterator r;		
		for(r=s.begin()+osiz;r!=s.end();++r)
			r->c=-r->c;			
		
		inplace_merge(s.begin(),s.begin()+osiz,s.end());
		compact();
		return *this;
	}

	void neg()
	{
		typename vector<summand>::iterator r;		
		for(r=s.begin();r!=s.end();++r)
			r->c =-r.c;			
	}

	inline poly operator -() const
	{
		poly ret=*this; ret.neg(); return ret;
	}


	inline poly operator - (const poly&b) const
	{
		poly ret=*this;
		ret-=b;
		return ret;
	}

	inline poly operator + (const poly&b) const
	{
		poly ret;
		add(b,ret);
		return ret;
	}


	void differentiate(int i)
	{
		assert(i<maxVar);	
		typename vector<summand>::iterator r,rr;

		for(r=s.begin(),rr=s.begin();r!=s.end();++r)
		{
			if(r->p.p[i]){
				r->c *= (Tcoeffs)r->p.p[i];
				--r->p.p[i];
				*rr=*r;
				++rr;
			}			
		}
		s.resize(rr-s.begin());		
	}

	void integrate(int i)
	{
		assert(i<maxVar);	
		typename vector<summand>::iterator r,rr;
		for(r=s.begin(),rr=s.begin();r!=s.end();++r)
		{
			++r->p.p[i];
			r->c /= r->p.p[i];
		}				
	}

	void integrate(int index,Tcoeffs lower,Tcoeffs upper)
	{
		assert(index<maxVar);	
		typename vector<summand>::iterator r,rr;
		Tcoeffs l,u;
		for(r=s.begin(),rr=s.begin();r!=s.end();++r)
		{
			int pw = r->p.p[index];
			++pw;

			l=lower;u=upper;
			for(int j=1;j<pw;++j)
				l*=lower,u*=upper;
			r->c *= u-l;
			r->c /= (Tcoeffs)pw;
			r->p.p[index]=0;
		}				

		sort(s.begin(),s.end());
		compact();
	}


	template<class T>
	void eval(const vector<T> values,T&ret) const
	{
		if(s.empty())
		{	ret=T(); return;}
			
		assert(values.size()>=maxVar);
		typename vector<summand>::const_iterator r;
		T buf;

		r=s.begin();
		ret=r->c;
		for(int i=0;i<maxVar;++i)
			for(int j=0;j<r->p.p[i];++j)
				ret*=values[i];

		for(++r;r!=s.end();++r)
		{
			buf=r->c;
			for(int i=0;i<maxVar;++i)
				for(int j=0;j<r->p.p[i];++j)
					buf*=values[i];
			ret+=buf;
		}				
	}

	template<class T>
	inline T eval(const vector<T> values) const
	{
		T ret;
		eval(values,ret);
		return ret;
	}



	void save(ostream&o)
	{
		o<<s.size()<<'\n';
		typename vector<summand>::iterator it;
		for(it=s.begin();it!=s.end();++it)
		{
			o<<it->c<<' ';
			for(int j=0;j<maxVar;++j)
				o<<((int)it->p.p[j])<<' ';
			o<<'\n';
		}		
	}

	void load(istream&i)
	{
		summand buf;
		s.clear();
		size_t siz;
		i>>siz;
		s.resize(siz);
		typename vector<summand>::iterator it;
		for(it=s.begin();it!=s.end();++it)
		{
			i>>it->c;
			int buf;
			for(int j=0;j<maxVar;++j)
			{i>> buf;it->p.p[j]=(unsigned char)buf;}
		}	
	}



};


template<int numVar,class Tcoeff>
ostream&print(std::ostream&o,const poly<numVar,Tcoeff> & p )
{
	typename vector<typename poly<numVar,Tcoeff>::summand>::const_iterator sit;

	o<<
		"calcIntegralPolys(const T*A)\n"
		"{\n"
		"T ret=0; \n";

	for(sit=p.s.begin();sit!=p.s.end();++sit)
	{
		o<<"ret+=";
		int i;
		for(i=0;i<numVar;++i)
			if(sit->p.p[i])break;

		if(i==numVar)
			o<<sit->c;
		else
		{
			if(sit->c==-1)
				o<<"-";
			else if(sit->c!=1)
				o<<sit->c<<"*";

			bool first=true;

			for(;i<numVar;++i)
			{
				for(int j=0;j<sit->p.p[i];++j)
				{
					if(!first)
					{	o<<"*";}
					first=false;
					o<<"A["<<i<<"]";
				}
			}

		}
		o<<";\n";

	}
	o<<"return ret;";
	return o;

}


template<int numVar,class Tcoeff>
ostream&operator<<(std::ostream&o,const poly<numVar,Tcoeff> & p )
{
	typename vector<typename poly<numVar,Tcoeff>::summand>::const_iterator sit;

	o<<" ( ";
	for(sit=p.s.begin();sit!=p.s.end();++sit)
	{
		if(sit!=p.s.begin())
			o<<" + ";

		int i;
		for(i=0;i<numVar;++i)
			if(sit->p.p[i])break;

		if(i==numVar)
			o<<sit->c;
		else
		{
			if(sit->c==(Tcoeff)-1)
				o<<"-";
			else if(sit->c!=(Tcoeff)1)
				o<<sit->c<<"*";

			bool first=true;

			for(;i<numVar;++i)
			{
				for(int j=0;j<sit->p.p[i];++j)
				{
					if(!first)
					{	o<<"*";}
					first=false;
					o<<"A["<<i<<"]";
				}
			}

		}

	}
	o<<" ) ";
	return o;
}




















#endif
#undef _NOTYPENAME_

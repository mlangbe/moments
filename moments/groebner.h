/*
	operations on  sets of polynomials used in buchberger algorithm to create groebner bases.

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
#ifndef GROEBNER_H
#define GROEBNER_H

#include<vector>
#include<iostream>

template<class P, class Su>
struct groebner_spec
{

	typedef P tpoly;
	typedef Su tsummand;

	/** helper methods that have to be defined for polnomial class P
	to be used in buchberger
	algo
	*/
	/** remove common divisors from a,b*/
	static void removeGCD(Su&a,Su &b);

	static void normalize(P&ret){
		ret/=LT(ret).c;
	}

	static bool isConstant(const Su &p);

	/** get summand highest in order*/
	static const Su& LT(const P &p);

		/** remove LT from e */
	static void removeLT(P&e)
	{
		e.s.pop_back();
	}


	/** remove LT from e and add it to rest*/
	static void transferLT(P&e,P&rest);

	static void setZero(P &p);
	static  bool isZero(const P &p);

	/** wrap summand in a single-summand polynomial*/
	static P& wrapSummand(P&ret,const Su &p);
};




template<class P, class Su=typename P::summand>
struct groebner:groebner_spec<P,Su>
{

#ifndef OLD_VC
	using groebner_spec<P,Su>::isZero;
	using groebner_spec<P,Su>::setZero;
	//using typename groebner_spec<P,Su>::tpoly;
	//using typename groebner_spec<P,Su>::tsummand;

	/** helper methods that have to be defined for polnomial class P
	to be used in buchberger
	algo
	*/
	/** remove common divisors from a,b*/
	using groebner_spec<P,Su>::removeGCD;

	using groebner_spec<P,Su>::normalize;

	using groebner_spec<P,Su>::isConstant;

	/** get summand highest in order*/
	using groebner_spec<P,Su>::LT;

		/** remove LT from e */
	using groebner_spec<P,Su>::removeLT;

	/** remove LT from e and add it to rest*/
	using groebner_spec<P,Su>::transferLT;


	/** wrap summand in a single-summand polynomial*/
	using groebner_spec<P,Su>::wrapSummand;

#else
#endif

	typedef P tpoly;
	typedef Su tsummand;


	///////////////////////////////////////////////////////////////////
	/*
	static void divide(P&result,P&e,const P&d);
	static void divide(std::vector<P >&result,P &e,const std::vector< P >&d);
	static void S(P&ret,const P&a,const P&b);
	static void buchberger(std::vector<P>&ret);
	*/

static void divide(P&result,P&e,const P&d)
{
	setZero(result);

	if(isZero(d))
		return;

	const tsummand &ltd= LT(d);
	

	tsummand bufs,lte,myltd;
	P buf,buf2;

	while(!isZero(e)){
		
		lte= LT(e);
		myltd = ltd;
		removeGCD(lte,myltd);
		//if ltd not contained in lte
		if(!isConstant(myltd))
		{
			break;
		}
		lte.c/=myltd.c;
		buf2=wrapSummand(buf,lte)*d;
		e-=buf2;
		result+=buf;
	}
}



static void reduce(std::vector<P>&ret)
{
	for(int i=ret.size()-1;i>=0;--i)
		normalize(ret[i]);

	std::vector<P> lts(ret.size());
	P buf;
	for(int i=0;i<ret.size()-1;--i)
	{
		wrapSummand(lts[i],LT(ret[i]));
	}
	

	for(int i=ret.size()-1;i>=0;--i)
	{
		buf=ret[i];
		remainder(buf,lts);
		if(isZero(buf))
		{
			setZero(ret[i]);	
			if(i>0)
				setZero(lts[i-1]);
		}
		else if(i>0) //add this lt to ltd and remove preceding one.
		{
			wrapSummand(lts[i-1],LT(ret[i]));
		}
	}

	ret.resize( remove_if(ret.begin(),ret.end(),isZero)-ret.begin() );

}

static void reduce2(std::vector<P>&ret)
{

	using namespace std;

	for(int i=ret.size()-1;i>=0;--i)
		normalize(ret[i]);

	P buf,oret;
	for(int i=0;i<ret.size();++i)
	{
	
		oret = ret[i];
		buf=ret[i];
	
		setZero(ret[i]);
		//ret[i] is now 0

		remainder(buf,ret);

		if(isZero(buf))
		{
			cout<<i<<"th term redundant."<<endl;
		}
		else{
			ret[i]=oret;
		}
	}
	//remove terms expressible by others
	ret.resize( remove_if(ret.begin(),ret.end(),isZero)-ret.begin() );

}



/** calculate remainder. return true if e changed*/
static bool remainder(P&e,P d)
{
	if(isZero(d))
		return false;

	const tsummand ltd= LT(d);
	removeLT(d);

	tsummand bufs,lte,myltd;
	P buf,buf2;

	bool change=false;
	while(!isZero(e)){
		
		lte= LT(e);
		myltd = ltd;
		removeGCD(lte,myltd);
		//if ltd not contained in lte
		if(!isConstant(myltd))
		{
			break;
		}
		removeLT(e);	

		lte.c/=myltd.c;
		buf2=wrapSummand(buf,lte)*d;

		e-= buf2;
		change=true;
	}

	return change;
}






static void divide(std::vector<P >&result,P &e,const std::vector< P >&d)
{
	P buf,rest;
	static const P ZERO;
	
	result.clear();
	result.resize(d.size(),ZERO);

	while(!e.s.empty())
	{
		bool hasdiv=false;

		for(int i=0;!isZero(e) && i<d.size();++i) 
		{
			divide(buf, e, d[i]); 
			
			if( !isZero(buf))
			{
				result[i]+=buf;
				hasdiv=true;
			}
		}

		if(isZero(e))
			break;
		//remove leading term from f and put it to r.
		if(!hasdiv)
		{
			transferLT(e,rest);	
		}
	}

	rest.s.swap(e.s);
}

static void remainder(P &e,const std::vector< P >&d)
{
	P buf,rest;

	while(!e.s.empty())
	{
		bool hasdiv=false;

		for(int i=0;!isZero(e) && i<d.size();++i) 
		{
			hasdiv |= remainder(e, d[i]); 
		}

		if(isZero(e))
			break;
		//remove leading term from f and put it to r.
		if(!hasdiv)
		{
			transferLT(e,rest);	
		}
	}

	rest.s.swap(e.s);
}



/**
S from buchberger algo
*/

static void S(P&ret,const P&a,const P&b)
{
	setZero(ret);
	if(isZero(a)||isZero(b))
	{
		return;
	}
	Su lca=LT(a),lcb=LT(b);
	removeGCD(lca,lcb);

	P bufa,bufb;

	wrapSummand(bufa,lca);
	wrapSummand(bufb,lcb);

	bufa*=b;
	bufb*=a;

	ret=bufa;
	ret-=bufb;
}


static void buchberger(std::vector<P>&ret)
{
	using namespace std;
	int onum=0,oonum=0;
	P s;
	std::vector<P> dummy;

	int monomials=0;
	for(int i=0;i<ret.size();++i)
	{
		monomials += ret[i].s.size();

	}

	int steps=0;
	int iter=0;
	do{
		reduce2(ret);
		onum=ret.size();
		//only have to test newly available pairs
		int esteps= onum *(onum-1)/2;

		int st0=0;
		for(int i=0;i<onum;++i)for(int j=i+1;j<onum;++j)
		{
			S(s,ret[i],ret[j]);		
			
			//cout<<"before:"<<s<<endl;
			
			remainder(s,ret);		

			if(!isZero(s)){
				//cout<<'\n'<<'+'<<'\n'<<endl;
				ret.push_back(s);
				monomials+=s.s.size();
			}
			
			++st0;
			cout<<"\riteration "<<iter<<" basic steps:"<<st0<<'/'<<esteps<<"  #polynomials:"<<ret.size()<<" #summands"<<monomials<<'\r';
			cout.flush();
		}
		++iter;
		cout<<endl;
	}while(ret.size()>onum);

	cout<<"\nready";
}




};


template<class P>
void buchberger(std::vector<P>&ret)
{
	groebner<P>::buchberger(ret);
}
#endif

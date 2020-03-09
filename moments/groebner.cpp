/*
	operations on  sets of polynomials used in buchberger algorithm.


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

#include "groebner_poly.hpp"
#include "infbruch.h"

#include "bruch.h"

template<>
struct groebner_spec<polyn,polyn::summand>
{
static
void normalize(polyn&a)
{
	float f=a.s.back().c;
	for(list<polyn::summand>::iterator 
		it=a.s.begin();it!=a.s.end();++it)
			it->c/=f;
}


/**
*remove common factors from a,b
*/
static
void removeGCD(polyn::summand&a,polyn::summand&b)
{

	//for every param id, the order
	typedef vector<polyn::factor>::iterator vpit;
	//remove common divisors from f1,f2

	vpit ita=a.f.begin(),itb=b.f.begin(), nita=ita, nitb=itb;
	
	while(itb!=b.f.end())
	{
		while(ita!=a.f.end() && ita->i < itb->i)
		{
			*nita=*ita;
			++ita;
			++nita;
		}
		if(ita==a.f.end())
		{
			break;
		}


		while(itb!=b.f.end() && itb->i < ita->i) 
		{
			*nitb=*itb;
			++itb;
			++nitb;
		}
		if(itb==b.f.end())
		{
			break;
		}


		if(ita->i==itb->i)
		{
			int oii=min( ita->p,itb->p );
			ita->p-=oii;
			itb->p-=oii;

			if(ita->p!=0)
				*(nita++)=*ita;

			if(itb->p!=0)
				*(nitb++)=*itb;

			++ita;
			++itb;
		}
	}

	while(itb!=b.f.end())
		*(nitb++)=*(itb++);

	while(ita!=a.f.end())
		*(nita++)=*(ita++);


	b.f.resize(nitb-b.f.begin());
	a.f.resize(nita-a.f.begin());
}

static bool isConstant(const polyn::summand&p)
{
	return p.f.empty();
}

static bool isZero(const polyn&p)
{
	return p.s.empty();
}

static void  setZero(polyn&p)
{
	return p.s.clear();
}


static const polyn::summand&LT(const polyn&p)
{
	return p.s.back();
}

static polyn& wrapSummand(polyn&ret,polyn::summand p)
{
		ret.s.clear();
		ret.s.push_back(p);
		return ret;	
}

static
void transferLT(polyn&e,polyn&rest)
{
	polyn buf;		
	buf.s.splice(buf.s.begin(),e.s,--e.s.end());	
	rest+=buf;
}

static
void removeLT(polyn&e)
{
	e.s.pop_back();
}


};



namespace gr{

void testdiv()
{
	polyn x(1.0,1),y(1.0,2),z(1.0,3),s(1.0,4),t(1.0,5);

	polynom::ioFormat=polynom::RAW;
	vector<polyn> d,result;
	d.push_back(x*x+y);

	polyn e=x*x*x+y,rem=e;
	
	groebner<polyn>::divide(result,rem,d);

	cout<<e<<" = "<<rem<<" + "<<result[0]<<" * "<<d[0]<<endl;

}


template<class T>
ostream&operator<<(ostream&o,const vector<T>&t)
{
	for(int i=0;i<t.size();++i)
		o<<t[i]<<endl;
	return o;
}

void test2()
{
	typedef polyn tpl;

	
	int id=0;
	tpl s[2]={tpl(1,id++),tpl(1,id++)};
	tpl x[3]={tpl(1,id++),tpl(1,id++),tpl(1,id++)};
	tpl c[30];
	int ii=97;
	for(int i=0;i<18;++i,ii*=97)
		c[i]=tpl(1,id++);


	tpl g[3];

	for(int i=0;i<3;++i)
	{
		//constant
		g[i]=x[i]- c[i];

		//linear
		g[i]-= s[0]*c[3+i]-s[1]*c[6+i];

		//quadratic
		g[i]-= s[0]*s[0]*c[9+i]-s[1]*s[1]*c[12+i] -s[1]*s[0]*c[15+i];
	}



	vector<tpl> res(g,g+3);

	vector<poly<36,infbruch > > res2(g,g+3);


	//vector<poly<36,bruch<long long,(~((long long unsigned)0) >> 1)> > >res2(g,g+3);
	//buchberger(res);

	cout<<"groebner of\n"<<res2<<endl;
	buchberger(res2);
	cout<<"result:"<<res2<<endl;

}

void testgroebner()
{
	polyn::expandPowers=false;
	polyn s(1.0,1),t(1.0,2),x(1.0,3),y(1.0,4),z(1.0,5);
	
	
	
	
	polyn base[]={x-s*t,y-s*t*t,z-s*s};
	
	vector<polyn> vx(base,base+3);

	vector<polyn> result=vx;
	buchberger(result);
	cout<<'\n'<<result<<endl;
	
	typedef poly<20,infbruch> tpl;

	cout<<"2nd try:"<<endl;
	vector<tpl > vx2(base,base+3);
	
	buchberger(vx2);
	cout<<'\n'<<vx2<<endl;


	polyn base2[]={x*x+y*y+z*z-1,x*x-y+z*z,x-z};
	vector<tpl > vx3(base2,base2+3);

	buchberger(vx3);
	cout<<'\n'<<vx3<<endl;

	
	test2();
	
	char c;
	cin>>c;
}


};

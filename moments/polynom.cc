/*
	operations on polynomials

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


#include<algorithm>
#include<cassert>
#include<string>

#include "tensWithSymm.h"

#include "polynom.h"

#include "tensorGraph.h"

using namespace std;

int polynom::ioFormat=polynom::RAW;

istream&operator>>(istream&in,polynom&p)
{
	if(p.ioFormat==p.RAW)
	{
		p.readRaw(in);
	}
	else{
		assert(false);
	}
	return in;
}


metainfo::~metainfo(){}

vector<metainfo*> * metainfo::types=0;


void polynom::printRaw(ostream&o) const
{
	o<<" ( ";
  for(int i=0;i<(int)p.size();i++){
		if(i>0) o<<" + ";
		vector<int>::const_iterator cit
			= p[i].begin(),citend=p[i].end();
		o<<*cit; 
		for(++cit;cit!=citend;++cit)
		{
			o<<" * "<< *cit;
		}
	}
	o<<" )";

	if(mi!=0)
	{
		o<<'m'<<' ';
		mi->save(o);
	}
	else
		o<<'e';
}

void polynom::readRaw(istream&in)
{
	vector<int> summand;
	char c; int i;
	in>>c;
	if( c!='(' ) return;
	do {
		summand.clear();
		do{ in>>i>>c; summand.push_back(i); } while(c=='*');
		p.push_back(summand);
	}	while(c!=')');
	in>>c;
	if(c=='m')
	{
		mi=metainfo::load(in);
	}
}

ostream&operator<<(ostream&o,const polynom&p)
{
	if(p.ioFormat==p.RAW)
	{
		p.printRaw(o);
		return o;
	}
	if(p.ioFormat==p.PRINTFORC)
	{p.printForC(o);return o;}

  for(int i=0;i<(int)p.p.size();i++){
    if(i!=0)o<<" + " ;

    int k=0;
    while(i+1<(int)p.p.size() && p.p[i]==p.p[i+1] ) k++,i++;
      if(k>0)
        o<<k+1<<'*';



    const vi & ppi =p.p[i];
    for(int j=0;j<(int)ppi.size();j++)
    {
      printindex(o,ppi[j]);
      int k=0;
      while( j+1<(int)ppi.size() && ppi[j]==ppi[j+1] ) k++,j++;
      if(k>0)
        o<<'^'<<k+1;

      if(j<(int)ppi.size()-1)
        o<<'*';
      
    }


  }

  return o;
}

void polynom::printForC(ostream&o) const
{
	
  for(int i=0;i<(int)p.size();i++){
    if(i!=0)o<<" + " ;

    int k=0;
    while(i+1<(int)p.size() && p[i]==p[i+1] ) k++,i++;
    
		if(k>0)
      o<<k+1<<'*';

    const vi & ppi =p[i];
    for(int j=0;j<(int)ppi.size();j++)
    {
      if(j>0)o<<'*';
			o<<"A["<<ppi[j]<<"]";      
    }


  }
}


int encodeIndices(const int*inds,int ord, int dim, int id)
{
  int ret=0;
  for(int i=ord-1;i>=0;--i)
    ret*=dim,ret+=inds[i];

  ret*=10;ret+=id;
  ret*=10;ret+=dim;
  ret*=10;ret+=ord;
  return ret;
}


void printindex(ostream&o,int i)
{
  int ord = i%10;
  i/=10;
  int dim = i%10;
  i/=10;
  o<<char('a'+(i%10));
  i/=10;

  int siz=1;
  for(int j=0;j<ord-1;++j)
    siz*=dim;
  for(;siz;siz/=dim){
    o<< (i/siz)%dim;
  }
}



//--------------------------------------------------------------------


polyn::polyn(const polynom &p)
{
	//buffer var
	summand sm;

	for(int i=0;i<p.p.size();++i)
	{
		sm.f.clear();
		//count coefficient.
		sm.c=1;
    while(i+1<(int)p.p.size() && p.p[i]==p.p[i+1] ) 
			++sm.c,++i;

		factor fk;
		const vector<int> &ppi=p.p[i];
		for(int j=0;j<ppi.size();++j)
		{
			//count power
			fk.p=1;
		  while(j+1<(int)ppi.size() && ppi[j+1]==ppi[j] ) 
				++fk.p,++j;
			fk.i=ppi[j];
			sm.f.push_back(fk);		
		}
		
		s.push_back(sm);
	}

	s.sort();
}

bool polyn::expandPowers=true;


ostream& operator <<(ostream&o,const polyn::summand&s)
{
	//should we expand powers (i.e x^3= x*x*x ?)
	if(s.f.empty())
		o<<s.c;
	else 
	{ 
		if(s.c==-1)o<<'-';
		else if(s.c!=1)o<< s.c <<'*'; 
	}
	for(int i=0;i<s.f.size();++i)
	{
		if(i>0)o<<'*';
		o<<"A["<<s.f[i].i<<"]";
		if(s.f[i].p>1)
		{
			if(polyn::expandPowers)
				for(int j=1;j<s.f[i].p;++j)
					o<<"*A["<<s.f[i].i<<"]";
			else
				o<<"^"<<s.f[i].p;			
		}
	}

	return o;
}

ostream&operator<<(ostream&o,const polyn&p)
{
	list<polyn::summand>::const_iterator sit;

	o<<" ( ";
	for(sit=p.s.begin();sit!=p.s.end();++sit)
	{
		if(sit!=p.s.begin())
			o<<" + ";
		o<<*sit;
	}
	o<<" ) ";
	return o;
}




void polyn::summand::printWithParamNames(
	ostream&o,
	const std::vector<std::string>&names,
	bool expandPowers,
	bool withcoeff
) const
{

	if(withcoeff)
	{
		if(f.empty())
			o<<c;
		else 
		{ 
			if(c==-1)o<<'-';
			else if(c!=1)
				o<< c <<'*'; 
		}
	}

	for(int i=0;i<f.size();++i)
	{
		if(i>0)o<<'*';
		o << names[f[i].i];
		if(f[i].p>1)
		{
			if(expandPowers)
				for(int j=1;j<f[i].p;++j)
					o<<'*'<<names[f[i].i];
			else
				o<<"^"<<f[i].p;			
		}
	}
}


bool cmpsumptr(const polyn::summand* a,
							 const polyn::summand* b)
{
	return a->c < b->c;
}

/*
print polynomial 
with common coefficients factored out
and parameters replaced with paramnames;
used for repeatTensor output.
*/
void polyn::printOptCommonCoeff
( ostream&o,
  const vector<string> &paramnames ) const
{
	//sort by factors:
	vector< const polyn::summand* > sums;
	list<polyn::summand>::const_iterator sit;
	for(sit=s.begin();sit!=s.end();++sit)
	{
		sums.push_back(&*sit);
	}

	sort(sums.begin(),sums.end(),cmpsumptr);

	for(int i=0;i<sums.size();++i)
	{
		double c= sums[i]->c;
		
		if(i>0)
			o<<" + ";

		if( c!=1 && c !=-1)
		{		
			//group equal factors
			int j=i+1;
			while(j<sums.size() && c == sums[j]->c)++j;				
			--j;

			//if more than one in a group:
			if( j > i )
			{
				o << c << " * ( ";

				sums[i]->printWithParamNames(o,paramnames,true,false); 		
				do
				{
					++i;
					o<<" + ";
					sums[i]->printWithParamNames(o,paramnames,true,false); 		
				}while(i<j);
			
				o<<" )";
			}
			else
				sums[i]->printWithParamNames(o,paramnames,true,true); 		
		}			
		else
		{
			sums[i]->printWithParamNames(o,paramnames,true,true); 			
		}

	}

}


//--------------------------------------------------------------------


void tensmetainfo::save_(ostream&out) const
	{
		out<<sourceTensors.size()<<'\n';
		for(unsigned i=0;i<sourceTensors.size();++i)
			out<<sourceTensors[i]<<' ';
		out<<'\n'<<pairs.size()<<'\n';

		for(unsigned i=0;i<pairs.size();++i)		
			out<<pairs[i].a<<' '<<pairs[i].b<<"  ";
		out<<'\n';		
	}

	void tensmetainfo::load_(istream&in)
	{
		size_t siz;
		in>>siz;
		sourceTensors.resize(siz);
		for(unsigned i=0;i<sourceTensors.size();++i)
			in>>sourceTensors[i];
		in>>siz;
		pairs.resize(siz);
		for(unsigned i=0;i<pairs.size();++i)		
			in>>pairs[i].a>>pairs[i].b;	
	}

	const char*tensmetainfo::name() const {return "tensmetainfo";}

	void tensmetainfo::print(ostream&o) const
	{
		for(int i=0;i<sourceTensors.size();++i)
			o<<sourceTensors[i]<<"A ";
		o<<':'<<pairs;
	}

	void tensmetainfo::init	( const tensSymmInfo&pi ,const vpr & pairs_)
	{
		grp=pi.groups();
		tgrp=pi.tgroups();

			pairs=pairs_;
			sourceTensors.resize(pi.tgroups().size(),0);
			int gid=0;
			for(int i=0;i<pi.tgroups().size();++i)
			{
				for(int j=0;j<pi.tgroups()[i];++j,++gid)
				{					
					sourceTensors[i]+=pi.groups()[gid];
				}
			}
	}


tmireg<tensmetainfo> tmireg_instance;

metainfo* tensmetainfo::clone() const{ 
	return new tensmetainfo(*this); 
}

tensmetainfo::~tensmetainfo(){};

//print it out in the "Einstein sum convention" format
void tensmetainfo::printEsumForTex(ostream&o) const
{
	//for every tensor index, the pair id
	vector< int > p;

	//init pairpertens array sizes
	size_t totord=0;
	for(size_t i=0;i<sourceTensors.size();++i)
	{
		totord+=sourceTensors[i];
	}

	p.resize(totord);

	//enter pair ids at every tensor index
	for(size_t i=0;i<pairs.size();++i)
	{
		p[pairs[i].a]=i;
		p[pairs[i].b]=i;
	}

	size_t pi=0;
	for(size_t i=0;i<sourceTensors.size();++i)
	{
		//o<<"\\A("<<sourceTensors[i]<<")";
		o<<'A';
		for(int j=0;j<sourceTensors[i];++j,++pi)
		{
			int pairid=p[pi];
			o<<"{}";
			if(pairs[pairid].a == pi)
				o<<'_';
			else
				o<<'^';
			o<< char('i'+pairid);				
		}
	}
}

//print it out in the "contraction of tensor product" format
void tensmetainfo::printPairsProdForTex(ostream&o) const
{

	o<<"\\sum_{";
	for(size_t i=0;i<pairs.size();++i)
	{
		o<<'('<<pairs[i].a<<','<<pairs[i].b<<')';
	}
	o<<'}';

	size_t pi=0;
	for(size_t i=0;i<sourceTensors.size();++i)
	{
		if(i>0)
			o<<"\\otimes";
		o<<"\\A("<<sourceTensors[i]-1<<")";
	}
}


//print it out as Graph in postscript format
void tensmetainfo::printGraphPS(ostream&o) const
{
	basicTensSymmInfo bts;
	bts.grp=grp;bts.tgrp=tgrp;
	tensorGraph tg(bts,pairs);
	tg.normalize();
	

	int maxid=0;
	for(int i=0;i<sourceTensors.size();++i)
	{
		if(maxid<sourceTensors[i])maxid=sourceTensors[i];
	}
	vector<vector<int> > gspt(maxid+1);

	int gid=0;
	for(int i=0;i<tgrp.size();++i)
	{
		vector<int> &g=gspt[sourceTensors[i]];
		g.clear();
		for(int j=0;j<tgrp[i];++j,++gid)
		{
			g.push_back(grp[gid]);
		}
	}

	tg.printPS(o,gspt,20,true);
}


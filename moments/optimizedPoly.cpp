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


#include "optimizedPoly.h"
#include<algorithm>

int optimizedPoly::vgl(optimizedPoly&b)
{
		children.sort();
		b.children.sort();

		list<pair<int,optimizedPoly> >::iterator i,j;
		for(i=children.begin(),j=b.children.begin();
				i!=children.end() && j!=b.children.end();
				++i,++j)
		{
			int v=vgl(*i,*j);
			if(v!=0)return v;
		}

		return (i==children.end())-(j==b.children.end());				
}


void optimizedPoly::init(list<polyn::summand>& s)
	{
		//cout<<'{'<<flush;
		int maxparam=0;
		list<polyn::summand>::iterator sit;
		//cout<<'('<<flush;
		for(sit=s.begin();sit!=s.end();++sit)
				for(int i=0;i<sit->f.size();++i)
				{ 
					int m = sit->f[i].i;
				  //cout<<m<<' '<<flush;
				  if(m>maxparam)maxparam=m;
				}
		//cout<<')'<<flush;
		
		//for every param the list of summands		
		map<int,list<polyn::summand> > groups;

		{ //livetime block of spp
			
			//summands per parameter
			vector<int> spp(maxparam+1,0);

			//cout<<'('<<flush;
			//count summand numbers for each parameter
			for(sit=s.begin();sit!=s.end();++sit)
			{						
				for(int j=0;j<sit->f.size();++j)
				{
					//cout<<'.'<<flush;
					spp[sit->f[j].i]++;			
				}		
			}
			//cout<<')'<<flush;

			//group summands by their parameters
			for(sit=s.begin();sit!=s.end();)
			{						
				int maxp=0;//spp of idmax
				int idmax=-1;//id of the param which has the biggest spp in this summand
				for(int j=0;j<sit->f.size();++j)
				{
					int nid=sit->f[j].i;
					int nmaxp=spp[nid];
					if(nmaxp> maxp)
					{ idmax = nid; maxp=nmaxp; }
				}

				if(idmax>=0){
					//put each summand into the list it belongs to
						list<polyn::summand> &l=groups[idmax];
					
						//cout<<'<'<<flush;
					list<polyn::summand>::iterator sit2=sit;
					++sit;

					l.splice(l.begin(),s,sit2);
						//cout<<'>'<<flush;
				}
				else
					++sit;
			}
		}

		//process the summands which are only a number:
		//(and for that are left in s)
		for( sit = s.begin();sit!=s.end();++sit)
		{
			children.push_back(pair<int,optimizedPoly>(
				sit->c,optimizedPoly()));			
		}

		//cout<<'%'<<flush;

		//for every group, do the necessary
		map<int,list<polyn::summand> >::iterator it;
		for(it=groups.begin();it!=groups.end();++it)
		{
			int pid=it->first;
			list<polyn::summand> &l = it->second;
			list<polyn::summand>::iterator lit;

			if(l.empty()) continue;

			//cout<<'&'<<endl;
			//"divide" all factors by pid
			for(lit=l.begin();lit!=l.end();++lit)
			{
				vector<polyn::factor> & f = lit->f;
				assert(!f.empty());
				//if(f.empty()){cout<<'_'<<flush;continue;}

				int j=0;
				polyn::factor cmp; cmp.i = pid;
				vector<polyn::factor>::iterator l=lower_bound(f.begin(),f.end(),cmp,polyn::lti());

				/*
				if(l==f.end())
				{
					cout<<'!'<<flush;
					for(int i=0;i<f.size();++i)
						cout<<f[i].i<<'^'<<f[i].p<<'*'<<flush;
				}
				*/
				assert(l!=f.end());

				--l->p;
				if(l->p==0)
					f.erase(l);
			}

			children.push_back(pair<int,optimizedPoly>(pid,optimizedPoly()));
			
			children.back().second.init(l);								
		}

		//cout<<'}'<<flush;

	}


	
	
	//produce sth C++ can eat
ostream& print (ostream&o, const optimizedPoly & p, int maxid=-1)
{
	if(p.children.empty())
	{o<<"(0)" ;return o;}
	


	//only one subpoly
	bool onlyOne= ++p.children.begin() ==p.children.end();

	if(!onlyOne)o<<'(';
	list<pair<int,optimizedPoly> >::const_iterator lit;
	for(lit=p.children.begin();lit!=p.children.end();++lit)
	{
		if(lit!=p.children.begin())
			o<<'+';

		const list<pair<int,optimizedPoly> > & lp = lit->second.children;

		if( lp.empty() )
		{
			o<<lit->first;
		}
		else{
		
			bool coeffOne=
				++lp.begin() ==lp.end() 
				&&
				lp.front().second.children.empty() 
				&&		
				lp.front().first==1
				;

			//don't draw coefficient that is one
			if(!coeffOne)print(o,lit->second,maxid)<<"*";

			if(maxid>=0 && lit->first >maxid)
				o<<"B["<<lit->first-maxid-1<<"]";
			else
				o<<"A["<<lit->first<<"]";
		}
	}
	if(!onlyOne)o<<')';

	return o;
}


	

ostream& operator << (ostream&o, const optimizedPoly & p)
{
	return print(o,p);
}
	
//---------------------------------------------------------
	
	//returns the maximum order
	int optimizedPolyCommonSubBase::recinit(optimizedPoly & p)
	{
		optimizedPoly::summands::iterator i;

		if(p.children.empty())
			return 0;
		
		int reti=0;

		for(i=p.children.begin();i!=p.children.end();++i)
		{
			//if i->first is a parameter and not a number:
			if( ! i->second.children.empty() )
			{ 
				if(maxid<i->first)maxid=i->first; 
				
				//if i->first is only multipied with a number, and it is 1:
				int o=recinit( i->second );				
				if(o>=prodOrders.size())prodOrders.resize(o+1);

				//don't insert products with a 1: 
				if( ! (o==1 && i->second.children.back().first==1))
					prodOrders[o].push_back(i);
	
				if(o>reti)reti=o;
			}
		}
		if(reti>=sumOrders.size())sumOrders.resize(reti+1);
		if(++p.children.begin()!=p.children.end())
			sumOrders[reti].push_back( &p );		

		return reti+1;
	}
	
//............................................................

	void optimizedPolyCommonSubBase::initBase(int maxorder)
	{
		

		for(int o=0;o<maxorder;++o)
		{

			if(prodOrders.size()>o)
			{
				prods &p=prodOrders[o];
				sort(p.begin(),p.end(),cmpProd());
				prods::iterator i=p.begin(),oi=i;
				while(i!=p.end())
				{
					optimizedPoly::summand &soi= **oi;
					++i;
					while(i!=p.end() &&  soi == **i )
						++i;

					//if at least two are equal:
					if(oi+1!=i)
					{
						cout<<"prods order"<<o<<":A["<<soi.first<<"]*"<<soi.second<<endl;

						np.push_back(optimizedPoly());
						np.back().children.push_back(soi);

						//replace **oi by a reference to its id in np
						soi.first =maxid+np.size();
						soi.second.children.clear();
						soi.second.children.push_back(
						pair<int,optimizedPoly>(1,optimizedPoly()));
						//copy single entry to others
						for(++oi;oi!=i;++oi)
							**oi = soi;
					}
					oi=i;
				}
			}

			if(sumOrders.size()>o)
			{
				sums & s =sumOrders[o];
				sort(s.begin(),s.end(),cmpSum());
				sums::iterator i=s.begin(),oi=i;
				while(i!=s.end())
				{
					optimizedPoly &optp=**oi;
					++i;
					while(i!=s.end() &&  **oi == **i )
						++i;

					if(oi+1!=i){
						cout<<"sums order"<<o<<":"<<optp<<endl;
						np.push_back(optp);

						//replace **oi by a reference to its id in np
						optp.children.clear();
						optp.children.push_back(optimizedPoly::summand());
						optimizedPoly::summand &soi= optp.children.back();
						soi.first =maxid+np.size();
						soi.second.children.clear();
						soi.second.children.push_back(
							pair<int,optimizedPoly>(1,optimizedPoly()));

						for(++oi;oi!=i;++oi)
							**oi = optp;
					}			
					oi=i;
				}

			}			

		}
	}

//---------------------------------------------------------

ostream&operator<<(ostream&o,optimizedPolyCommonSub &s)
{	
	o<<"//T B["<<s.np.size()<<"];\n";
	for(int ord=0;ord <= s.np.size();++ord)
	{
		optimizedPoly&p = ord<s.np.size()? s.np[ord] : s.op;
		if( ord<s.np.size())
			o<<"B["<<ord<<"]=";
		else
			o<<" return ";

		if(p.children.empty())
		{o<<"0;\n";continue;}

		print(o,p,s.maxid)<<';'<<'\n';
	}

	return o;
}




	void optimizedPolyCommonSub::init(const optimizedPoly& p)
	{
		op=p;
		prodOrders.reserve(100);
		sumOrders.reserve(100);
		np.reserve(1000);
		maxid=0;//(maxid will be set in recinit)
		
		int maxorder = recinit(op); 

		initBase(maxorder);
	}


	
//---------------------------------------------------------
//methods/functions for optimizedPolyCommonSubMulti


ostream & optimizedPolyCommonSubMulti::print
(ostream&o, const vector<metainfo*>*mi) const
{	
	if(mi)
		assert(mi->size()==op.size());
	o<<"T B["<<np.size()<<"];\n";
	for(int i=0;i <np.size();++i)
	{
		o<<"B["<<i<<"]=";
		::print(o,np[i],maxid)<<';'<<'\n';
	}
	
	//the return values
	for(int i=0;i<op.size();++i)
	{
		if(mi && (*mi)[i] )
			o<<"\n//"<<*(*mi)[i]<<'\n'<<'\n';
		o<<"M["<<i<<"]=";
		::print(o,op[i],maxid);
		o<<';'<<'\n';
	}
	
	return o;
}

ostream&operator<<(ostream&o,optimizedPolyCommonSubMulti&s)
{	return s.print(o);}




	void optimizedPolyCommonSubMulti::init(const vector<optimizedPoly>& pp)
	{
		//set the original polynomes to np
		op = pp;

		maxid=0;//(maxid will be set in recinit)
		int maxorder=0;
		for(int i=0;i<pp.size();++i)
		{
			int o = recinit(op[i]);
			if(o>maxorder)maxorder=o;
		}

		initBase(maxorder);
	}


	
	
	
	
	
	
	

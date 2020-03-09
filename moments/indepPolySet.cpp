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


#include "indepPolySet.h"
#include <float.h>

//fill p with prime numbers;numthere:number thats  are already there
void getPrimes(vector<int>&p,int numthere)
{
  assert(numthere<=(int)p.size());
  vector<int>::iterator current=p.begin(),j;
  int num=2;
  if(numthere){
    num=p[numthere-1]+1;
    current+=numthere;
  }

  for(;current!=p.end()&& num>0 ;++num)
  {
    bool isprime=true;
    for(j=p.begin(); j!=current && *j * *j <= num; ++j)
      if( num % *j == 0 ) 
        { isprime=false;break;}

    if(isprime)
      *(current++)=num;
  }
}

void indepPolySet::setParamValues(const std::vector<mynum>& paramValues_)
	{
		usePresetParamValues=true;paramValues=paramValues_;	

		if(polys.size()!=0)
		{
			cout<<"was not empty, recalculating polynomes and matrix"<<endl;
			mat.clear();

			tpolyset polys2;
			polys2.swap(polys);
			
			tpolyset::iterator it;
			for(it=polys2.begin();it!=polys2.end();++it)
			{
				addpolynom(*it);	
			}
		}

	}

	void indepPolySet::printMat(std::ostream&out)
	{
		for(int i=0;i<mat.size();++i){
			for(int j=0;j<mat[i].size();++j)
			{
				out<<mat[i][j]<<' ';
			}
			out<<endl;
		}
	}
 bool indepPolySet::addpolynom(const polynom&x)
  {
		/* test if we hve the same already there*/

		if(polys.find(x)!=polys.end())
		{
			//cout<<"\nl"<<endl;
			return false;
		}
		//cout<<"not yet present.adding it."<<endl;

    /** add parameters in x to parameter list*/

		
    //set rounding model to round_down so that
    //interval arithmetics in myinterval uses outward rounding
#ifndef GNUCC
		unsigned int ocfp = _controlfp(0,0);
		_controlfp(_RC_DOWN,_MCW_RC);
#endif		
    if(x.getMaxParamID() >= (int)paramThere.size()){

      paramThere.resize(x.getMaxParamID()+1);

			if(!usePresetParamValues){
				//add representative values for the new params
				int numthere=(int)primes.size();
				paramValues.resize(paramThere.size());
				primes.resize(paramThere.size());
				getPrimes(primes,numthere);
				for(int i=numthere;i<(int)primes.size();++i)
        paramValues[i]=mynum(infnum(primes[i])); //sqrt((double)primes[i]);
			}
			else{
				//if we don't have enough preset values, fill with random numbers
				if(paramValues.size()<paramThere.size()){
					int i=paramValues.size();
					paramValues.resize(paramThere.size()*2);
					for(;i<paramValues.size();++i)
						paramValues[i]=(mynum)rand();			
				}
			}

    }
  
		//cout<<"fillin param list"<<endl;
    vvi::const_iterator j;
    vi::const_iterator k;
    for(j=x.p.begin();j!=x.p.end();++j)
      for(k=j->begin();k!=j->end();++k)
        if(!paramThere[*k])
            paramIds.push_back(*k), 
            paramThere[*k]=true;

    //adapt matrix size   
    vector<vector<mynum > >::iterator mit;
    for(mit=mat.begin();mit!=mat.end();++mit)
      mit->resize(paramIds.size());

    //--------------------------------------------------
    //cout<<"evaluating derivatives of"<<x
		//	<<"with "<<paramIds.size()<<"parameters"<<endl;
    //prepare candidate row for the matrix
    vector<mynum> cand(paramIds.size());
    for(int i=0;i<(int)paramIds.size();++i)
    {cand[i] = mynum(x.evalderiv<mynum>(paramValues,paramIds[i]));    
//     cout<<"d/d"<<paramIds[i]<<" "<<cand[i]<<endl; 
//		 cout<<" i "<<i<<" paramIds.size() "<<paramIds.size()<<endl;
    }

    //--------------------------------------------------
    //apply upper rows on new row    
    for(int i=0;i<(int)mat.size();++i){
      vector<mynum> &m=mat[i];
      mynum c=cand[i],d;
      for(int j=i; j< (int)m.size();++j){
        d=c; d*= m[j];
        cand[j] -= d;      
      }
    }

//cout<<" before.."<<flush;
    //find element not null
    int maxi=0;
    while(maxi<cand.size() && !cand[maxi] )
      ++maxi;
	//	cout<<" after:maxi="<<maxi<<endl;

    //if maximum value too low, new polynom seems to be dependant
    if(maxi==cand.size())
    {
      /*
        cout.precision(16);cout<<"maxval:"<<maxval<<endl;
        for(int i=0;i<matsiz;++i)
        {
          vector<mynum> &m =mat[i];
          for(int j=0;j<m.size();++j)
            cout<<m[j]<<' ';

          cout<<endl;
        }
        */
      //restore float control word
#ifndef GNUCC
        _controlfp(ocfp,_MCW_RC);
#endif
        return false;
    }
    else{
     
			/*
        cout<<endl;
        for(int j=0;j<cand.size();++j)
          cout <<(!cand[j] ? " ": " +" )<< cand[j];
        cout<<endl;
      */
      
    }


    int matsiz=(int)mat.size();
    //normalize row(so that max. norm is 1)
    mynum maxval=cand[maxi];
    for(int i=matsiz;i<(int)cand.size();++i)
		{
			//cout<<i<<endl;
      cand[i]/=maxval;
		}

    //polynome is independent, so add it(plus it's description)
		//cout<<"insert.."<<flush;
		polys.insert(x);
    mat.resize(matsiz+1);
    mat.back().swap(cand);
		//cout<<"ready. interchange..."<<flush;

    
    //interchange param positions:
    //(so that on the diagonal is not a null)
    if(maxi!=matsiz){
      int x=paramIds[maxi];paramIds[maxi]=paramIds[matsiz]; paramIds[matsiz]=x;
      mynum xch;
      for(int i=0;i<(int)mat.size();++i){
        vector<mynum> & mi=mat[i];
        xch=mi[maxi];mi[maxi]=mi[matsiz];mi[matsiz]=xch;
      }
    }
		//cout<<"ready."<<endl;
#ifndef GNUCC
    _controlfp(ocfp,_MCW_RC);
#endif
    return true;
  }


  int indepPolySet::join(const indepPolySet&xx)
  {
    int ret=0;
		tpolyset::const_iterator it=xx.polys.begin();
    for(;it!=xx.polys.end();++it)
    {
      ret+=addpolynom(*it);
    
    }
    return ret;  
  }

	void indepPolySet::save(ostream&o,const tpolyset&a)
	{
		o<<a.size()<<'\n';
		for(tpolyset::const_iterator it=a.begin();it!=a.end();++it)
		{
			it->printRaw(o);
		}	
	}

	void indepPolySet::load(istream&in,tpolyset&a)
	{
		size_t siz;
		in>>siz;
		
		for(size_t i=0;i<siz;++i)
		{
			polynom p;
			p.readRaw(in);
			a.insert(p);
		}	
	}

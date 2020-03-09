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

#include "optimizedFold.h"

int repeatTensor::nextid=0;

repeatTensor::~repeatTensor()
{}

void repeatTensor::createTindexPerComp(vector<int>&idxpert) const
{	
	if(comps.empty())
		idxpert.resize(numcomps);
	else
		idxpert.resize(comps.size());

	for(int i=0;i<t.t.size();++i)
			idxpert[t.t[i]] = i;
}

string& repeatTensor::writeTindex(string&ret,int ii) const
 {
	  ret="";
		for(int j=0;j<ord();++j)
		{
			ostringstream os;
			os<<' '<<(ii % dim());
			ret=os.str()+ret;
			ii /=dim();
		}	
		return ret;
	}


polyn repeatTensor::operator[](const int* inds) const
{
	return get(inds);
}

void repeatTensor::init(const tens<polyn>&tt)
{
	  id=nextid;
	  ++nextid;
	
		tvisitor tv(tt,t);
		
		visitTensorIndices<tvisitor> vt(tt,tv);
		
		comps.resize(tv.nextid);
		
		for(map<polyn,int>::iterator it= tv.s.begin();
		    it!=tv.s.end();
		   ++it)
				comps[it->second] = it->first;   				
	
		numcomps=comps.size();

}


repeatTensor::repeatTensor(const tens<polyn>&tt,const string& nam)
		:tens<polyn>(tt),t(tt.ord(),tt.dim()),comment(nam)
	{			
		numcomps=-1;
		init(tt);
	}





repeatTensor::repeatTensor(const tens<polynom>&tt,const string& nam)
		:tens<polyn>(tt),t(tt.ord(),tt.dim()),comment(nam)
	{			
		numcomps=-1;
		struct polynom2poly:public tens<polyn>
		{
			const tens<polynom>&t;
			polynom2poly(const tens<polynom>&tt)
			:tens<polyn>(tt),t(tt)
			{}

			polyn operator[](const int* inds) const
			{
				return t[inds];
			}
		};

		polynom2poly ppt(tt);

		init(ppt);
	}
                                                                       
//print independent components as polynomes in 
//the source tensor components
void repeatTensor::printForC(ostream&o,const char*varname) const
{
	  string thisname=id2name(id);	  
	  o<<
	  "\n"
	  "//"<<comment<<"\n"
		"T "<<thisname<<"["<<comps.size()<<"];\n"
	  ;

		vector<int> tindexpercomp;		
		createTindexPerComp(tindexpercomp);
	  string buf;

		int compsinline=0;
		for(int i=0;i<comps.size();++i){
			o<<"/*"<<writeTindex(buf,tindexpercomp[i])<<"*/ ";
			o<<thisname<<"["<<i<<"]=";
			o<<comps[i]<<";\n";
		}						
}


//--------------------------------------------------------------
polyn repeatTensorAdd::operator[](const int* inds) const
{
	return polyn(1,t.t.get(inds)+addme);
}

//----------------------------------------------------------------
          
void repeatTensorUneval::init(const tens<polyn>&a , const vector<polyn>&bel)
  
  { 

		//this case is used if not enough mem is there for full optimization
		if(bel.empty())
		{
			cout<<"The bad case!  "<<endl;
			repeatTensor rt(a);
			compsUneval.swap(rt.comps);
			numcomps = compsUneval.size();
			t.t.swap(rt.t.t);
			return;
		}

		creationVisitor cv(a, t.t.begin() , bel );
		
		visitTensorIndices<creationVisitor> vti(t,cv);		

		creationVisitor::tevalps::iterator it;

		//the components to evaluated polynomes
		comps.resize(cv.id);
		compsUneval.resize(cv.id);
		
		for(it=cv.evalps.begin();it!=cv.evalps.end();++it)
		{
			int idx=it->second->second;
			comps[idx] = it->first;
			compsUneval[idx]=it->second->first;
		}			
	}


repeatTensorUneval::~repeatTensorUneval(void)
{
}

//-------------------------------------------------------------------


void RTpartialFold::printForC(ostream&o,const char*varname) const
{
	  string aname=id2name(aid);
	  string thisname=id2name(id);
	  
	  o<<
	  "\n"
	  "/*  sum_{"<<pairs<<"} ( "<<aname<<" )  */\n"
	  "/* ( = "<<comment<<" )  */\n";
		if(!varname)
			o<<"T "<<thisname<<"["<<comps.size()<<"];\n";
		
	  vector<string> paramNames(acomps);
	  
	  for(int i=0;i<acomps;++i)
	  {
			ostringstream oss;
			oss<<aname<<"["<<i<<"]";
			paramNames[i]=oss.str();
	  }
	  
		vector<int> tindexpercomp;		
		string buf;
		createTindexPerComp(tindexpercomp);

	  for(int i=0;i<compsUneval.size();++i)
	  {
			o<<"/*"<<writeTindex(buf,tindexpercomp[i])<<"*/ ";
			if(compsUneval.size()==1 &&varname)
				o<<varname<<"=";
			else
	  		o<<thisname<<"["<<i<<"]=";

			compsUneval[i].printOptCommonCoeff(o,paramNames);

	  	o<<";\n";
	  }
	   
	
}


RTpartialFold::	RTpartialFold(const repeatTensor * a,
															const vi&perm_,const vpr&inds,
															 const string&nam)
		:
		repeatTensorUneval(a->ord()-inds.size()*2,a->dim(),nam),
		pairs(inds),perm(perm_)
	{
		aid = a->id;
		acomps = a->comps.size();		
		repeatTensorAdd aa(*a,0);
		indexPerm<polyn> paa(aa,perm);
		partialFold<polyn> pf(aa,inds);
		init(pf,a->comps);	
		if(comment=="")
		{
			ostringstream s;
			s<<"sum_{"<<inds<<"}("<<a->comment<<")";
			s.str(comment); 	
		}
	}

RTpartialFold::~RTpartialFold(){}

//------------------------------------------------------------


void RTpartialFoldProd::printForC(ostream&o, const char*varname) const
{
	  string aname=id2name(aid),bname=id2name(bid);
	  string thisname=id2name(id);
	  
	  o<<
	  "\n"
	  "/*  sum_{"<<pairs<<"} "<<aname<<" "<<bname<<"  */\n"
	  "/*  (="<<comment<<")  */\n";
		if(!varname)
			o<<"T "<<thisname<<"["<<comps.size()<<"];\n"
	  ;
		
	  vector<string> paramNames(acomps+bcomps);
	  vector<string>::iterator pniter=paramNames.begin();
	  
	  for(int i=0;i<acomps;++i,++pniter)
	  {
			ostringstream oss;
	  	oss <<aname<<"["<<i<<"]";
			*pniter=oss.str();
	  }
	  for(int i=0;i<bcomps;++i,++pniter)
	  {
			ostringstream oss;
	    oss <<bname<<"["<<i<<"]";
			*pniter=oss.str();
	  }
	  
		vector<int> tindexpercomp;		
		createTindexPerComp(tindexpercomp);
	  string buf;
	  for(int i=0;i<compsUneval.size();++i)
	  {
			o<<"/*"<<writeTindex(buf,tindexpercomp[i])<<"*/ ";
			
			if(compsUneval.size()==1 &&varname)
				o<<varname<<"=";
			else
	  		o<<thisname<<"["<<i<<"]=";

			compsUneval[i].printOptCommonCoeff(o,paramNames);
	  	o<<";\n";
	  }
	   
	

}


double RTpartialFoldProd::estimateEvalPolyMem(
		const repeatTensor* a,
		const repeatTensor* b)
{

		//try to estimate size of resulting polynomes in comps:
		//try to get the maximum size of one polynome,
		//then multiply by dim^(ord+indpairs.size())
	  // dim^ord for the number of entries entry times a size of 
	  //dim^indpairs.size() for a single entry
		
		if(a->comps.empty()||b->comps.empty())
			return 0;

		size_t maxa=0;
		for(int i=0;i<a->comps.size();++i)
		{
			size_t s=a->comps[i].s.size();
			if(maxa < s)	maxa = s;
		}
		
		size_t maxb=0;
		for(int i=0;i<b->comps.size();++i)
		{
			int s=b->comps[i].s.size();
			if(maxb < s)	maxb = s;
		}

		size_t fnuma = a->comps[0].s.front().f.size(); 
		size_t fnumb = b->comps[0].s.front().f.size(); 
		size_t nssiz = 
			(fnuma+fnumb)*sizeof(polyn::factor)
			+sizeof(polyn::summand)
			+16;
		
		double ret=double(maxa) * maxb * nssiz;
		
		for(int i=0;i<ord()+pairs.size();++i)
			ret*=dim();
		
		//ret *= a->comps.size() * b->comps.size();

		return ret;
}

	RTpartialFoldProd::RTpartialFoldProd
		(const repeatTensor* a,
		 const vi&perma,
		 const repeatTensor* b,
		 const vi&permb,
			const vpr &indpairs)
			:	repeatTensorUneval(
					int(a->ord()+b->ord()-indpairs.size()*2),a->dim(),""
				),
				pairs(indpairs),
				perma(perma),permb(permb)
	{		
		aid=a->id;bid=b->id;
		
		acomps=a->comps.size();
		
		bcomps=b->comps.size();

		repeatTensorAdd	aa(*a,0),bb(*b,a->comps.size());
		
		indexPerm<polyn> paa(aa,perma),pbb(bb,permb);

		assert(a->dim()==b->dim());
		tensProd<polyn> aabb(paa,pbb);		
		
		partialFold<polyn> fold(aabb,pairs);

		vector<polyn> bel;
		double eepm=estimateEvalPolyMem(a,b);
		if(eepm>=1e9)
		{
			//if too much mem, don't evaluate polynomes
			cerr<<"Not using polynomes in original components:\n"
				"would need too much memory: ("<<eepm*1e-9<<"GB)"<<endl;
		
		}
		else if(a->comps.empty() || b->comps.empty() )
		{
			//if too much mem, don't evaluate polynomes
			cerr<<"Not using polynomes in original components:\n"
				"because they are not given "<<endl;		
		}
		else
		{
			bel.resize(a->comps.size()+b->comps.size());
		
			vector<polyn>::iterator belaend
			=std::copy ( a->comps.begin(),a->comps.end(),bel.begin() );
			std::copy	( b->comps.begin(),b->comps.end(),belaend	);		
		}

		init(fold,bel);

		if(comment=="")
		{
			ostringstream s;
			s<<"sum_{"<<pairs<<"}("<<a->comment<<")("<<b->comment<<")";
			comment=s.str(); 	
		}

	}

RTpartialFoldProd::~RTpartialFoldProd(){}

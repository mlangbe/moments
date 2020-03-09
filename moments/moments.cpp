/*
	implementations of operations on polynomial-valued tensors.

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


#include "stdafx.h"
#include<vector>
#include<map>
#include<cassert>
#include<iostream>
#include<algorithm>
#include<fstream>
#include<cmath>
#include<sstream>
#include <float.h>
#include<time.h>

#include "createSymm3DTens.h"
#include "symmtensor.h"
#include "tensor.h"
#include "polynom.h"
#include "optimizedPoly.h"
#include "vvpr.h"
#include "infbruch.h"
#include "bruch.h"
#include "myinterval.h"
#include "pointCloud.h"
#include "computeMoments.h"
#include "Momkdtree.h"
#include "momentSet.h"
#include "tensWithSymm.h"
#include "indepPolySet.h"
#include "optimizedFold.h"
#include "tensorGraphSplitTree.h"

#include "moments.h"

using namespace std;


struct tensmetainfo;








//create tensor with entries represented as ids
struct createTensId 
{

  tensor<polynom> x;

  static int nextid;

  int 
    d,
    so //number of indices in which tensor is symmetric
    ;

  int *inds,*inds2;//buffer var for symmetrize
  int depth;//buffer var for symmetrize, gives index in inds

public:
  operator tensor<polynom>& (){
    return x;
  }

	operator tensorWithSymm(){

		tensorWithSymm xx(x);
		if(so>0)
			xx.grp.push_back(so);
		
		if(x.ord()-so>0)
			xx.grp.push_back(x.ord()-so);

		xx.tgrp.push_back(xx.grp.size());

    return xx;
  }

  createTensId ( int order,int symmorder=0,int dim = 3 ,int id=0)
    :x(order,dim),so(symmorder),d(dim)
  { 
    if(!id)
      id=nextid;


    inditer ii(dim,order);
    do{
      x[*ii] = polynom(1,encodeIndices(*ii,order,dim,id)) ;  
    }while(++ii!=0);

    //set all entries whose #so first indices are not ordered
    //to those with the indices ordered
    if(symmorder>1){
      inds = new int [order];
      inds2= new int [order];
      depth = order-1; 
      symmetrize();
      cout<<endl;
      delete[]inds;
      delete[]inds2;
    }
    do{cout<<x[*ii]<<endl;}while(++ii!=0);




    nextid=id+1;
  }
  
  void symmetrize()
  {
    if(depth){
      for(int i=0;i<d;++i){
        inds[depth]=i;
        --depth;
        symmetrize();
        ++depth;
      }
    }
    else
      for(int i=0;i<d;++i){
        inds[depth]=i;
        copy(inds,inds+x.ord(),inds2);
        sort(inds2,inds2+so);

        cout<<"[";
        for(int j=0;j<x.ord();++j)
        {cout<<inds[j];}
        cout<<" -> ";
        for(int j=0;j<x.ord();++j)
        {cout<<inds2[j];}
        cout<<"]";

        x[inds]=x[inds2];        
      }    
  }



};

int createTensId::nextid=0;








bool fullRank(vector<vector<double> > & m,vector<int> &perm,double eps=1e-6)
{
  int nc=(int)m[0].size(),nr=(int)m.size();
  perm.resize(nc);

  for(int i=0;i<nc;i++)
    perm[i]=i;

  for(int i=0;i<nr;++i){
    double maxp=0;
    int rm=-1,rc=-1;

    //groesstes elem suchen
    for(int j=i;j<nr;++j)
      for(int k=i;k<nc;++k)
      {
        double fm=fabs(m[j][k]);
        if(fm>maxp)
          {rm=j;rc=k; maxp=fm;}
      }

        
    if(maxp<eps)return false;


    //zeilen vertauschen
    if(rm!=i)
      m[i].swap(m[rm]);
    
    //spalten vertauschen
    if(rc!=i){
      int x=perm[i];perm[i]=perm[rc];perm[rc]=x;
      for(int j=0;j<nr;++j)
        {double xch=m[j][i];m[j][i]=m[j][rc];m[j][rc]=xch;}
    }



    maxp=1./m[i][i];
    

    for(int j=i+1;j<nc;++j)
      m[i][j]*=maxp;
    m[i][i]=1;

    for(int j=i+1;j<nr;++j)
    {      
      double mm = m[j][i];
      m[j][i]=0;
      for(int k=i+1;k<nc;++k)
         m[j][k]-=m[i][k] * mm; 
    }
  }

  return true;
}





//apply permutations perm to paired index ids in x
void applyperm(const vector<int> & per, vpr &x)
{
  vpr::iterator i;
  for(i=x.begin();i!=x.end();++i)
  {
    i->a=per[i->a];
    i->b=per[i->b];
  }

  sort(x.begin(),x.end());
}

//check if there is one index range [i*granul,(i+1)*granul)
//where all index pairs in x stay in the same range
bool checkgranul(const vpr&x,int granul)
{
  if(granul < int(2*x.size()) )
  {
    vector<int> nump(x.size()*2/granul,0);

    for(int i=0;i<(int)x.size();++i){
      int ag= x[i].a/granul;
      int bg= x[i].b/granul;
      //if both indoces in same part: increase nump for this part
      if(ag==bg) ++nump[ag];

    }

    for(int i=0;i<(int)nump.size();++i)
    {
      if(nump[i] == granul/2)
        return true;
    }
  
  }

  return false;

}
//check if x is symmetric to an entry of y (using permutations in s)
bool checksym(const symm &s,const vpr&x,vvpr &y)
{

  sort(y.begin(),y.end());

  symm::const_iterator si;
  
  vpr permd;

  //try out all permutations if they lead to an element of y
  for(si=s.begin();si!=s.end();++si)
  {
    permd=x;
    applyperm(*si,permd);
    if(binary_search(y.begin(),y.end(),permd))
    { 
      vpr::const_iterator i;

      cout<<endl<<'<';
      for(i=x.begin();i!=x.end();++i)
        cout<<'('<<i->a<<' '<<i->b<<')';
      cout<<" <-> ";

      for(i=permd.begin();i!=permd.end();++i)
        cout<<'('<<i->a<<' '<<i->b<<')';
      cout<<"> ";

      vector<int>::const_iterator ii;
      cout<<'(';
      for(ii=si->begin();ii!=si->end();++ii)
        cout<<*ii<<' ';
      cout<<")"<<endl;
      
      return true;
    }
  }

  return false;

}


















/*
*get the pairings for a tensor that is a (tensor-product-)power
* of a total symmetric tensor.
* it leaves out pairings that can be computed already from 
* a lower power.(because there are sets of tensors unconnected)
*/
	struct getpairs2_s
	{
		int tot,ord,num;
		vvpr&ret;
		vpr vprbuf; //a buffer for addpairing
		
		//number of pairs inside the tensors
		vector<int> numinside;

		//number of pairs between the tensors
		vector<int> numbetween;

		//number of indices left in each tensor
		vector<int> numleft;

		//recursion depth
		int depth;

		//the possible pairings between num tensors
		vector< pair<int,int> > pairings;



		//groups of tensor that are still symmetric after having
		//paired inside the tensors (and so getting tensors of smaller order)
		//vector<vector<int> > stillsymm;
	public:
	
		getpairs2_s(int ord1,int numt,vvpr&r)
			:ret(r)
		{
			ord=ord1;num=numt;tot=ord*num;							
			assert(numt>0);
			//even total order required
			assert((tot&1)==0);

			if(numt==1)
			{
				for(int i=0;i<ord1;i+=2)
				{
					pr x={i,i+1};
					vprbuf.push_back(x);
				}
				r.push_back(vprbuf);
				return;
			}
			
			numinside.resize(num);
			numleft.resize(num);

			numbetween.resize(num*(num-1)/2);

			pairings.reserve(numbetween.size());
			for(int i=0;i<num;++i)
				for(int j=i+1;j<num;++j)
					pairings.push_back(pair<int,int>(i,j));
			
			depth=0;
			recurse();
		}
	
	private:

		void recurse()
		{			
			int i=0;
			if(depth)i+=numinside[depth-1];
			
			//we stop at ord be cause we don't want it totally inside one group
			for(int i=0; i<ord; i+=2)
			{
				numinside[depth]=i;
				if(depth<num-1)
				{ 
					//cout<<'('<<flush;
					++depth;
					recurse();
					--depth;
					//cout<<')'<<flush;
				}
				else 
				{
					//cout<<"numleft after inside pairing:";
					for(int j=0;j<numinside.size();++j)
					{		
						numleft[j]=ord-numinside[j];
						//cout<<numleft[j]<<' ';
					}
					//cout<<endl;
					//test if numleft is sane
					//(all indices will find their partner in other tensors)

					int maxnumleft = numleft[1];
					for(int j=2;j<num;++j)
						maxnumleft += numleft[j];

					//if the largest number of indices left in one
					//tensor is <= the number of indices left in all other tensors
					//(numleft[0] is the largest because we constructed it this way)
					if(numleft[0] <= maxnumleft )
					{

						int odepth=depth;
						depth=0;
						subrecurse();
						depth=odepth;
						//cout<<'!'<<endl;
					}
					//else cout<<'%'<<endl;
				}
			}		
		}




		//recurse over the interconnections (between tensors)
		//
		void subrecurse()
		{
			pair<int,int>& tids=pairings[depth];			
			int &a = numleft[tids.first],
					&b = numleft[tids.second];

			int oa=a,ob=b;

			if(depth == pairings.size()-1)
			{

				if(a==b){
					//pair all connections left
					numbetween[depth]=a;
					//test if all tensors satisfied:
					a=b=0;
					int i;
					for(i=0;i<numleft.size();++i)
						if(numleft[i]!=0)break;

					if(i==numleft.size()){
						if(num<=3 || testconnected() )
							addpairing();			
					}
				}
				else{	
				//	cout<<'$'<<flush;
				}

				a=oa,b=ob;
				return;
			}

			//determine max.number of inds in this connection
			int minnum = a;
			if(b<minnum)minnum=b;
			
			for(int i=0;i<=minnum;++i,--a,--b)
			{
				numbetween[depth]=i;
				++depth;
				//cout<<'<'<<i<<flush;
				subrecurse();
				//cout<<'>'<<flush;
				--depth;
			}
			a=oa,b=ob;
		}

		/** test if all 
		* tensors are connected
		*(so there are no unconnected sets)
		*/
		bool testconnected() const
		{
			// for every tensor the other
			//tensors it has a common pair with
			vector<vector<int> > x(num);

			for(int i=0;i<pairings.size();++i)
			{
				if(numbetween[i]==0)
					continue;

				int a=pairings[i].first,b=pairings[i].second;
				x[a].push_back(b);x[b].push_back(a);
			}

			//..................................

			//has tensor i already been reached ?
			vector<bool> reached(num,false);

			//list of elems reached,but not yet processed(neighbors added)
			vector<int> unprocessed ;
			unprocessed.reserve(num);

			
			//collect now all tensors indirectly reachable from the first
			unprocessed.push_back(0);
			reached[0]=true;

			int numreached=1;
			while(!unprocessed.empty())
			{					
					vector<int> & conns = x[ unprocessed.back() ];
					unprocessed.pop_back();
					
					vector<int>::iterator 
						i=conns.begin(),iend=conns.end();

					while(i!=iend)
					{
						if(!reached[*i])
						{	++numreached; unprocessed.push_back(*i); reached[*i]=true; }
						++i;
					}
			}
/*
			if(numreached<num){
				cout<<"\nreached:";
				for(int i=0;i<num;++i)
					cout<< (reached[i]);
				cout<<endl;
			}
*/
			//have we reached all tensors ?
			return (numreached==num);
		
		}


		//add a pairing to ret (compute it from numinside,numbetween)
		void addpairing(){
			

			//cout<<'Y'<<flush;
			vprbuf.clear();

			//add the pairings inside
			for(int i=0;i<num;++i)
			{
				int j=i*ord,jend=j+numinside[i];

				for(;j<jend;j+=2)
				{
					pr ppp={j,j+1};
					vprbuf.push_back(ppp);
				}
			}

			for(int i=0;i<num;++i)
				numleft[i]=ord-numinside[i];

			//add the pairings between tensors
			for(int i=0;i<pairings.size();++i)
			{
				int a=pairings[i].first,
						b=pairings[i].second;

				pr ppp={ord-numleft[a] +a*ord , ord-numleft[b] + b*ord};
				for(int j=0;j<numbetween[i];++j,++ppp.a,++ppp.b)
				{
					vprbuf.push_back(ppp);					
				}
				numleft[a]-=numbetween[i];numleft[b]-=numbetween[i];
				
			}

			ret.push_back(vprbuf);

		}


	};

	
void getpairs2( int ord1,
							 int numt ,vvpr&ret)
{

	getpairs2_s(ord1,numt,ret);
}



inline void testgp(int o,int m)
{

	vvpr ret;

	getpairs2(o,m,ret);
	cout<<endl;
	cout<<"order "<<o<<" times "<<m<<": "
		<<ret.size()<<" possibilities:"<<endl;	
	
	if(ret.size()>10){
		cout<<ret[0]<<endl;
		cout<<ret[1]<<endl;
		cout<<"..."<<endl;
		cout<<ret.back()<<endl;
	}else{
		for(int i=0;i<ret.size();++i)
			cout<<ret[i]<<endl;
	}


}

inline void testgetpairs2()
{
	testgp(2,2);
	testgp(2,3);

	testgp(3,2);
	testgp(3,4);
	testgp(3,6);

	testgp(4,2);
	testgp(4,3);
	testgp(4,4);
	testgp(4,5);

}


/**
* compute the possible distinct pairings
* of a tensor with ords.size() total symmetric index groups
* by goin through every possible number of connections between them
* and output a pair for every possibility
*/
	struct getpairs3_s
	{
		int num,tot;
		vvpr&ret;
		vpr vprbuf; //a buffer for addpairing
		

		//the orders/size of the index groups
		const vector<int>& ords;
		
		//the orders/size of the source tensor index groups
		vector<int>tords;

		//number of pairs between the tensors
		//(one entry for every entry of pairings)
		vector<int> numbetween;

		//number of indices left in each tensor
		vector<int> numleft;

		//buffer for addpairs
		vector<int> vibuf;

		//recursion depth
		int depth;

		//the possible pairings between num tensors
		vector< pair<int,int> > pairings;

 	public:
	
		getpairs3_s(const vector<int>& oords,vvpr&r,const vector<int>*otords)
			:ords(oords),ret(r)
		{
			if(otords)tords = *otords;
			else tords.push_back(ords.size());


			num=ords.size();
			tot=0;
			for(int i=0;i<ords.size();++i)
				tot+=ords[i];

			//even total order required
			assert((tot&1)==0);

			if(ords.size()==1)
			{
				for(int i=0;i<ords[0];i+=2)
				{
					pr x={i,i+1};
					vprbuf.push_back(x);
				}
				r.push_back(vprbuf);
				return;
			}
			


			numleft.resize(num);
			for(int i=0;i<num;++i)
				numleft[i]=ords[i];

			numbetween.resize(num*(num-1)/2);

			pairings.reserve(numbetween.size());
			for(int i=0;i<num;++i)
				for(int j=i+1;j<num;++j)
					pairings.push_back(pair<int,int>(i,j));
			
			depth=0;
			subrecurse();
		}
	
	private:

		//recurse over the interconnections (between tensors)
		//
		void subrecurse()
		{
			pair<int,int>& tids=pairings[depth];			
			int &a = numleft[tids.first],
				&b = numleft[tids.second];

			int oa=a,ob=b;


			//determine max.number of inds in this connection
			int minnum = a;
			if(b<minnum)minnum=b;

			for(int i=0;i<=minnum;++i,--a,--b)
			{
				numbetween[depth]=i;
				//cout<<'<'<<i<<flush;
				if(depth < pairings.size()-1)
				{
					++depth;
					subrecurse();
					--depth;
				}
				else
				{

					//test if all numleft's are even.
					//don't test if all are lower
					//than their initial value
					int j;
					for(j=0;j<numleft.size();++j){
						if(( numleft[j]&((int)1) ) !=0 )break;
						//if(numleft[j]== ords[j] )break;
					}

					if(j==numleft.size()){
						if( testconnected() )
							addpairing();			
					}
				}
			}

			//cout<<'>'<<flush;
			a=oa,b=ob;
		}

		/** test if all 
		* tensors are connected
		*(so there are no two unconnected sets of tensors)	
		*/
		bool testconnected() const
		{
			if(tords.size()==1)return true;
			//for every symmetric index group the index of the 
			//source tensor it is belonging to
			vector<int> srcperind(num);			

			//matrix with connections between src tensors
			vector<bool> conns(tords.size()*tords.size(),false);

			vector<int>::iterator src=srcperind.begin();
			for(int i=0;i<tords.size();++i)
				for(int j=0;j<tords[i];++j)
						*(src++) = i;

			// for every tensor the other
			//tensors it has a common pair with
			vector<vector<int> > x(tords.size());

			for(int i=0;i<pairings.size();++i)
			{
				if(numbetween[i]==0)
					continue;

				int a = srcperind[pairings[i].first], 
						b = srcperind[pairings[i].second];
			

				conns[a*tords.size()+b]=true;
				conns[b*tords.size()+a]=true;
			}

			//..................................

			//has tensor i already been reached ?
			vector<bool> reached(tords.size(),false);

			//list of elems reached,but not yet processed(neighbors added)
			vector<int> unprocessed ;
			unprocessed.reserve(tords.size());

			
			//collect now all tensors indirectly reachable from the first
			unprocessed.push_back(0);
			reached[0]=true;

			int numreached=1;
			while(!unprocessed.empty())
			{					
					vector<bool>::iterator 
						conn  = conns.begin()+unprocessed.back()*tords.size();

					unprocessed.pop_back();					

					for(int i=0;i<tords.size();++i)
					{
						if( conn[i] && !reached[i] )
						{	++numreached; 
							unprocessed.push_back(i); 
							reached[i]=true; 
						}
					}
			}
			/*
			if(numreached<tords.size()){
				cout<<"\nto discard:";
				for(int i=0;i<pairings.size();++i)
					cout<<" "<<numbetween[i]<<" x ("<<pairings[i].first<<' '<<pairings[i].second<<")";
				cout<<endl;
			}
			*/
			//have we reached all tensors ?
			return (numreached==tords.size());
		
		}


		//add a pairing to ret (compute it from numleft,numbetween)
		void addpairing(){
			

			//cout<<'Y'<<flush;
			vprbuf.clear();

			//a buffer storing the current index for each symmetric group
			vibuf.resize(num);

			vibuf[0]=0;
			for(int i=1;i<num;++i)
				vibuf[i]=vibuf[i-1]+ords[i-1];

			//add the pairings inside
			for(int i=0;i<num;++i)
			{
				int j=vibuf[i],jend=j+numleft[i];
				
				for(;j<jend;j+=2)
				{
					pr ppp={j,j+1};
					vprbuf.push_back(ppp);
				}
				vibuf[i]+=numleft[i];
			}


			//add the pairings between tensors
			for(int i=0;i<pairings.size();++i)
			{
				int a=pairings[i].first,
						b=pairings[i].second;

				pr ppp={vibuf[a], vibuf[b]};
				for(int j=0;j<numbetween[i];++j,++ppp.a,++ppp.b)
				{
					vprbuf.push_back(ppp);					
				}

				vibuf[a]=ppp.a; vibuf[b]=ppp.b;
				
			}

			ret.push_back(vprbuf);

		}


	};


//for the export(so we dont'have to export the struct in the headers)
void getpairs3(const vector<int>&ords,vvpr&ret,const vector<int>*tords=0)
{
	getpairs3_s(ords,ret,tords);
}


template<class T>
class vgen
{
	vector<T> buf;

public:
  vgen(){}
	vgen & operator,(const T&x)
	{ buf.push_back(x);return *this;}

	operator vector<T> & ()
	{return buf;}

};


class vgen0{
public:

	template<class T> 
	vgen<T> operator ,( const T&y)
	{
	  vgen<T> x;
		return (x,y);
	}

} vg0;


inline void testgetpairs3()
{
  vgen<vector<int> > vg,vg2;
	vector< vector<int> > testv = 
	( vg, 
	(vg0,2,4) ,
	(vg0,1,2,4,1), 
	(vg0,2,3,3), 
	(vg0,2,2,4,4),
	(vg0,2,4,2,4) )
	,
	testv2=(vg2,
	(vg0,2),
	(vg0,2,2),
	(vg0,3),
	(vg0,1,2,1),
	(vg0,2,2))
	;

	for(int i=0;i<testv.size();++i)
	{
		vvpr ret;

		getpairs3(testv[i],ret,&testv2[i]);

		for(int j=0;j<testv[i].size();++j)
			cout<<testv[i][j]<<"A ";
		cout<<":"<<endl;
		
		for(int j=0;j<ret.size();++j)
		{
			cout<<"  "<<ret[j]<<endl;
		}	
	}
}




void getpairs( int th ,vvpr&ret , 
              const symm * s=0, //if symm is given:exclude pairs symmetric to existing ones
              int granul=0 //if granul is given: exclude index parings where both indices are oin the same matrix
              //(and so this pair can be simulated by a lower-order one)
              )
{
  assert((th &1)==0);
  int nump=1;
  for(int x=th-1;x>0;x-=2)
    nump*=x;

  ret.resize(0);
  ret.reserve(nump);

  struct rec{
    vector<bool> used;
    vvpr * ret;
    vvpr::value_type vpt;
    vvpr::value_type::iterator depth;

    const symm * s;
    int granul;

    void recurse()
    {
      //find first free place
      int a=0;
      while(used[a])++a;      
      used[a]=true;     
      depth->a=a;      

      for(int b=a+1;  ;++b){
        //find next free place        
        while(b<(int)used.size() && used[b] )++b;
        
        //if no more free place, break
        if(b==used.size())break;

        depth->b = b;

        if(depth==vpt.end()-1)
        {

          if( (!s || !checksym(*s,vpt,*ret))  && (!granul || !checkgranul(vpt,granul)) ){            
            ret->push_back(vpt);
            //cout<<vpt<<endl;
          }
          else{
           // cout<< ".";
          }

          break;
        }
        else
        {

          used[b]= true;
          ++depth;
          //do the recursion
          recurse();
          //undo step:
          --depth;
          used[b]= false;        
          //----
        }
      }

      //undo step:
      used[a]=false;
      //----
    }
  
  } r;

  r.ret = &ret;
  r.s = s;
  r.vpt.resize(th/2);
  r.depth = r.vpt.begin();
  r.used.resize(th); 
  r.granul=granul;

  r.recurse();

}


int fak(int n)
{
  int ret=1;
  for(;n>0;--n)
    ret*=n;
  return ret;
}


struct createperm{

  void recurse(){   
      for(int i=0;i<(int)wh.size();++i){
        //find next free place
        while(i< (int)wh.size() && wh[i]>=0 ){
          ++i;
        }

        if(i>=(int)wh.size()) break;

        wh[i]=depth;
        ++depth;

        if(depth==wh.size()){
          
          si->resize(wh.size()*gs);
          for(int j=0;j<(int)wh.size();++j)
            for(int k=0;k<gs;k++)
             (*si)[ j*gs + k ] = wh[j]*gs +k;                   
            
          ++si;
          
        }
        else{
          recurse();
        }

        --depth;
        wh[i]=-1;
      }
    }

    vector<int> wh ;
    int depth,gs;
    symm::iterator si;

    createperm(symm&s, int groupsiz,int ngroups)
      :gs(groupsiz),wh(ngroups,-1)
    {
      s.resize( fak(ngroups) );
      si = s.begin();
      depth = 0;
      cout<<"o"<<endl;
      recurse();
    }

};


//test if the polynoms in p are independant
bool areIndependent(const vector<polynom> & p)
{
  vector<polynom>::const_iterator i;
  vvi::const_iterator j;
  vi::const_iterator k;

  //get parameters there:

  int maxid=0;
  for(i=p.begin();i!=p.end();++i)
    {int nmi=i->getMaxParamID();if(maxid<nmi)maxid=nmi; }

  vector<bool> paramThere(maxid+1);
  vector<int> params;
  params.reserve(maxid+1);
  
  for(i=p.begin();i!=p.end();++i)
    for(j=i->p.begin();j!=i->p.end();++j)
      for(k=j->begin();k!=j->end();++k)
        if(!paramThere[*k])
          params.push_back(*k), paramThere[*k]=true;

  if(params.size()<p.size())
    return false;

  //build test values:
  vector<double> paramValues(maxid+1);
  vector<int> primes(params.size());
  getPrimes(primes);

  for(int i=0;i<(int)params.size();i++)
    paramValues[params[i]] = sqrt( double(primes[i]) );


  vector<vector<double> > mat(p.size());

  for(int i=0;i<(int)p.size();++i)
  { 
    mat[i].resize(params.size());
    for(int j=0;j<(int)params.size();++j)
      mat[i][j] = p[i].evalderiv(paramValues,params[j]);
  }

  vector<int> perm;
  return fullRank(mat,perm);      
 
}

inline bool operator ! 
(const myinterval<long double> &x)
{return x.contains(0);}


void testdependencytest()
{

  polynom a[]={"ab+e","bc+e","abcde","x","y"};
  vector<polynom> x(a,a+5);
  cout<<(areIndependent(x)?"true":"false")<<endl;

  indepPolySet xx;
  for(int i=0;i<(int)x.size();++i)
    if(!xx.addpolynom(x[i]))
      cout<<i<<"dependant!"<<endl;
}

void testgetpairs()
{
 polynom x;

  vvpr ret;
  symm s;


  tensor<polynom> 
    p1 = createTensId(1,0,3)
    ,p2 = createTensId(2,1,3)
    ,p3 = createTensId(3,2,3)
    ;
  
  tensProd<polynom> 
    p2p2(p2,p2),
    p2p2p2(p2p2,p2);

  tensor<polynom> pi=p2;

  int granul=2;
  //  createperm(s,6,1);
  for(int np=2;np<=8;np+=2,pi*=p2){
    ret.clear();
    getpairs(np,ret,0,granul);

    map<polynom,vpr> retpolys;
  
    for(int i=0;i<(int)ret.size();++i)
    {
      fold<polynom> f(pi, ret[i]);
      pair<polynom,vpr> p( f ,  ret[i]);    

      cout<<p.second;
      pair<map<polynom,vpr>::iterator,bool> 
        pp = retpolys.insert(p);

      if(pp.second)
        cout<<':'<<endl<<p.first<<endl<<endl;
      else
        cout<<": ="<<pp.first->second<<endl;

    }


    cout<<"number:"<<retpolys.size()<<endl;
  }

}

/** get the dimension of the space of rotations for a specific dimension
*and a specific tensor order */
int numrotparms(int d,int o)
{  
  if(o==0)
    return 0;
  if(o==1)
    return d-1;
  
	return d * (d-1) / 2;


}

/**
*\param base: the tensor to compute rot.inv. moments from
*\param xx: a set of independant multivariate polynomes
*\param nmoments: number of moments before we are content
*\param accepted:
* the pair sets that led to a polynom that could be accepted
*/
void getRotInvPolys(const tens<polynom>& base,
										indepPolySet &xx,
										int nmoments=0,
										vvpr *
										accepted=0,ostream*out=0)
{
  vvpr ret;
	tensProd<polynom> basebase(base,base);
	
	const tens<polynom> * t = &base;
	//we need an even order:
  if((base.ord()&1)!=0)
    t = &basebase;
	

  tensProd<polynom> 
		pp(*t,*t),
		ppp(*t,pp),
		pppp(*t,ppp);

	const tens<polynom> *p[]={t,&pp,&ppp,&pppp };

  int nadded=0;
	for(int i=0;
		i<4  && 
		(nmoments==0||nadded < nmoments) ;++i){

			const tens<polynom> &pi=*p[i]
			;
    
    if(out)*out<< " getting possible index pairings order"<<pi.ord()<<endl;
		ret.clear();
		getpairs(pi.ord(),ret, 0,i>0? t->ord() : 0 );

    for(int j=0;j<(int)ret.size();++j)
    {
      if(out)*out<<j+1<<'/'<<ret.size()<< ": calculate polynome \r"<<flush;
      fold<polynom> f(pi, ret[j]);
      if(out)*out<<j+1<<'/'<<ret.size()<< ": testing independency\r"<<flush;
      if(xx.addpolynom(f))
			{
				if(accepted)
					accepted->push_back(ret[j]);
        ++nadded;
				if(out)*out<<endl<<ret[i]<<endl;
			}
    }
  }
}



void testpoly()
{
	
	cout<<"testing polynom"<<endl;
	polynom p("cc+aab+ee");

	polyn pp(p);
	
	infnum testv[5]={2,3,5,7,11};
	vector<infnum> test(&testv[0],&testv[0]+5);

	for(int i=0;i<5;++i)
		cout<<p.evalderiv(test,i).tofloat()<<endl;


	cout<<"integration:"<<endl;
	infnum upper=4;
	for(int i=0;i<5;++i)
		cout<<pp.evalint(test,i,upper).tofloat()<<endl;

}

/**
*for total symmetric tensors.
*\param base: the tensor to compute rot.inv. moments from
*\param xx: a set of independant multivariate polynomes
*\param nmoments: number of moments before we are content
*\param accepted:
* the pair sets that led to a polynom that could be accepted
*/
void getRotInvPolysSymm(const tens<polynom>& base,
										vector<indepPolySet> &xx,
										int nmoments=0,
										ostream*out=0,int maxord=0)
{
	cout<<"expecting max. "<<nmoments<<"moments"<<endl;
  vvpr ret;
	tensProd<polynom> basebase(base,base);
	
	const tens<polynom> * t = &base;
	//we need an even order:
  if((base.ord()&1)!=0)
    t = &basebase;
	

  tensProd<polynom> 
		pp(*t,*t),
		ppp(*t,pp),
		pppp(*t,ppp),
		ppppp(*t,pppp),
		pppppp(*t,ppppp);

	const tens<polynom> *p[]={t,&pp,&ppp,&pppp ,&ppppp,&pppppp};



  int nadded=0;
	for(int i=0;
		i<6 && 
		(nmoments==0||nadded < nmoments) 
		&& (maxord==0|| p[i]->ord()<=maxord);++i){

			const tens<polynom> &pi=*p[i]
			;
    
			if(out)*out<<endl<<"Getting possible index pairings of tensor order"<<pi.ord()<<endl;
		ret.clear();
		getpairs2(base.ord(),pi.ord()/base.ord(),ret);


		tensmetainfo mi;
		int id = base.ord(),num=pi.ord()/base.ord();
		mi.sourceTensors=vector<int>(num,id);

    for(int j=0;j<(int)ret.size()&&(nmoments==0||nadded < nmoments);++j)
    {

			mi.pairs=ret[j];

      polynom f = fold<polynom>(pi, ret[j]);
			f.setmetainfo(mi);

      if(out)*out<<j+1<<'/'<<ret.size()<< ": testing independency "<<ret[j]<<flush;


			int naccept=0;
			for(int k=0;k<xx.size();++k){
				naccept+= (int)xx[k].addpolynom(f);
			}
			if(naccept)
			{
        ++nadded;
				if(out){
					*out<<":"<<naccept<<"x independent."<<endl;
					for(int k=0;k<xx.size();++k)
					{
						//xx[k].printMat(*out);
						//*out<<endl;
					}
				}
			}
      else
        if(out)*out<<'\r'<<flush;
    }
  }
}

/**
* for tensors which are total symmetric in the index groups
* whose sizes are given in ords
*\param ords: the groups sizes of total symmetric indices in base
*\param base: the tensor to compute rot.inv. moments from
*\param xx: a set of independant multivariate polynomes
*\param accepted:
* the pair sets that led to a polynom that could be accepted
*/
void getRotInvPolysSymm2(tensWithSymm & base,
										vector<indepPolySet> &xx,
										int nmoments=0,
										vvpr *
										accepted=0,ostream*out=0)
{
	
	cout<<"expectin max. "<<nmoments<<"moments"<<endl;
	cout<<"Tensor Structure:";
	base.writeStructure(cout);
	cout<<endl;
			
  vvpr ret;
	
	tensProdWithSymm basebase(base,base);
	
	tensWithSymm * t = &base;
	
	//we need an even order:
  if((base.t().ord()&1)!=0)
    t = &basebase;

  tensProdWithSymm pp(*t,*t),ppp(*t,pp),pppp(*t,ppp);
	const tensWithSymm *p[]={t,&pp,&ppp,&pppp };

  int nadded=0;
	for(int i=0;
		i<4 && 
		(nmoments==0||nadded < nmoments) 
		&& p[i]->t().ord()<=14;
	++i){

			const tensWithSymm &pi=*p[i] ;
	
			cout<<" Tensor structure:";pi.writeStructure(cout);
    cout<< " getting possible index pairings..."<<flush;
		ret.clear();
		getpairs3(pi.groups(),ret,&pi.tgroups());
		cout<< "got"<<ret.size()<<" pairings"<<endl;

    for(int j=0;j<(int)ret.size() && (nmoments==0||nadded<nmoments) ;++j)
    {
      fold<polynom> f(pi, ret[j]);
			cout<<j+1<<'/'<<ret.size()<< ": testing independency..."<<flush;

			polynom fp=f;
			tensmetainfo mi;		
			mi.init(pi,ret[j]);
				
			
			fp.setmetainfo(mi);

			int naccept=0;
			for(int k=0;k<xx.size();++k)
			{
				//if(out!=0)*out<<"trying to add polynome..."<<flush;
				naccept += int(xx[k].addpolynom(fp));
				//if(out!=0)*out<<k<<":"<<naccept<<endl;
			}
			if(naccept){
				if(accepted)accepted->push_back(ret[j]);
				++nadded;
				cout<<":"<<naccept<<"x accepted = "<<nadded<<'/'<<nmoments<< endl;
				if(out)
					*out<<ret[j]<<endl;
				else
					cout<<ret[j]<<endl;
			}
			else
				cout<<'\r'<<flush;
    }
  }
	cout<<nadded<<"moments added"<<endl;
}




void createMomSetClass(ostream&out,
											 int id,
											 const string&valuetype,
											 const set<polynom>& p)
{
	int n=p.size();

	ostringstream os;
	os<<"momentSet<"<<valuetype<<","<<n<<","<<id<<">::";
	string prefix=os.str();

	out<<"static const "<<prefix<<"meta * "
		<<prefix<<"metainfo ={\n"
		<<" 35, 3,4,"
		<<" {\n";

	set<polynom>::const_iterator it;
	for(it=p.begin();it!=p.end();++it)
	{
		if(it!=p.begin())
			out<<" ,\n";

		if(!it->mi)
		{cout<<"no mi!"<<endl;out<<"{}";continue;}
		tensmetainfo &mi=*((tensmetainfo*)it->mi);
		int dord=0;
		for(int i=0;i<mi.sourceTensors.size();++i)
			dord+=mi.sourceTensors[i];
		dord += 3;

		out<<" {" << mi.sourceTensors.size() <<","<<dord 
			<<",\""<< mi<<" \" } \n";
	}
	out<<" }\n}\n";

	out<<"\nvoid "<<prefix<<"compute("<<valuetype
		<<"*M,const "<<valuetype<<"*A){\n";
		
	int i=0;
	for(it=p.begin();it!=p.end();++it,++i)
	{
		polyn pol(*it);
		optimizedPoly opol(pol);
		if(!it->mi)
		{cout<<"no mi!"<<endl; }
		else
			out<<"/*\n* "<<*it->mi<<"\n*/\n";
		out<<"M["<<i<<"]="<<opol<<";\n";
	}

	out<<"}\n";
	
}


void getMomentsScalarField(ostream&out,ostream& polysout)
{
	int maxord=12;

	tensor<polynom> p[]={
		createSymm3DTens(2,4),
		createSymm3DTens(3,4),
		createSymm3DTens(4,4)
	};


	typedef indepPolySet::mynum mynum;

	//example values for the polynomes
	vector<vector<mynum> >
		examplevals(2);


	for(int i=0;i<examplevals.size();++i)
	{
		vector<mynum> &v=examplevals[i];
		//v.resize(10000);
		for(int j=0;j<v.size();++j)
		{
			v[j] = (mynum)rand();
			//cout<<v[j]<<' ';
		}
		cout<<endl;
	}

	int nadded=0;
	/**maximum number of independant invariant moments*/
	int maxnum = -3;
	for(int i=0;i<3;++i)maxnum+=p[i].getNumDistinct();

	vector<indepPolySet> whole(examplevals.size());


	set<polynom> all;

		
	ostringstream des;

	for(int i=0;i<3;++i)
	{	
		vector<indepPolySet> pure(examplevals.size());
		for(int j=0;j<pure.size();++j)
			pure[j].setParamValues(examplevals[j]);

		des<<p[i].ord()<<"A :"<<endl;
		out<<p[i].ord()<<"A :"<<endl;
		cout<<p[i].ord()<<"A :"<<endl;
		
		getRotInvPolysSymm(p[i],pure,p[i].getNumDistinct()-3,&cout,maxord); 

		cout<<'\n'<<pure[0].p().size()<<" moments found. "<<endl;
		const indepPolySet::tpolyset &p=pure[0].p();
		indepPolySet::tpolyset::const_iterator it ;
		for(it=p.begin();it!=p.end();++it)
		{
			out<<"  "<<it->mi<<endl;
		}

		
		for(int j=0;j<examplevals.size();++j)
		{
			set<polynom>::const_iterator it;
			for(it=pure[j].p().begin();it!=pure[j].p().end();++it)
			{
				int numj=whole[j].addpolynom(*it);
				if(numj)
				{
					all.insert(*it);
				}
				else
					cout<<"x"<<flush;
			}
		}

		

		cout<<"\n#moments now:"<<whole[0].p().size()<<endl;
	}

	//the tensors to be multiplied
	int comps[] ={ 0,2, -1 , 0,0,2,-1, 0,0,0,2,-1, 0,2,2,-1,  0,1,1, -1, 0,0,1,1, -1,  1,1,2, -1,-1};

	int compid=0;
	vvpr pairs;
	while(comps[compid]>=0 && all.size()<maxnum)
	{
		vector<int> cids;
		vector<int> ords;
		while(comps[compid]>=0){
			cids.push_back(comps[compid]);
			ords.push_back(p[comps[compid]].ord());
			out<<ords.back()<<"A ";
			cout<<ords.back()<<"A ";
			++compid;
		}
		++compid;
		out<<endl; cout<<endl;


		tensProd<polynom> pp (p[cids[0]],p[cids[1]]);
		tensProd<polynom> * ppp=0, *pppp=0;
		
		if(ords.size()>= 3) ppp=new tensProd<polynom>(pp,p[cids[2]]);
		if(ords.size()>= 4) pppp=new tensProd<polynom>(*ppp,p[cids[3]]);

		const tens<polynom> *tp[]={0,0,&pp,ppp,pppp};

		pairs.clear();
		getpairs3(ords,pairs);

		tensmetainfo mi;
		mi.sourceTensors = ords;
		for(int k=0;k<pairs.size();++k)
		{
			polynom f = fold<polynom>(*tp[ords.size()],pairs[k]);
			mi.pairs=pairs[k];
			f.setmetainfo(mi);

			cout<<'\r'<<k+1<<"/"<<pairs.size()<<": testing "<<pairs[k]<<flush;
			int numj=0;
			for(int j=0;j<examplevals.size();++j)
					numj+= whole[j].addpolynom(f);				

			if(numj)
			{
				cout<<":"<<numj<<"x independent"<<endl;
				all.insert(f);
				out<<pairs[k]<<endl;
			}
			else
				cout<<"x"<<flush;
		}
		cout<<"\n#moments now:"<<whole[0].p().size();

		if(pppp)delete pppp;
		if(ppp)delete ppp;


			
	}  

	set<polynom>::iterator it;
	out<<"all polynomes:"<<endl;
	polysout<<" { ";
	for(it=all.begin();it!=all.end();++it)
	{

		out<<*it->mi;it->printForC(out);out<<endl;
		if(it!=all.begin())
		polysout<<" , ";
		polysout<<*it;		
	}
	polysout<<" } ";

	polysout.flush();out.flush();
	createMomSetClass(cerr,0,"mytype",all);

	ofstream f("momsetClassScalDim3Ord4.txt");
	createMomSetClass(f,0,"mytype",all);
	f.close();

	f.open("momScal_raw.txt");
	indepPolySet::save(f,whole[0].p());
	f.close();

}




void getMomentsVectorField(ostream&polysout,indepPolySet::tpolyset*outpoly=0)
{
  tensorWithSymm
    p1 = createSymm3DTens2(0,2,1)
    ,
    p2 = createSymm3DTens2(1,2,1)
    ,
    p3 = createSymm3DTens2(2,2,1)
    ;

	tensProdWithSymm
		p1p1(p1,p1),
		p1p1p2(p1p1,p2),
		p3p3(p3,p3),
		p2p3p3(p2,p3p3),
		p1p3(p1,p3),
		p1p3p2(p1p3,p2);

	int ndp1=p1.t().getNumDistinct(),
			ndp2=p2.t().getNumDistinct(),
			ndp3=p3.t().getNumDistinct();
			



	int numparams = ndp1+ndp2+ndp3;
  cout<<"num distinct params:"<<numparams<<endl;

  //number of parameters for rotation:
  int numrot = p1.t().dim() * (p1.t().dim()-1) /2;

  numparams-=numrot;
  cout<<"num expected params:"<<numparams<<endl;

	vector<indepPolySet::mynum> testv[2];
	for(int j=0;j<2;++j){
		testv[j].resize(numparams);
		for(int i=0;i<testv[j].size();++i)
			testv[j][i]=rand();
	}

  vector<indepPolySet> xzy(2),xx(2),yy(2),zz(2),xy(2),yz(2),xz(2),whole(2);
	for(int i=0;i<2;++i)
	{
		xx[i].setParamValues(testv[i]);
		yy[i].setParamValues(testv[i]);
		zz[i].setParamValues(testv[i]);
		xy[i].setParamValues(testv[i]);
		yz[i].setParamValues(testv[i]);
		xz[i].setParamValues(testv[i]);
		xzy[i].setParamValues(testv[i]);
		whole[i].setParamValues(testv[i]);
	}

	getRotInvPolysSymm2(p1,xx,1,0,&polysout); 
   cout<<"merging..."<<endl;
  for(int i=0;i<2;++i) whole[i].join(xx[i]);
	 getRotInvPolysSymm2(p2,yy,ndp2-3,0,&polysout); 
    cout<<"merging..."<<endl;


  whole[0].join(yy[0]);  
	getRotInvPolysSymm2(p3,zz,ndp3-3,0,&polysout); 
    cout<<"merging..."<<endl;
  for(int i=0;i<2;++i) whole[i].join(zz[i]);
  
	if(whole[0].p().size()<numparams){
    getRotInvPolysSymm2(p1p1p2,xy,ndp1+ndp2-3,0,&polysout), 
    cout<<"merging..."<<endl;
    for(int i=0;i<2;++i) whole[i].join(xy[i]);
	}
	if(whole[0].p().size()<numparams){
    getRotInvPolysSymm2(p2p3p3,yz,ndp2+ndp3-3,0,&polysout), 
    cout<<"merging..."<<endl;
		for(int i=0;i<2;++i) whole[i].join(yz[i]);
	}
	if(whole[0].p().size()<numparams){
    getRotInvPolysSymm2(p1p3,xz,ndp1+ndp3-3,0,&polysout), 
    cout<<"merging..."<<endl;
    for(int i=0;i<2;++i) whole[i].join(xz[i]);
	}
	if(whole[0].p().size()<numparams){
    getRotInvPolysSymm2(p1p3p2,xzy,ndp1+ndp3+ndp2-3,0,&polysout), 
    cout<<"merging..."<<endl;
    for(int i=0;i<2;++i) whole[i].join(xzy[i]);
	}

  cout<<"num elems:"
		<< whole[0].p().size()<<" "<<whole[1].p().size()<<endl;




/*
  indepPolySet xx,yy,zz,xy,yz,xz,whole;

	getRotInvPolys(p1,xx,1); 
   cout<<"merging..."<<endl;
   whole.join(xx);
  getRotInvPolys(p2,yy,p2.getNumDistinct()-3); 
    cout<<"merging..."<<endl;

  whole.join(yy);  
	getRotInvPolys(p3,zz,p3.getNumDistinct()-3); 
    cout<<"merging..."<<endl;

  whole.join(zz);

	if(whole.p().size()<numparams){
    getRotInvPolys(p1p1p2,xy), 
    cout<<"merging..."<<endl,
    whole.join(xy);
	}
  if(whole.p().size()<numparams)
    getRotInvPolys(p2p3p3,yz), 
    cout<<"merging..."<<endl,
    whole.join(yz);
  if(whole.p().size()<numparams)
    getRotInvPolys(p1p3,xz), 
    cout<<"merging..."<<endl,
    whole.join(xz);

  cout<<"num elems:"
      << whole.p().size()<<endl;
	*/

	if(outpoly)
		*outpoly = whole[0].p();
}

template<class T>
inline int countMults(const T&mytype)
{
		ostringstream oss;
		oss<<mytype<<ends;

		string s=oss.str();
		int nummult=0;
		for(const char* cc=s.c_str();*cc!=0;++cc)
			nummult += (*cc=='*');

		return nummult;
}

//read in polynomes, optimze them, write them.
inline void doOptimize()
{
	ifstream myin("polysScalarDim3Ord4.txt");
	ofstream myout("polysScalarDim3Ord4_opt.txt");
	ofstream myout2("polysScalarDim3Ord4_nopt.txt");
	char c;
	myin>>c;
	if(c!='{')
		throw;

	myout <<"\nvoid computeMomenteDim3Ord4_opt ( MOMTYP*M,const MOMTYP*A){\n"<<endl;
	myout2 <<"\nvoid computeMomenteDim3Ord4_nopt ( MOMTYP*M,const MOMTYP*A){\n"<<endl;

	int numOptmults=0,numunoptmults=0;
	int i=0;
	do{	
		polynom p;
		myin>>p>>c;
		polyn pp(p);
		cout<<"polynome with "<<pp.s.size()<<"summands: "<<endl;
		optimizedPoly opp(pp);
		myout<<"M["<<i<<"]= "<<opp<<";\n";
		myout2<<"M["<<i<<"]= "<<pp<<";\n";
		int mopp=countMults(opp),mpp=countMults(pp);
		cout<<"#multiplications: before optimize:"<<countMults(pp)<< " after: "
			<<mopp<<endl;

		numOptmults+=mopp;
		numunoptmults+=mpp;
		++i;
	}while(c==',' );

	myout <<"}"<<endl;
	myout2 <<"}"<<endl;

	myout.close();
	myout2.close();
	myin.close();

	cout<<"total number mults: before optim: "<<numunoptmults<< " after optim:"<<numOptmults<<endl;
}

bool compareVisually(const polyn&a,const polyn&b) 
	{
		list<polyn::summand>::const_iterator ait,bit;

		for(ait=a.s.begin(),bit=b.s.begin();
				ait!=a.s.end() && bit!=b.s.end();
				++ait,++bit)
		{
			bool different=(
				ait->c!=bit->c
				||
			  ait->f.size()!=bit->f.size());

			if(!different)
				for(int i=0;i<ait->f.size();++i)
					different |=
						(ait->f[i].p!=bit->f[i].p 
						|| ait->f[i].i!=bit->f[i].i );
			
			if(different)
			{
				cout<<"difference found :"<<endl;
				cout<<*ait<<" != "<<*bit<<endl<<endl;				

				return false;
			}

			
		}
		
		return ait==a.s.end() && bit==b.s.end();
	
	}

/*
//test if the 2 functions are equal
void testMomFuns()
{
	polyn A[35],M1[28],M2[28];

	for(int i=0;i<35;++i)
	{
		A[i]=polyn(1,i);
	}
	cout<<"1st comp..."<<endl;
	computeMomenteDim3Ord4_nopt(M1,A);
	cout<<"2nd comp..."<<endl;
	computeMomenteDim3Ord4_opt(M2,A);
	cout<<"comparing..."<<endl;
int i;
	for(i=0;i<28;++i)
		if(!compareVisually(M1[i],M2[i])){
			cout<<M1[i]<<"\n!=\n"<<M2[i]<<endl;
			cout<<"M1["<<i<<"] != M2["<<i<<"]"<<endl;
			break;
		}

  if(i==28)
	cout<<"seems to be correct"<<endl;

}
*/

void testMomFunSpeed()
{
	int num= 100*1000;
	static const int N=4;

	int numc=getAsiz(N);

	cout<<"number components:"<<numc<<endl;
	int numm=28;
	mytype *A=new mytype[numc*num];
	mytype *M=new mytype[numm*num];

	for(int i=0;i<numc*num;++i)
	{
		A[i]=rand();
	}

	
	clock_t st=clock();
	for(int j=1;j<=100&&(clock()-st) < 2*CLOCKS_PER_SEC ;j++)
	{
		mytype* Ai=A,*Aend= A+numc*num,*Mi=M;
		for(;Ai!=Aend;Ai+=numc,Mi+=numm){
			computeMomenteDim3Ord4_opt(Mi,Ai);
		}

		double nspt=1e9*(clock()-st)/CLOCKS_PER_SEC/j/num;
		cout<<"time per tset :"<< nspt	
		<<" ns time per component:"<<nspt/numc<<" ns"<<endl;
	}


	st=clock();
	for(int j=1;j<=100&&(clock()-st) < 2*CLOCKS_PER_SEC ;j++)
	{
		mytype* Ai=A,*Aend= A+numc*num,*Mi=M;
		for(;Ai!=Aend;Ai+=numc,Mi+=numm){
			computeMomenteDim3Ord4_nopt(Mi,Ai);
		}

		double nspt=1e9*(clock()-st)/CLOCKS_PER_SEC/j/num;
		cout<<"not opt: time per tset :"<< nspt	
		<<" ns time per component:"<<nspt/numc<<" ns"<<endl;
	}

}




inline void createTotalFolds()
{
	ofstream out("totalFoldPartialSpecs.txt");
	for(int i=1;i<=4;++i)
		createTotalFold(out,4,i,"double");
	out.close();
}

void testclustertree();

void testIntegrateTsets();

void momStatsmain(int argc, const char*argv[]);


inline void testOptPolyCommonSub()
{
	polyn aa[]={
		polynom("aa+bb+cca"),
		polynom("ab+bb+ca"),
		polynom("ab+b+ca")
	};

	vector<optimizedPoly> av(aa,aa+sizeof(aa)/sizeof(aa[0]));

	optimizedPolyCommonSubMulti ocs(av);

	cout<<ocs;

	vector<polyn> p(3),ret(3);
	p[0]=polyn(1,0);
	p[1]=polyn(1,1);
	p[2]=polyn(1,2);

	ocs.eval(p,ret);

	for(int i=0;i<av.size();++i)
		if(av[i].eval(p)!=ret[i])
		{cout<<av[i]<<" != "<<ret[i]<<endl;break;}

	cout<<"OK!"<<endl;




}



inline void testSavePolynom()
{
	polynom x("aa+bb+caa");

	tensmetainfo tmi;
	tmi.sourceTensors.resize(1,1);
	tmi.pairs.resize(1);
	tmi.pairs[0].b=1;
	x.setmetainfo(tmi);

	ofstream of("testsavepoly.txt");

	x.printRaw(of);
	of.close();
	ifstream in("testsavepoly.txt");
	polynom y;
	
	y.readRaw(in);



}


/**
convert a string to a representation
that re-creates the string when compiled by the c-compiler
*/
 
struct cify
{
	const char*x;
	cify(const char*xx):x(xx){}
};
ostream&operator<<(ostream&o,const cify&xx)
{
	const char*x=xx.x;
	o<<'"';
	for(;;)
	{
		switch(*x)
		{
		case 0:
			o<<'"';return o;
		case '\n':
			o<<"\\n\"\n\"";break;
		case '\r':
			o<<"\\r";break;
		case '\t':
			o<<"\\t";break;
		case '\\':
			o<<"\\\\";break;
		default:
			if(!isprint(*x))
				o<<"\\0"<<oct<<unsigned(*x);
			else
				o<<*x;			
		}
		++x;
	}
}



void createVecMomSpecializationHeader
(ostream&of,const vector<metainfo*> &mis,int nparams,int id,const char*comment)
{
	of<<
		"include\"momentSet.h\"\n"
		"typedef long double mytype;\n"
		"typedef momentSet<mytype,"<<mis.size()<<","<<id<<"> tmomset;\n"
		"\n"
		"const int tmomset::n;\n"
		"\n"
		"const tmomset::meta tmomset::metainfo ={\n"
		<<
		nparams<<",3,3,{\n";

	for(int i=0;i<mis.size();++i)
	{
		tensmetainfo &tmi=*((tensmetainfo*)mis[i]);
		int domord=0;
		for(int j=0;j<tmi.sourceTensors.size();++j)
			domord+=tmi.sourceTensors[j]-1;
		
		ostringstream oss;
		oss<<tmi;
		of<<"  {"<<tmi.sourceTensors.size()<<','<<domord<<","
			<<cify(oss.str().c_str())<<"}\n";
		if(i<mis.size()-1)
			of<<"  ,\n";
		else
			of<<" }\n";
	}
	if(comment)
	{
		of<<",\n"<<cify(comment)<<'\n';
	}
	of<<"};\n";
}


void createVecMomSpecialization(ostream& of,const indepPolySet::tpolyset&vecpolys,int id, bool optimize,bool docommonsubopt,const char*comment)
{

	vector<metainfo*> mis;
	
	vector<optimizedPoly> pols;
	
	for(indepPolySet::tpolyset::const_iterator it= vecpolys.begin();
			it!=vecpolys.end();
			++it)
	{
		mis.push_back(it->mi);
		pols.push_back(optimizedPoly(*it));		
	}

	optimizedPolyCommonSubMulti opcsm(pols);

	vector<polyn> ret1(pols.size()),ret2(pols.size()),
		params(opcsm.getMaxId()+1);

	for(int i=0;i<params.size();++i)
		params[i]=polyn(1,i);


	opcsm.eval(params,ret1);

	int i=0;
	for(;i<pols.size();++i)
	{
		ret2[i]=pols[i].eval(params);
		if(ret2[i]!=ret1[i])
		{
			cout<<"ERROR at number:"<<i<<endl;
			cout << ret2[i] << " != " << ret1[i] << endl;		
			break;
		}
	}
	
	if(i==pols.size())
		cout<<"opcsm was O.K!"<<endl;

	createVecMomSpecializationHeader(of,mis,params.size(),id,comment);
		
	of<<"void tmomset::compute(mytype*M,const mytype*A){\n";
	if(docommonsubopt)
		opcsm.print(of,&mis);
	else
	{
			indepPolySet::tpolyset::const_iterator it= vecpolys.begin();
			for(int i=0;i<pols.size();++i,++it)
			{
				if(mis[i])of<<"\n//"<<*mis[i]<<"\n\n";
				of<<"M["<<i<<"]=";
				if(optimize)
					of<<pols[i];
				else
					of<<polyn(*it);

				of<<";\n";									
			}
	
	}
	of<<"}\n";
}


//create the specialization 
//of momentSet for moments of vectorfield
//docommonsubopt: optimize for common subexpressions?
void createVecMomSpecialization(ostream& of,int id,bool optimize,bool docommonsubopt,const char*comment)
{
	indepPolySet::tpolyset vecpolys;

	ifstream mif("momvecraw.txt");

	if(mif.bad())
	{	
		getMomentsVectorField(cout,&vecpolys);
		ofstream mof("momvecraw.txt");
		indepPolySet::save(mof,vecpolys);
		mof.close();
	}
	else
	{
		vecpolys.clear();	
		indepPolySet::load(mif,vecpolys);
	}
	

	createVecMomSpecialization(of,vecpolys,id,optimize,docommonsubopt,comment);

}


/*
*print moments corresponding to the polynomes in tpolyset
*in einstein-sum-convention-format
*/
void printMomEsumForTex(ostream&out,const indepPolySet::tpolyset& vecpolys)
{	

	indepPolySet::tpolyset::const_iterator it;
	
	for(it=vecpolys.begin();it!=vecpolys.end();++it)
	{
		dynamic_cast<tensmetainfo*>(it->mi)->printEsumForTex(out);
	}

}

/*
*print moments corresponding to the polynomes in tpolyset
* as a picture of a graph
*/
void printMomGraphPS(ostream&out,const indepPolySet::tpolyset& vecpolys)
{	

	indepPolySet::tpolyset::const_iterator it;
	
	for(it=vecpolys.begin();it!=vecpolys.end();++it)
	{
		dynamic_cast<tensmetainfo*>(it->mi)->printGraphPS(out);
	}

}



void testOptimizedFold()
{
 tensorWithSymm
    p1 = createSymm3DTens2(0,2,1)
    ,
    p2 = createSymm3DTens2(1,2,1)
    ,
    p3 = createSymm3DTens2(2,2,1)
    ;

	refp<repeatTensor>
		rp1(new repeatTensor(p1,"p1")), 
		rp2(new repeatTensor(p2,"p2")),
		rp3(new repeatTensor(p3,"p3"));

	pr x[]={{0,1}};
	vpr pr1(x,x+1);

	pr x2[]={{0,2}};
	vpr pr2(x2,x2+1);

	pr x3[]={{0,3}};
	vpr pr3(x3,x3+1);

	pr x4[]={{2,5}};
	vpr pr4(x4,x4+1);

	pr x5[]={{1,4},{2,5},{0,6}};
	vpr pr5(x5,x5+2);


	vector<int> unitperm(20);
	for(int i=0;i<20;++i)
		unitperm[i]=i;

	refp<RTpartialFold>
		pf(new RTpartialFold(rp3,unitperm,pr1));
	refp<RTpartialFoldProd> 
		pfp(new RTpartialFoldProd(rp1,unitperm,rp2,unitperm,pr1));
	rp1->printForC(cout);
	rp2->printForC(cout);
	rp3->printForC(cout);
	pf->printForC(cout);
	pfp->printForC(cout);


	refp<RTpartialFoldProd> pfrp3
		=new RTpartialFoldProd(pf,unitperm,rp3,unitperm,pr1);
	pfrp3->printForC(cout);

/*

	RTpartialFoldProd pfrp4(pfrp3,rp3,pr2);
	pfrp4.printForC(cout);
	//((repeatTensor&)pfrp3).printForC(cout);

	RTpartialFoldProd pfrp5(rp3,rp3,pr3);
	pfrp5.printForC(cout);

	RTpartialFoldProd pfrp6(pfrp5,pfrp5,pr4);
	pfrp6.printForC(cout);

	RTpartialFoldProd pfrp7(pfrp6,pfrp6,pr5);
	pfrp7.printForC(cout);
	*/

}






void createTensorGraphAndCcodeForVectorMoms(	const indepPolySet::tpolyset& vecpolys )
{


	indepPolySet::tpolyset::const_iterator it;

	vector< vector<int> > gsizpertype(6);
	for(int i=1;i<6;++i)
	{
		vector<int> & g=gsizpertype[i];
		if(i>1)
			g.push_back(i-1);
		g.push_back(1);
	}




	tensorGraphSplitTree2 tgs2(gsizpertype,3);

	for(it=vecpolys.begin();it!=vecpolys.end();++it)
	{
		cout<<*it->mi<<endl;	

		tensmetainfo&tmi = *dynamic_cast<tensmetainfo*>(it->mi);

		tensorGraph tg;
		
		tg.nodelabels = tmi.sourceTensors;
		vector<int> endpointlabels;
		vector<int> tensorids;
		
		for(int i=0;i<tmi.sourceTensors.size();++i)
		{
			for(int j=0;j<tmi.sourceTensors[i]-1;++j)
				endpointlabels.push_back(0);		

			endpointlabels.push_back( tmi.sourceTensors[i]>1 ? 1: 0 );
			
			for(int j=0;j<tmi.sourceTensors[i];++j)
				tensorids.push_back(i);
		}

		for(int i=0;i<tmi.pairs.size();++i)
		{
			tensorGraph::edge e;
			pr&p=tmi.pairs[i];
			e.d.a= tensorids[p.a];
			e.d.b= tensorids[p.b];
			e.d.u= endpointlabels[p.a];
			e.d.v= endpointlabels[p.b];
			tg.edges.push_back(e);
		}

		tgs2.findOptimalSplit(tg);


	}


	vector<const repeatTensor*> tensPerType(6);
	for(int i=1;i<=3;++i)
	{
		tensorWithSymm
    p = createSymm3DTens2(i-1,2,1);
		tensPerType[i] = new repeatTensor(p);			
	}

	ofstream of("vectorgraphs.ps");
	ofstream ofc;

	tgs2.printPS(of);	
	
	ofc.open("vectorTens.txt");
	tgs2.printForC(ofc,tensPerType);
	ofc.close();

	of<<"\nshowpage\n";
	
	tgs2.joinSingleParentNodes();
	
	tgs2.printPS(of);

	ofc.open("vectorTens2.txt");
	tgs2.printForC(ofc,tensPerType);
	of.close();

}

void createTensorGraphAndCcodeForVectorMoms()
{
	ifstream mif("momvecraw.txt");

	indepPolySet::tpolyset vecpolys;
	
	vecpolys.clear();	
	indepPolySet::load(mif,vecpolys);
	createTensorGraphAndCcodeForVectorMoms(vecpolys);
}

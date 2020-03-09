/*
	some simple clustering algorithms.

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
#pragma once
#ifndef KMEANS_H
#define KMEANS_H
#include<cassert>
#include<vector>


#ifndef OLD_VC
#include<algorithm>
#define _TYPENAME_ typename
#else
#define _TYPENAME_
#endif


/*
#include "TsTree2.h"
struct kmeans3d{


	std::vector<T> ret;
	TsTree2 t;

	kmeans3d(std::vector<T>&ret_,const std::vector<mynum>& data,int k=0)
		:t(data)
	{
		ret.swap(ret_);

		if(k!=0){
			ret.resize(3*k)
			for(int i=0;i<k;++i)
			{
				int id = ((data.size()/3)*rand()/RAND_MAX)*3;
				std::copy(data.begin()+id,
					data.begin()+id+3,
					ret.begin()+3*i);
			}
		}







		ret.swap(ret_);
	}
};

*/


namespace cluster
{

template<class T>
struct eukliddist
{
	T operator()(const T*a,const T*b, int dim)
	{
		T ret=0;T buf;
		for(int i=0;i<dim;++i)
		{
			buf=a[i]-b[i];ret+=buf*buf;
		}
		return ret;
	}
};



/**helper for cluster_boxes*/
struct ilt
{
	int dim;
	ilt(int d):dim(d){}
	bool operator() (int*a,int*b)
	{
		return std::lexicographical_compare(a,a+dim,b,b+dim);
	}
};

template<class mynum>
	struct stats{

	private:
		mynum avg,var;
	
	public:
		mynum min,max;

		//number of values
		int n;

		stats(){n=0;};

		//add values
		void add(mynum v)
		{
			if(n==0){
				min=max=v; avg=var=0; 			
			}
			else
			{
				if(min>v)min=v; 
				if(max<v)max=v;
				avg+=v;var+=v*v;
			}
			++n;
		}

		void operator +=(const stats&o){
			if(n==0)*this=o;
			else{
				avg+=o.avg;var+=o.var;
				if(o.min<min)min=o.min;
				if(o.max>max)max=o.max;
				n+=o.n;
			}
		}
		
		mynum average() const {return avg/n;}
		mynum variance()const { return (var -avg *avg/n)/(n-1);}


	};

//returns number of clusters
template<class T>
int cluster_boxes(std::vector<int> & cellperpt,
									 const std::vector<T>& data,int dim,
									 //minimum number of clusters
									 int k
									 )
{
	using namespace std;
	int n = data.size()/dim;
	vector<int> idata(data.size());
	vector<int*> idataptr(n);
	cellperpt.resize(n);
	int* idataiter = &idata[0];
	for(int i=0;i<n;++i,idataiter+=dim)
		idataptr[i]=idataiter;

	_TYPENAME_
	vector<T>::const_iterator cit;
	_TYPENAME_
	vector<int>::iterator iit;

	vector< stats<T> > s(dim);

	cit=data.begin();
	for(int i=0;i<n;++i)
		for(int j=0;j<dim;++j,++cit)
			s[j].add(*cit);

	vector<T> avg(dim);

	T maxdist=0,buf;
	for(int i=0;i<dim;++i)
	{
		avg[i]=s[i].average();
		buf=s[i].max-avg[i];
		if(buf>maxdist)maxdist=buf;
		buf=avg[i]-s[i].min;
		if(buf>maxdist)maxdist=buf;
	}
		
	vector<T> sub(dim);

	T boxwidth = 2*maxdist;
	T dimsqrt2 = pow(.5,1.0/dim);
	int numcells;
	do{
		boxwidth *= dimsqrt2;

		cout<<"trying boxwidth="<<boxwidth<<flush;

		T m=1.0/boxwidth;
		for(int i=0;i<dim;++i)
			sub[i] = avg[i] - .5*boxwidth;

		cit=data.begin();iit=idata.begin();
		
		for(int i=0;i<n;++i)
		{
			for(int j=0;j<dim;++j,++cit,++iit)
				*iit = int(floor(m * ( *cit - sub[j] )) );
		}

		ilt iltd(dim);
		std::sort(idataptr.begin(),idataptr.end(),iltd);

		numcells=0;
		int i=0;
		cellperpt[( idataptr[0] - &idata[0] )/dim]=0;
		while(i<n){
			++i;
			while(i<n && !iltd(idataptr[i-1],idataptr[i]))
			{
				cellperpt[(idataptr[i] - &idata[0] )/dim]=numcells;				
				++i;
			}
			++numcells;		
		}		

		cout<<" numcells was "<<numcells<<endl;

	}while(numcells<k);


	return numcells;
}












//unoptimized kmeans with k*n*dim effort per iteration
template<class T,class d >
void kmeans_simple(std::vector<T> &ret, std::vector<int> & cellperpt,const std::vector<T>& data,int dim,int k=0,int maxiter=100,T mind=.1)
{
	using namespace std;

	d dist;

	_TYPENAME_
	vector<T>::iterator retit,minretit,oretit;
	_TYPENAME_
	vector<T>::const_iterator datit;
	
	assert(data.size()%dim==0);

	if(k>0)
	{
		ret.resize(dim*k);
		for(retit=ret.begin();retit!=ret.end();retit+=dim)
		{
			int id = ((data.size()/dim)*rand()/RAND_MAX);
			datit = data.begin()+(id*dim);
			copy(datit,datit+dim,retit);
		}
	}
	else
		k= ret.size()/dim;

	assert(k>0);
	assert(ret.size()%dim==0);
	assert(data.size()%dim==0);

	cellperpt.resize(data.size()/dim);
	vector<T> oret;
	vector<int> numpercell(k);
	for(int it=0;it<maxiter;++it)
	{

		cout<<"iteration "<<it+1<<"/"<<maxiter<<'\r'<<flush;
		oret=ret;		
		fill(ret.begin(),ret.end(),0);
		fill(numpercell.begin(),numpercell.end(),0);
		T buf,mindist;

		int dati=0;
		for(datit=data.begin();datit!=data.end();datit+=dim,++dati)
		{

			oretit=oret.begin();
			mindist = dist(&*datit,&*oretit,dim);
			minretit = oretit;
			for(oretit+=dim;oretit!=oret.end();oretit+=dim)
			{
				buf=dist(&*datit,&*oretit,dim);
				if(buf<mindist)
				{ 
					mindist=buf;
					minretit=oretit;
				}
			}

			//add current point to center of nearest ret.point
			int id=(minretit-oret.begin());
			retit=ret.begin()+id;
			for(int i=0;i<dim;++i)
				retit[i] += datit[i];			

			id/=dim;
			cellperpt[dati]=id;
			numpercell[id]++;
		}

		//norm new means by number of points in voronoi cell:
		cout<<"cells left of "<<k<<":";
		retit=ret.begin();
		for(int i=0;i<k;++i)
		{
			T m = 1.0 / numpercell[i];
			if(numpercell[i]!=0)
				cout<<i<<' '<<flush;

			for(int j=0;j<dim;++j,++retit)
				*retit *= m;
		}


		oretit=oret.begin();retit=ret.begin();
		T maxdist=0;
		for(int i=0;i<k;++i,oretit+=dim,retit+=dim)
			if(numpercell[i])
		{
			buf=dist(&*retit,&*oretit,dim);
			cout<<"distance"<<i<<":"
				<< sqrt(buf);										
			if(buf>maxdist)maxdist=buf;
		}

		if(maxdist<mind)
			break;

		cout<<endl;
	}
}



};//end namespace cluster

#undef _TYPENAME_
#endif

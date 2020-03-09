/*
	some basic  statistic methods for vectors
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

#ifndef VECSTATS_H
#define VECSTATS_H
#ifndef GNUCC
#pragma once
#endif

#include<ios>
#include<iostream>
#include<algorithm>
#include<float.h>
#include<tnt/tnt.h>
#include<tnt/jama_eig.h>
#ifndef OLD_VC
#include<stdexcept>
#endif



template<int n,class T=long double>
struct vecstats
{
	int num;
	T avg[n],sigma[n],min[n],max[n];
	T cov[n][n];

	void write(std::ostream&o) const
	{
		o<<'\n'<<n<<"\n\n";
		for(int i=0;i<n;++i)
			o<<avg[i]<<' ';
		o<<"\n\n";
		for(int i=0;i<n;++i)
			o<<sigma[i]<<' ';
		o<<"\n\n";
		for(int i=0;i<n;++i)
			o<<min[i]<<' ';
		o<<"\n\n";
		for(int i=0;i<n;++i)
			o<<max[i]<<' ';
		o<<"\n\n";
		for(int i=0;i<n;++i,o<<'\n')
			for(int j=0;j<n;++j)
				o<<cov[i][j]<<' ';
		o<<"\n"<<num<<'\n';
	}

	void read(std::istream&in)
	{
		int dummy;
		in>>dummy;
		if(dummy!=n)
		{
		 in.setstate(
		 std::ios::failbit
		 |std::ios::badbit
		 );
		 return;
		}
		for(int i=0;i<n;++i)
			in>>avg[i];
		for(int i=0;i<n;++i)
			in>>sigma[i];
		for(int i=0;i<n;++i)
			in>>min[i];
		for(int i=0;i<n;++i)
			in>>max[i];
		for(int i=0;i<n;++i)
			for(int j=0;j<n;++j)
				in>>cov[i][j];
		in>>num;

	}

	void initConvert(T convert[n][n]) const
	{
		  /* 
	Convert moments so that
	their covariance matrix is one:

  The covariance Matrix of a set of 
	vectors transformed with matrix M is
	C = M C' M^T, if C' is the covariance matrix of the untransformed ones.
	With an eigenvector decomposition
	C can be expressed as V D V^-1, with V being the orthonormal
	eigenvector matrix (here named ev)
	and D being a diagonal matrix holding the eigenvalues on the diagonals
	(the vector of diagonals is named ew in the code)
	D can be expressed as a product E E^T, with E being a diagonal matrix
	with E_ii = sqrt(ew_i)
	so it holds : C = V E I E^T V^T , with I being the unit matrix.
	If we set M= V E  and I=C' then we have C = M I M^T. 
	So, to transform from vectors with cov.-mat. C'=I to vectors 
	with cov. mat C you have to multiply with M. 
	To transform from vectors havin cov.mat. C to ones having cov.mat I,
	you have to transform with M^-1 = E^-1 V^-1  
	with V^-1=V^T, E^-1 being diagonal and E^-1_ii = 1/sqrt(ew[i]).
	*/
	
	
	for(int i=0;i<n;++i)
		for(int j=0;j<n;++j)
		{
			if(!finite(cov[i][j]))
			{
				cerr<<"Error: cov["<<i<<"]["<<j<<"]="<<cov[i][j]<<endl;cerr.flush();
				#ifndef OLD_VC
				throw std::runtime_error("non-finite values when calc. eigenvalues");
				#else
				throw new exception("non-finite values when calc. eigenvalues"); 
				#endif
			}
		}
	TNT::Array2D<T> a2d(n,n,const_cast<T*>(&cov[0][0]));
	JAMA::Eigenvalue<T> eigenv(a2d);

	TNT::Array2D<T> ev;
	TNT::Array1D<T> ew;

	eigenv.getRealEigenvalues(ew);
	eigenv.getV(ev);

	//for(int i=0;i<ew.dim();++i)
	//	if(dbgout)(*dbgout)<<"ew["<<i<<"]="<<ew[i]<<endl;

	

	for(int i=0;i<n;++i){
		int revi = n-1-i;
		T m = sqrt(1.0/ew[revi]);
		for(int j=0;j<n;++j)
			convert[i][j] = ev[j][revi]*m;
	}


	}


};

template<int n,class T>
std::ostream &operator << (std::ostream&o, const vecstats<n,T>&v)
{	v.write(o);return o; }

template<int n,class T>
std::istream &operator >> (std::istream&i, vecstats<n,T>&v)
{ v.read(i);return i;}



template<int n,class T=long double>
class vecstatscollector
{
	T sum[n],min[n],max[n];
	T cov[n][n];

	int num;	

public:


	vecstatscollector(){num=0;}

	void reset(){num=0;}

	void add(const T x[n])
	{
		if(num==0)
		{
			for(int i=0;i<n;++i)
			{
				std::fill_n(cov[i],n,0);
				sum[i]=0;
				min[i]=max[i]=x[i];
			}
		}
		
		for(int i=0;i<n;++i)
		{
			T xi=x[i];
			if(min[i]>xi)min[i]=xi;
			if(max[i]<xi)max[i]=xi;
			sum[i]+=xi;

			T*c=cov[i];

			for(int j=0;j<=i;++j)
				c[j]+=xi*x[j];			
		}	

		++num;
	}

	void add(const vecstatscollector &c)
	{
		if(num==0)
		{
			*this=c;
		}
		else{
		
			for(int i=0;i<n;++i)
			{
				if(min[i]>c.min[i])min[i]=c.min[i];
				if(max[i]<c.max[i])max[i]=c.max[i];
				sum[i]+=c.sum[i];

				T *cc=&cov[i][0];
				const T *nc=&c.cov[i][0];
				for(int j=0;j<=i;++j)
					cc[j]+=nc[j];			
			}	

			num+=c.num;
		}
	}


	operator vecstats<n,T>() const
	{
		return getstatscov(true);
	}

	vecstats<n,T> getstatscov(bool corr=false) const
	{	
		vecstats<n> s;
		for(int i=0;i<n;++i)
			s.avg[i]=sum[i]/num;

		std::copy(max,max+n,s.max);
		std::copy(min,min+n,s.min);

		if(num>0)
		{
			for(int i=0;i<n;++i){
				s.cov[i][i] =( cov[i][i]- s.avg[i]*s.avg[i]*num)/(num-1);
				s.sigma[i]=sqrt(s.cov[i][i]);
				if(corr)
					s.cov[i][i]=1;
			}


			for(int i=0;i<n;++i)
				for(int j=0;j<i;++j)
					s.cov[i][j] = s.cov[j][i]
					= (cov[i][j] - s.avg[i]*s.avg[j]*num)
						/
						( (num-1) * (corr ? s.sigma[i]*s.sigma[j] : 1) ) ;			
		}		
		else
		{
			std::fill_n(s.sigma,n,0);
			std::fill_n(&s.cov[0][0],n*n,0);			
		}

		s.num=num;

		return s;
	}

};

#endif

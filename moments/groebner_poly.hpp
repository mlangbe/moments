/*
	specializations for polynomial sets used in buchberger algorithm.

	
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


#include "groebner.h"

#include "poly.h"
/** specializations for poly<n,T>:
*/


template<int n, class T>
struct groebner_spec< poly<n,T> ,typename poly<n,T>::summand >
{


	static void normalize(poly<n,T>&ret){
		ret/=LT(ret).c;
	}


	/**
	*remove common factors from a,b
	*/
	static 
		void removeGCD(typename poly<n,T>::summand&a,typename poly<n,T>::summand &b)
	{
		unsigned char* ap=a.p.p,*bp=b.p.p,*apend=ap + n ;
		for(;ap!=apend;++ap,++bp)
		{
			if(*ap<*bp)
			{
				*bp-=*ap;
				*ap=0;
			}
			else 
			{
				*ap-=*bp;
				*bp=0;
			}
		}

	}


	static
		bool isConstant(typename poly<n,T>::summand &p)
	{
		const unsigned int* i =&p.p.i[0], *iend = i+ poly<n,T>::maxVarInt;
		for(;i!=iend;++i)
			if(*i!=0)
				return false;

		return true;
	}

	static
		bool isZero(const poly<n,T> &p)
	{
		return p.s.empty();
	}

	static
		void setZero(poly<n,T> &p)
	{
		return p.s.clear();
	}

	static
		const typename poly<n,T>::summand & LT(const poly<n,T> &p)
	{
		return p.s.back();
	}

	static
		poly<n,T>& wrapSummand(poly<n,T>&ret,const typename poly<n,T>::summand &p)
	{
		ret.s.clear();
		ret.s.push_back(p);
		return ret;
	}

	static
		void transferLT(poly<n,T>&e,poly<n,T>&rest)
	{
		poly<n,T> buf; 
		rest+=wrapSummand(buf,LT(e));
		e.s.pop_back();
	}
	static
		void removeLT(poly<n,T>&e)
	{
		e.s.pop_back();
	}
};



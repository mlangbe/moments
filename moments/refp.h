/*
	simple ref-counting class

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


#ifndef REFP_H
#define REFP_H

//a reference-countable object
class refCountable
{
	mutable int count;

protected:
	refCountable()
	{count=0;}

	template<class T>
	friend class refp;
public:
	int numref() const {return count;}

};


//simple refrence counting, not thread-safe,
//no subclasses allowed
template<class T>
class refp
{
	const T * p;

public:

	refp()
	{p=0;}
	
	refp(T*pp)
	{
		if(pp!=0)
			p=pp, ++p->count;
		else
			p=0;	
	}
	
	operator bool () const {return p!=0;}


	operator T* () const {return (T*)p;}

	void assign(T*pp)
	{
		++pp->count;

		if(p!=0)
		{
			--p->count;
		  if(!p->count) delete ((T*)p);
		}
		p = pp;			
	}

	template<class S>
	refp<T>& operator = (S * pp)
	{
		assign(pp);
		return *this;
	}

	refp<T>& operator = (T * pp)
	{
		assign(pp);
		return *this;
	}


	refp<T>& operator =(const refp<T>&pp)
	{
		++pp->count;

		if(p!=0)
		{
			--p->count;
		  if(!p->count) delete ((T*)p);
		}
		p = pp;		
		return *this;
	}

	~refp()
	{ 
		if(p)
		{
			--p->count;
			if(!p->count)
				delete ((T*)p);			
		}
	}		

	T* operator -> () const
	{
		return (T*)p; 
	}

	T& operator * () const
	{
		return *((T*)p); 
	}
};


/*
this would also work instead of using refCountable objects;
but would increase mem for all new ops
*/

/*
inline void* operator new(size_t x)
{
	size_t* pp = (size_t*) malloc(x+sizeof(size_t));
	*pp = 0;
	return pp+1;
}

inline void* operator new[](size_t x)
{
	size_t* pp = (size_t*) malloc(x+sizeof(size_t));
	*pp = 0;
	return pp+1;
}

inline void operator delete(void*x)
{
	size_t * pp =(size_t*)x;
	free(pp-1);
}

inline void operator delete[](void*x)
{
	size_t * pp =(size_t*)x;
	free(pp-1);
}
*/

#endif

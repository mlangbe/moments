/*
    hashmap for moment invariants.

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

#include<functional>
#include<algorithm>
#include<iostream>

#ifndef OLD_VC

namespace std{

template<typename T>
struct identity
{
	T operator()(const T&x)
	{return x;}
};

#define _NOTYPENAME_

}
#else

#define _NOTYPENAME_ typename

#endif












//an other hashset impl, copied a bit from java, but only allowing forward iterators.
template<class tkey,class tvalue,class thashvalue=std::identity<tkey>, class tequals=std::equal_to<tkey>  > 
class myhashmap
{

	/** max.fill percentage before tabsiz is doubled*/
	static float maxfill;

	/** min fill percentage: tabsiz is halved if below*/
	static float minfill;

private:

	typedef std::pair<const tkey,tvalue> T;

	struct entry 
	{
		inline entry(const T&kv,size_t h)
			:next(0),hash(h),values(kv)
		{
		}

		inline entry(const tkey&k,const tvalue&v,size_t h)
			:next(0),hash(h),values(k,v)
		{
		}
		
		entry * next;

		size_t hash;

		T values;

		friend class myhashmap<tkey,tvalue,tequals,thashvalue>;


	};

public:


	class const_iterator{
	
	protected:

		entry ** ctab, **tabend;
		entry *e;


		const_iterator(const myhashmap&m,bool end)
		{
			ctab = const_cast<entry**>(m.table); tabend = const_cast<entry**>(ctab)+m.tablesize;
			if(end)
			{
				e=0;
				ctab=tabend;
			}
			else
			{
				gotonextfilled();
			}
		}

		inline const_iterator(entry**ctab_,entry**tabend_,entry*e_)
			:ctab(ctab_),tabend(tabend_),e(e_)
		{
		}

	public:
		inline const_iterator(const const_iterator&ci)
			:ctab(ci.ctab),tabend(ci.tabend),e(ci.e)
		{
		}

	protected:


		friend class myhashmap;

		void gotonextfilled()
		{
		
				while( ctab<tabend && ! *ctab )
					++ctab;

				if(ctab<tabend)
					e=*ctab;
				else
					e=0;
		
		}

		void plusplus()
		{
			if(e!=0)
				e=e->next;

			if(e==0)
			{
				++ctab;
				gotonextfilled();
			}
		}

	public:

		inline const_iterator& operator ++()
		{
			plusplus();
			return *this;
		}


		inline const T & operator * () const
		{
			return e->values;
		}
		
		inline bool operator ==(const const_iterator&x) const
		{
			return e==x.e && ctab==x.ctab;
		}

		inline bool operator !=(const const_iterator&x) const
		{
			return e!=x.e || ctab!=x.ctab;
		}

	};


	class iterator:public const_iterator
	{
		#ifndef OLD_VC
		using const_iterator::plusplus;
		using const_iterator::e;
		#endif

		inline iterator(myhashmap&m,bool end)
			:const_iterator(m,end)
		{}
		
		inline iterator(entry**ctab,entry**tabend,entry*e)
			:const_iterator(ctab,tabend,e)
		{}


		friend class myhashmap;
	public:
		inline iterator& operator ++()
		{
			plusplus();
			return *this;
		}

		inline T & operator * () const
		{
			return e->values;
		}
	};




protected:


	typedef entry* pentry;

	const tequals equals;
	const thashvalue hash_value;


	entry** table;
	size_t tablesize;
	size_t size;



	static inline size_t hasher(size_t h)
	{
		h ^= (h >> 20) ^ (h >> 12);
        return h ^ (h >> 7) ^ (h >> 4);
	}


	entry** findpentry(entry**entpp,const tkey&x,size_t h) const
	{
		while(*entpp){
			entry & e=**entpp;
			if( h==e.hash && equals(x,e.values.first) )
				return entpp;

			entpp = &e.next;
		}
		return entpp;
	}

	entry* findentry(entry*e,const tkey&x,size_t h) const
	{
		while(e){
			if( h==e->hash && equals(x,e->values.first) )
				return e;

			e = e->next;
		}
		return e;
	}



	inline entry**findpentry2(const tkey&x, size_t h)
	{
		entry **entpp = table + ( h & (tablesize-1) );
		return findpentry(entpp,x,h);
	}

	inline const entry* const *findpentry2(const tkey&x, size_t h) const
	{
		entry **entpp = const_cast<entry**>(table) + ( h & (tablesize-1) );
		return findpentry(entpp,x,h);
	}
	

public:
	const_iterator find(const tkey&x) const
	{
		size_t h=hasher(hash_value(x));
		
		entry **tentpp = const_cast<entry**>(table) + ( h & (tablesize-1) );
		
		entry*entp=findentry(*tentpp,x,h);

		if(entp==0)
			return cend();
		
		return const_iterator(tentpp,table+tablesize,entp);
	}


	bool containsKey(const tkey& x) const
	{
		if(size==0)
			return 0;

		size_t h=hasher(hash_value(x));
		entry **tentpp = const_cast<entry**>(table) + ( h & (tablesize-1) );

		return findentry(*tentpp,x,h) !=0;
	}

	bool remove(const tkey&x)
	{
		if(size==0)
			return 0;

		entry**entpp = findpentry2(x,hasher(hash_value(x)));
		
		if(*entpp==0)
			return false;

		entry *entp= *entpp;
		entry *entpnext=entp->next;

		delete entp;
		
		*entpp = entpnext;	

		--size;

		if(size < tablesize* minfill )
			changetablesize(tablesize/2);

		return true;

	}


	bool insert(const T&x)
	{

		if(tablesize==0)
			changetablesize(4);

		size_t h=hasher(hash_value(x.first));

		entry**entpp = findpentry2(x.first,h);

		if(*entpp!=0)
			return false;

		*entpp=new entry(x,h);
	
		++size;

		if(size > tablesize*maxfill )
			changetablesize(tablesize*2);

		return true;
	}



	inline myhashmap()
	{
		size=0;table=0;tablesize=0;
	}

protected:
	inline iterator begin()
	{
		return iterator(*this,false);	
	}

	inline iterator end()
	{
		return iterator(*this,true);
	}
public:
	inline const_iterator cbegin() const
	{
		return const_iterator(*this,false);	
	}

	inline const_iterator cend() const
	{
		return const_iterator(*this,true);
	}


	void lengthHistogram(std::ostream&out)
	{
		static const int  maxlen=5;

		int sumlen=0;
		int maxl=0;
		int lens[maxlen]={0,0,0,0,0};

		for(entry**oents=table;oents!=table+tablesize;++oents){
			int len=0;
			entry *e= *oents;

			while(e){ sumlen+=len;++len; e=e->next;}
			if(len>maxl)maxl=len;
			lens[len>=maxlen?maxlen-1:len]++;
		}
		
		for(int i=0;i<maxlen;++i)
			out<<"length "<<i<<":"<< (lens[i]*100./tablesize)<< " %" << std::endl;

		out<< "mean path len:"<<sumlen*1.0/size<<"max path len:"<<maxl-1<<std::endl;

	}

private:

	void changetablesize(size_t newtabsiz)
	{
		if(newtabsiz==0)
		{
			tablesize=0;
			delete[] table;
			table=0;
			return;
		}

		entry**newtable = new pentry[newtabsiz];
		std::fill(newtable,newtable+newtabsiz,(pentry)0);

		if(table)
		{
			for(entry**oents=table;oents!=table+tablesize;++oents)
			{
				entry* oent = *oents;
				while(oent)
				{			
					entry** nent = newtable + ( oent->hash & (newtabsiz-1) );

					entry* oentnext=oent->next;
					oent->next = *nent;
					*nent = oent;

					oent=oentnext;
				}
			}
			delete[] table;
		}
		table=newtable;
		tablesize=newtabsiz;

		
	}


};

template<class tkey,class tvalue, class thashvalue,class tequals>
float myhashmap<tkey,tvalue,thashvalue,tequals>::maxfill=1.0f;

template<class tkey, class tvalue, class thashvalue,class tequals>
float myhashmap<tkey,tvalue,thashvalue,tequals>::minfill=
		myhashmap<tkey,tvalue,thashvalue,tequals>::maxfill*.25f;

template<class tkey, class thashvalue=std::hash<tkey> ,class tequals=std::equal_to<tkey> >
class myhashset
{
public:
	struct empty{};

	typedef myhashmap<tkey,empty,thashvalue,tequals> tdata;

	typedef typename tdata::const_iterator datcit;


	tdata data;



	class iterator: public datcit
	{
		
	public:

		inline iterator(const _NOTYPENAME_ datcit &x)
			:datcit(x)
		{}


		inline iterator(const iterator &x)
			:datcit(x)
		{}

		inline const tkey& operator *() const
		{
			return _NOTYPENAME_ datcit::operator*().first;
		}
	};	


	inline iterator begin(){ return iterator(data.cbegin()); }
	inline iterator end(){ return iterator(data.cend()); }


	inline iterator find(const tkey&k) const
	{
		return iterator(data.find(k));
	}

	inline bool contains(const tkey& k) const
	{
		return data.containsKey(k);
	}
	
	inline bool remove(const tkey&k)
	{
		return data.remove(k);
	}

	inline bool insert(const tkey&k)
	{
		return data.insert(std::pair<tkey,empty>(k,empty()) );
	}

};

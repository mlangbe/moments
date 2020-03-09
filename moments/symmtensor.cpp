/*
	operations to handle different symmetry types of tensors.
	(especially iterate over the non-repeated components)

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

#include<vector>
#include<map>
#include<algorithm>
#include<cassert>
#include<iostream>
#include<time.h>
#include"tensoriter.h"

#ifndef OLD_VC
#define _TEMPLATE_ template<>
#else
#define _TEMPLATE_
#endif



/** an iterator that
* iterates through all independant components
* in a symmetric tensor
* used primarily for building index map
*/
class symminditer
{
protected:
	//----------------------
	unsigned indbits;

	unsigned ordm1;

	unsigned dimm1;

	//dim is power of 2
	bool ispow2;
	//--------------------


 
	///the number of the component
	unsigned idx;

	///the indices of the tensor, sorted
	///(last one in least-significant bits)
	unsigned inds;

public:

	static int symmetrytype() {return SYMM;}

	///number of bits per index
	unsigned nbits;


	void init(unsigned dim, unsigned order)
	{
		ordm1=order-1; dimm1=dim-1;
		nbits=1;
		while(dimm1>>nbits)
			++nbits;
		indbits=(1<<nbits)-1;			


		//throw an err if not enough bits
		assert( ((unsigned)1)<<(order*nbits-1) !=0 );
		
		idx=0;inds=0;
	}

	symminditer()
	{nbits=0;}	

	/*use with care: it will not be checked if number+indexset	match*/
	void set(unsigned number,unsigned indexset)
	{ idx=number;inds=indexset;}

	unsigned getnumber() const { return idx;}

	unsigned getindexset() const
	{ return inds; }

	unsigned getindex(int id) const
	{	return ( inds >> ((ordm1-id)*nbits) ) & indbits;}

	//is the iterator at its maximum ?
	bool ismax() const
	{ return (inds>>(ordm1*nbits)) == dimm1;  }

	unsigned packinds(const unsigned*inds) const
	{
		unsigned ret=0;
		const unsigned*indsend=inds+ordm1;
		for(;inds!=indsend;++inds)
			ret<<=nbits,ret|=*inds;
		return ret;
	}

	void unpackinds(unsigned indexset,unsigned*inds) const
	{
		int nb= nbits*ordm1;
		unsigned *i= inds+ordm1+1;
		do{
			--i;
			*i = indexset&indbits;
			indexset>>=nbits;
		}while(i!=inds);

	}

	const symminditer& operator ++ ()
	{ 

		++idx;	
		unsigned nb = 0;

		//look how many indices are at max.value
		while( ((inds>>nb)& indbits ) == dimm1 )
			nb += nbits;

		//for sake of efficiency:
		if(nb==0)
		{ ++inds; return *this; }

		//increase the first one not at max value
		inds += 1<<nb;

		//now set all lower indices to 0
		inds &= ~ (( 1<<nb) -1 );

		//set all lower indices to the one just increased
		nb = inds & (indbits << nb);
		do{
			nb>>=nbits;
			inds|=nb;
		}while(nb);

		return *this;
	}

	const symminditer& operator --()
	{ 
		--idx;	
		//look how many indices are at min.value
		//= are at the value of the next upper index

		unsigned nb = 0,buf=(inds>>nbits)^inds;
		while((buf& indbits ) ==0 ) 
			nb+=nbits,buf>>=nbits;

		//for sake of efficiency:
		if(nb==0)
		{ --inds; return *this; }

		//decrease the first one not at min value
		inds -= 1<<nb;

		//now set all lower indices to 0
		inds &= ~ (( 1<<nb) -1 );

		//set all lower indices to dim-1
		nb = dimm1 << nb;
		do{
			nb>>=nbits;
			inds|=nb;
		}while(nb);

		return *this;
	}

};




class antisymminditer:public symminditer
{
public:
	static int symmetrytype() {return ANTISYMM;}

	void init(unsigned dim, unsigned order)
	{
		assert(order<=dim);

		ordm1=order-1; dimm1=dim-1;
		nbits=1;
		while(dimm1>>nbits)
			++nbits;
		indbits=(1<<nbits)-1;			

		//throw an err if not enough bits
		assert( ((unsigned)1)<<(order*nbits-1) !=0 );
		
		idx=0;
		inds=0;
		
		for(unsigned i=0;i<=ordm1;++i)
			inds=(inds<<nbits) | i;
	}
	
	const antisymminditer& operator ++ ()
	{ 

		++idx;	
		unsigned nb = 0;

		//look how many indices are at max.value
		while( ((inds>>nb)& indbits ) == dimm1 )
			nb += nbits;

		//for sake of efficiency:
		if(nb==0)
		{ ++inds; return *this; }

		//increase the first one not at max value
		inds += 1<<nb;

		//now set all lower indices to 0
		inds &= ~ (( 1<<nb) -1 );

		//set all lower indices to the one just increased,
		//+the index number (every one must be one bigger)
		unsigned d = (inds>>nb) & indbits;
		do{
			nb-=nbits;
			++d;
			inds|= (d<<nb);
		}while(nb);

		return *this;
	}

	const antisymminditer& operator --()
	{ 
		--idx;	
		//look how many indices are at min.value
		//= are at the value of the next upper index
		unsigned nb = 0,buf = (inds>>nbits)^inds;
		while( (buf & indbits ) == 1 	) 
			nb+=nbits, buf>>=nbits;

		//for sake of efficiency:
		if(nb==0)
		{ --inds; return *this; }

		//decrease the first one not at min value
		inds -= 1<<nb;

		//now set all lower indices to 0
		inds &= ~ (( 1<<nb) -1 );

		//set all lower indices to dim-1
		nb = dimm1 << nb;
		do{
			nb>>=nbits;
			inds|=nb;
		}while(nb);

		return *this;
	}

	bool ismax()
	{
		return (inds>>(ordm1*nbits)) == dimm1-ordm1;
	}

};




unsigned sortpackedinds(unsigned indices,int bitw)
{
	/**make insertion-sort */

	const unsigned 
		bw=bitw,
		bits=(1<<bw)-1;

	unsigned 
		inds = indices , 
		ret = inds&bits,
		nret = 1;

	inds>>=bw;


	while(inds)
	{

		unsigned a,b,x;	
		x = inds&bits;


		/*find bitpos where to insert x*/
		
		for(a=0,b=ret; b!=0; a+=bw,b>>=bw)		
			if( ( b & bits ) <= x )
				break;
			
		
		ret =   
			/* the bits after the position*/
			( ret & ((1<<a)-1) ) 
			/* the bits before the position, shifted by bw*/
			|( ( ret & ~((1<<a)-1) )  <<bw )	
			| (x<<a); /*the bits to insert*/

		inds>>=bw;
		++nret;
	}

	return ret;
}



struct indmap{

	int dim,ord,numcomp;
	indmap(){ dim=0;ord=numcomp=-1; }
	virtual int symmetrytype() const =0;
	virtual void getindices(unsigned number,unsigned*inds) const=0;
	virtual unsigned getnumber(const unsigned*inds) const=0;

};



struct idxmap:public indmap
{

	int symmetrytype() const {return ASSYMM;}

	idxmap(){}

	void init(int dimen,int order)
	{
		dim=dimen;ord=order;
		numcomp=dim;
		for(int i=1;i<ord;++i)
			numcomp*=dim;
	}

	idxmap(int dimen,int order)
	{init(dimen,order);}
	
	void getindices(unsigned number,unsigned*inds) const
	{
		unsigned*i=inds+ord;
		do{
			--i;
			*i=number%dim;
			number/=dim;
		}while(i!=inds);
	}
	
	unsigned getnumber(const unsigned*inds) const
	{
		unsigned ret=0;
		const unsigned *iend = inds+ord,*i;
		for(i=inds;i!=iend;++i){
			ret*=dim,ret+=*inds;
		}
		return ret;		
	}

	typedef std::map<int,std::map<int, idxmap > > tcache;
	static  tcache cache;


	static idxmap& getmap(int dim, int ord)
	{
		idxmap & im = cache[dim][ord];
		if(im.dim==0)
			im.init(dim,ord);	
		return im;
	}
	

};


template<class inditer>
struct indexmap:public indmap
{
	std::map<unsigned,unsigned> indices2number;
	std::vector<unsigned> number2indices;	
	
	int symmetrytype() const {return inditer::symmetrytype();}
	
	inditer iter;

	indexmap(){}

	void init(int dimen,int order)
	{
		dim=dimen;ord=order;
		iter.init(dim,ord);
		for(;;){
			number2indices.push_back(iter.getindexset());
			indices2number[iter.getindexset()]=iter.getnumber();
			if(iter.ismax())
				break;
			++iter;			
		}
		numcomp=number2indices.size();
	}


	void getindices(unsigned number,unsigned*inds) const
	{
		iter.unpackinds(number2indices[number],inds);
	}

	unsigned getnumber(const unsigned*inds) const
	{
		std::map<unsigned,unsigned>::const_iterator ret
			=indices2number.find
			( sortpackedinds(iter.packinds(inds),iter.nbits) );
		return 
			(ret==indices2number.end()? -1 : ret->second);
	}

	typedef std::map<int,std::map<int, indexmap<inditer> > > tcache;
	static  tcache cache;


	static indexmap<inditer>& getmap(int dim, int ord)
	{
		indexmap<inditer> & im = cache[dim][ord];
		if(im.dim==0)
			im.init(dim,ord);	
		return im;
	}

};




typedef indexmap<symminditer> mapsymm;
typedef indexmap<antisymminditer> mapantisymm;


idxmap::tcache idxmap::cache;
_TEMPLATE_
mapsymm::tcache mapsymm::cache;
_TEMPLATE_
mapantisymm::tcache mapantisymm::cache;


struct indextree:public indmap
{
	int	multitype;


	//a buffer for the indices that are symmetric
	//(when totalnum  is a multiple of num indices in children)
	//its size should be: totalnum/ sum_i( children[i].totalnum );
	mutable std::vector<unsigned> indbuf;

	/*the subtrees*/
	std::vector<indextree> children;

	/** how often is the subtre in children replicated ?*/
	unsigned multi;

	/*product of components of children,
	numcomp= pow(childrennumcomp,multi) */
	unsigned childrennumcomp;

	/** sum of al children orders, ord=multi*childrenord */
	unsigned childrenord;

	/** the mapping for that symmetr type*/
	indmap * mymap;


	int symmetrytype() const{return TREE;}

	indextree()
	{
		dim=0;numcomp=0;
	}
	

	void init(int dimen,int order,int typ)
	{
		dim=dimen,ord=order,multitype=typ,multi=1;
		childrennumcomp=1,childrenord=1;
		
		if(multitype!=ASSYMM)
		{
			if(multitype==SYMM)
				mymap=&mapsymm::getmap(dim,ord);
			else
				mymap=&mapantisymm::getmap(dim,ord);

			numcomp = mymap->numcomp;	
		}
		else
		{
			mymap=0;
			numcomp=dim;
			for(unsigned i=1;i<ord;++i)
				numcomp*=dim;
		}


	}

	indextree(int dimen, int order,int typ=0)
	{
		init(dimen,order,typ);
	}

	void init(std::vector<indextree>&it,int mult=1,int typ=0)
	{
		multi=mult;
		dim=it[0].dim;
		childrenord=it[0].ord;
		multitype=typ;
		for(unsigned i=1;i<it.size();++i)
		{ childrenord+=it[i].ord;assert(it[i].dim==dim);}		
		ord = childrenord * multi;		

		children.swap(it);

		std::vector<indextree>::iterator j;
		childrennumcomp=1;
		for(j=children.begin();j!=children.end();++j)
			childrennumcomp *= j->numcomp;

		mymap=0;
		if(multitype!=ASSYMM)
		{
			indbuf.resize(multi);
			if(multitype==SYMM)
				mymap=&mapsymm::getmap(childrennumcomp,mult);
			else
				mymap=&mapantisymm::getmap(childrennumcomp,mult);
			
			numcomp=mymap->numcomp;
		}
		else
		{
			numcomp = childrennumcomp;
			for(unsigned i=1;i<multi;++i)
				numcomp *= childrennumcomp;
		}

	}

	indextree(std::vector<indextree>&it,int multi=1,int typ=0)
	{init(it,multi,typ);}


	const char * init(int dimen,const char*x)
	{	
		dim=dimen,ord=0,multitype=0,numcomp=0;		
		children.clear();
		char typechar;

		assert(sscanf(x,"%c",&typechar));
		x++;
		int multi=1;
		int n=0;
		sscanf(x,"%d%n",&multi,&n);
		x+=n;

		switch(typechar){
			case 'n':multitype=ASSYMM;break;
			case 's':multitype=SYMM;break;
			case 'a':multitype=ANTISYMM;break;
			default: return 0;
		}	

		//subtrees are there
		if(*x=='(')
		{
			//overjump (
			++x;
			while( *x != ')' ){
				assert(*x!=0);
				children.push_back(indextree());
				x = children.back().init(dimen,x);							
			}
			//overjump )
			++x;

			init(children,multi,multitype);
		}
		else
			init(dimen,multi,multitype);

		return x;
	}
	
	indextree(int dimen,const char*x)
	{init(dimen,x);}


	/** convert the tensor  indices given in indices to 
	one index in an array containing the independent tensor components*/
	unsigned getnumber(const unsigned*indices) const
	{
		assert(dim!=0);
		unsigned b=0;
		/** a leaf node */
		if(children.size()==0)
		{			
			return mymap->getnumber(indices);
		}

		std::vector<unsigned>::iterator ib=indbuf.begin();
		std::vector<indextree>::const_iterator j;

		const unsigned * i=indices,*iend=indices+ord; 
		while(i!=iend)
		{
			
			for(j=children.begin();j!=children.end();++j)
			{
				b *= j->numcomp;
				b += j->getnumber(i);
				i += j->ord;						
			}			

			if(multitype!=ASSYMM)
			{
				*ib=b; b=0;++ib;				
			}
		}

		if(multitype==ASSYMM)
			return b;
		else
			return mymap->getnumber(&indbuf[0]);

	}

	/** convert the component number number 
	to the indices in indices*/
	void getindices(unsigned number,unsigned *indices) const
	{
		/** a leaf node */
		if(children.size()==0)
		{			
			if(multitype==ASSYMM){
				unsigned *i=indices+ord;
				do{	
					--i; 
					*i = number % dim; number /=dim;
				}while(i!=indices);
			}
			else
				mymap->getindices(number,indices);
		}
		else{
			std::vector<indextree>::const_iterator j;				
			std::vector<unsigned>::iterator k;

			unsigned * i = indices+ord;
			if(multitype!=ASSYMM){
				mymap->getindices(number,&indbuf[0]);
				k= indbuf.end();
			}
			do{					
				if(multitype!=ASSYMM)
					number = *(--k);

				j=children.end();
				do{ --j; i-=j->ord;
					j->getindices(number % j->numcomp , i );						
					number/= j->numcomp;
				}while(j!=children.begin());
			}while(i!=indices);			
		}

	}


	
	
};














void testsymmtensor()
{
	using namespace std;
	static const int ord=10,dim=4;

	/*
	int n=10<<20;


	cout<<hex<<sortpackedinds(0xFEDCBA98,4)<<dec<<endl;

	int test2[32]={1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,};

	cout<<"testing speed"<<endl;


	int entries = 8; 
	int nbits=32/entries;
	clock_t t0=clock();
	for( int i =0;i<n;++i)
		sortpackedinds(0xFEDCBA98 >>(32-nbits*entries) ,nbits);

	cout<<(clock()-t0)<< " for mine"<<endl;
	t0=clock();

	for( int i =0;i<n;++i)
		std::sort(test2,test2+entries);
	cout<<(clock()-t0)<< " for std"<<endl;



	mapsymm & mymap 
		= mapsymm::getmap(dim,ord);

	cout<<mymap.numcomp<<endl;

	unsigned n=mymap.numcomp;
	for(unsigned i=0;i<n;i++)
	{
		unsigned inds[20];
		mymap.getindices(i,inds);
		cout<<'<';
		for(int j=0;j<ord;++j)
			cout<<' '<<inds[j];
		cout<<'>';
	}
	cout<<endl;



	indextree a(3,"s2(s2)");	


	for(unsigned i=0;i<a.numcomp;++i)
	{
		unsigned inds[20];
		a.getindices(i,inds);
		cout<<'<';
		for(int j=0;j<a.ord;++j)
			cout<<inds[j];
		cout<<'>';
		if(i%4==3)cout<<'\n';
		
		if(i%12==11)cout<<'\n';
	}
	cout<<endl;

*/

	{
	symmii<3,2,2> a;
	symmii<3,3,2> b;
	symmii<3,4,2> c;
	symmii<3,5,2> d;


	while(!a.ismax())++a;
	while(!b.ismax())++b;
	while(!c.ismax())++c;
	while(!d.ismax())++d;
	
	cout<<a.idx+1<<" " <<b.idx+1<<" " <<c.idx+1 <<" " << d.idx+1 <<endl;

	return;

	}
	typedef symmii<2,5,1> ta;
	typedef assymmiicomposed<ta,ta> tb;
	typedef symmiicomposed< tb,2 > tc;
	
	typedef antisymmiicomposed< ta,4 > td;


	td xx;

	//cout<<hex<<xx.imax<<' '<<hex<<xx.imin<<endl;
	do
	{
		
		unsigned inds[td::ord];
		xx.getindices(inds);

		cout<<'<';
		//cout<<oct<<xx.inds<<'|';
		for(int j=0;j<xx.ord;++j)
		{cout<<inds[j];if(j%5==4)cout<<' ';}
		cout<<'>';

		if(xx.ismax()) break;
		++xx;

	}while(xx.idx<80);

	



}

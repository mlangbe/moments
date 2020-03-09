/*
    Object database using sets of moment invariants as representation.

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

#include<set>
#include<map>
#include<cmath>
#include<algorithm>
#include<numeric>
#include<time.h>
#include<iostream>
#include<tnt/jama_eig.h>
#include<array>
#include<functional>

#include "Objectdb.h"
#include "TsTree.h"
#include "tensorSet.h"
#include "momentSet.h"
#include "approxmath.h"
using namespace std;

#define _TYPENAME_ typename



std::ostream * &Objectdb::dbgout = outputflow::dbgout;

Objectdb::Objectdb(void)
{


}

Objectdb::~Objectdb(void)
{
}


void Objectdb::kdtreeEntryFromTset(
	mynum* ret,
	faceMomSet & bufm,
	const tensorSet3D<> &ts,
	mynum radius
	)
{
	bufm.compute(ts.A);
	bufm.scaleValue(1.0/ts.A[0]);
	bufm.scaleDomain(1.0/radius);

	for(int j=0;j<faceMomSet::n;++j)
	{
		mynum v=bufm.M[j];
		//take vord's root of v, 
		//so the histogram is less distorted
		v=(v>0?1:-1) * pow(fabs(v),
			mynum(1)/tmom::metainfo.me[j].vord);

		*(ret++)=v;
	}	
}



bool Objectdb::dosortouttsets=true;

//to help linker: keep local methods/functions in an own namespace even if not inlined.
namespace Objectdb_impl
{
	/** forwarding decls from Objectdb*/
	typedef Objectdb::mynum mynum;
	typedef Objectdb::tmom tmom;
	typedef Objectdb::momsetMeta momsetMeta;
	std::ostream* &dbgout= Objectdb::dbgout; 
	bool &dosortouttsets=Objectdb::dosortouttsets;



	/** now local decls*/
	typedef array<int,tmom::n> ivec;






	// help struct for sortouttsets
	struct a0lt{
		mynum m;
		a0lt(mynum mm):m(mm){}
		bool operator ()(tensorSet3D<>&t)
		{ return t.A[0]<m;	}
	};


	/** sort out tensor sets taken from the border of the object
	*(with less than the mean number of points)
	*/	
	static inline void sortouttsets(vector<tensorSet3D<> > & in,mynum factor=1)
	{
		mynum medpt=0;
		for(int i=0;i<in.size();++i)
			medpt+=in[i].A[0];
		medpt/=in.size();
		in.resize( std::remove_if(in.begin(),in.end(),a0lt(medpt*factor))
			-in.begin() );
	}


	/*
		get the invariant moments for 
		balls with radius ra around all
		leaves of the octree of the points of the pointcloud
		describing the object
		\retval m
			a vector containing for every ball the invariant moments
		\retval pt 
			a vector containing the center of gravity for every pointset
			from which a moment set was computed
		\retval stats:
		 a set of statistics (variance, mean,...)
		 for every component of the moment sets
		 \param t: 
		 an octree of the pointcloud having the tensorset
		 0A,1A,2A,3A,4A of the points inside the node 
		 stored in each node
		 \param ra:
			the radius of the balls
		 \param minh:
			minimum height at which to stop the 
			summing up of the tensorsets for each ball
	*/
	static int getMoms(
		//the invariant moments
		vector<mynum>& m,
		//the corresponding points
		vector<mynum>& pt,
		//statistik
		vecstatscollector<tmom::n,mynum> &s,
		//Baum, der die tensorsets enthaelt				
		const TsTree&t,
		//radius of the balls
		mynum ra,
		//minimum height to go down in the tree
		int minh=0
		)
	{
		vector<tensorSet3D<> > in,bo;
		if(dbgout)(*dbgout)<<"getting tsets..."<<flush;
		t.numbersummed=0;t.ptssummed=0;
		t.getBallTsets(in,bo,minh,(int)(ra/t.getScale()));
		
		if(dbgout)(*dbgout)<<"ready: #tsets summed:"<<t.numbersummed<<" #pts summed:"<<t.ptssummed;

		if(dbgout)(*dbgout)<<"computing "<<in.size()<<" moments ..."<<flush;
		tmom bufm;

		if(dosortouttsets)
			sortouttsets(in);
		/**/
		if(dbgout)(*dbgout)<<bo.size()<<" before "<<in.size()<<" after "<<endl;

		m.resize(in.size() * faceMomSet::n);
		pt.resize(in.size()* 3);
		vector<mynum>::iterator mit=m.begin(),ptit=pt.begin();
		for(int i=0;i<in.size();++i,mit+=tmom::n,ptit+=3)
		{			
			in[i].getCenter(&*ptit);
			ptit[0]=-ptit[0];ptit[1]=-ptit[1];ptit[2]=-ptit[2];
			in[i].translate(&*ptit);
			ptit[0]=-ptit[0];ptit[1]=-ptit[1];ptit[2]=-ptit[2];

			Objectdb::kdtreeEntryFromTset(&*mit, bufm ,in[i], ra);
			s.add(&*mit);
		}	

		if(dbgout)(*dbgout)<<"ready. "<<flush;
		return in.size();

	}




	/** distance in maxnorm */
	static inline mynum maxd(const mynum*a,const mynum*b,int n)
	{		
		mynum md=0;
		for(int i=0;i<n;++i)
		{ 
			mynum d=fabs(a[i]-b[i]);
			if(md<d)md=d;
		}
		return md;
	}

	/** distance in 2-norm(Euclidean */
	static inline mynum eucd(const mynum*a,const mynum*b,int n)
	{		
		mynum md=0,d;
		for(int i=0;i<n;++i)
		{ 
			d=a[i]-b[i];
			md+=d*d;
		}
		return sqrt(md);
	}



	/*
	*throw out all which have different object 
	*id's but are too close(closer than 1.0/sc[i])
	*\param m  the moment sets
	*\param metas:
	*for every moment set the object it belongs to
	*\param s: for every moment the statistics
	*\param sc: the scaling factor for every 
	* moment before rounding
	*\param added: the remaining moment sets
	*of the remainin moment sets
	*\param verycommon: if not zero, 
	* points to an array which will be filled 
	* with moment id's which are very common among the objects
	* (belongs to more than 20 objects)
	*\param supportmulti:
	* leave one moment for every Object in the bucket.
	*/
	static void throwOutMoments(	
		const vector<mynum> & m,
		const vector<momsetMeta> & meta,
		const mynum sub[tmom::n],
		const mynum sc[tmom::n],
		vector<int>& added,
		vector<int>*verycommon=0,
		bool supportmulti=false
		)
	{
		typedef map< vector<int> ,vector<int> > mymap;

		int nmomsets=m.size()/tmom::n;

		vector<int> buf(tmom::n);
		mymap momi;

		added.resize(m.size()/tmom::n);
		for(int i=0;i<added.size();++i)
			added[i]=i;

		if(dbgout)(*dbgout)<<added.size()<<" moment sets before."<<endl;
	
		vector<int>::const_iterator ait;
		//throw out moments in 29 steps
		for(mynum delt=0;delt<1;delt+=1.0/(tmom::n+1))
		{
			//fill map:
			for(ait=added.begin();ait!=added.end();++ait)
			{
				vector<mynum>::const_iterator 
					mit	= m.begin() + *ait * tmom::n;

				
				for(int j=0;j<tmom::n;++j,++mit)
				{	
					mynum v= (*mit - sub[j] );
					v*= sc[j];
					//if(dbgout)(*dbgout)<<v<<' '<<flush;
					buf[j]=(int)floor(v+delt);			
					//if(dbgout)(*dbgout)<<buf[j]<<' '<<flush;
				}

				momi[buf].push_back(*ait);
			}

			added.clear();
			
			mymap::iterator mapit;
			set<int> objbuf;
			for(mapit=momi.begin();mapit!=momi.end();++mapit)
			{
				objbuf.clear();
				vector<int>::const_iterator momids;
				vector<int> & ids = mapit->second;
				
				int nobj=0;//number of objects in current bucket
				for(momids= ids.begin();momids!=ids.end();++momids)
					//if a new object is encountered:
					if(objbuf.insert(meta[*momids].objectid).second)
						if(supportmulti)
							added.push_back(*momids);
						else
							++nobj;

				
				if(!supportmulti && nobj==1)
					added.push_back(ids.front());												

				if(verycommon && nobj>20)
					verycommon->push_back(ids.front());


			}
			//free old map, so we can reuse it
			momi.clear();

			if(dbgout)(*dbgout)<<added.size()<<" moment sets afterwards."<<endl;
	
		}

	}
	

	
	/**helper for throwOutMoments2*/
	struct ilt
	{
		ilt(){}
		bool operator() (int*a,int*b)
		{
			return lexicographical_compare(a,a+tmom::n+1,b,b+tmom::n+1);
		}
	};

	/**helper for throwOutMoments2*/
	struct ieq
	{
		ieq(){}
		//compare first tmom::n entries for equality
		bool operator() (int*a,int*b)
		{
			int*aend=a+tmom::n;
			while(a!=aend){
				if(*a!=*b)return false;				
				++a;++b;
			}
			return true;
		}
	};



/*
	*throw out all which have different object 
	*id's but are too close(closer than 1.0/sc[i])
	*\param m  the moment sets
	*\param metas:
	*for every moment set the object it belongs to
	*\param sub: value to subtract from every moment
	* for every moment the statistics
	*\param sc: the scaling factor for every 
	* moment before rounding
	*\param added: the remaining moment sets
	*of the remainin moment sets
	*\param verycommon: if not zero, 
	* points to an array which will be filled 
	* with moment id's which are very common among the objects
	* (belongs to more than 20 objects)
	*/

	static void throwOutMoments2(		
		const vector<mynum> & m,
		const vector<momsetMeta> & meta,
		const mynum * sub,
		const mynum * sc,
		vector<int>& added
		)
	{
		int nmomsets=m.size()/tmom::n;

		//an array with the rounded moments plus the moment id as last entry
		vector<int> imoms(m.size()+nmomsets);
		//pointers on the integer moment sets in imoms
		added.resize(nmomsets);

		vector<int*> imomlist(nmomsets);

		for(int i=0;i<added.size();++i)
		 added[i]=i;


		if(dbgout)(*dbgout)<<nmomsets<<" moment sets before."<<endl;
	
		//throw out moments in three steps
		for(mynum delt=0;delt<1;delt+=0.25)
		{
			
			imomlist.resize(added.size());
			
			//fill rounded moments
			vector<int>::iterator imomsit= imoms.begin();
			vector<int*>::iterator imomlistit=imomlist.begin();
			vector<int>::const_iterator ait;
				
			for(ait=added.begin();	
					ait!=added.end();
					++ait)
			{
				vector<mynum>::const_iterator 
					mit	= m.begin() + *ait * tmom::n;
				*(imomlistit++) = &(*imomsit);

				if(sub&&sc)
					for(int j=0;j<tmom::n;++j)
					{	
						mynum v= (*(mit++) - sub[j] );
						v*= sc[j];
						//if(dbgout)(*dbgout)<<v<<' '<<flush;
						*(imomsit++)=(int)floor(v+delt);			
						//if(dbgout)(*dbgout)<<buf[j]<<' '<<flush;
					}
				else
					for(int j=0;j<tmom::n;++j)
					{	
						*(imomsit++)=(int)floor(*(mit++) + delt );			
					}
					*(imomsit++) = *ait;
			}

			sort(imomlist.begin(),imomlist.end(),ilt());
			
			added.clear();			
			vector<int*>::iterator mapit;
			mapit=imomlist.begin();
			while(mapit<imomlist.end())
			{
				int * lasti= *mapit;				
				int firstobj= meta[lasti[tmom::n]].objectid;
				++mapit;
				int nobj=1;				
				while(mapit<imomlist.end() && ieq()(lasti,*mapit))
				{
					nobj+=( meta[ (*mapit)[tmom::n] ].objectid != firstobj);
					mapit++;
				}
				if(nobj==1)
					added.push_back(lasti[tmom::n]);												
			}
			if(dbgout)(*dbgout)<<added.size()<<" moment sets afterwards."<<endl;

		}

	}






	static void clusterMoms(
		const vector<mynum> &m,
		const mynum*sub,
		const mynum*sc,
		vector<int>&added
		)
	{
		int nmomsets=m.size()/tmom::n;

		vector<int> buf(tmom::n),buf2(tmom::n);
		set<vector<int> > momi,momi2;
	
		int num0=0, num1=0,num2=0;
		vector<mynum>::const_iterator mit=m.begin();
		for(int i=0;i<m.size()/tmom::n;++i)
		{
			bool is0=true;
			bool is1=true,is2=true;
			//scale down by std. deviation
			for(int j=0;j<tmom::n;++j)
			{	
				mynum v= (*mit - sub[j] );
				v*= sc[j];
				//if(dbgout)(*dbgout)<<v<<' '<<flush;
				buf[j]=(int)floor(v+0.5);			
				buf2[j]=(int)floor(v);
				//if(dbgout)(*dbgout)<<buf[j]<<' '<<flush;
				++mit;
				is0 &= (fabs(v)<0.5);
				is1 &= (fabs(v)<1);
				is2 &= (fabs(v)<2);
			}
			
			num0+= is0;
			num1+= is1;
			num2+= is2;

			//cluster it:
			//only take those moments which are not yet in the sets
			//momi,momi2
			if( momi.insert(buf).second && momi2.insert(buf2).second)
			{
				added.push_back(i);
			}
			
		}

		if(dbgout)(*dbgout)<<"clustered it :"<<added.size()<<" of "
			<<nmomsets<<" left"<< endl;
		
		if(dbgout)(*dbgout)<<num0<<"were < 1/2 std deviation (in all dims)"<<endl;
		if(dbgout)(*dbgout)<<num1<<"were < 1 std deviations (in all dims)"<<endl;
		if(dbgout)(*dbgout)<<num2<<"were < 2 std deviations (in all dims)"<<endl;

	}





	struct cmpindices
	{
		const vector<mynum> & value;
		cmpindices(const vector<mynum>& val)
			:value(val){}

		bool operator() (int	a,int b)
		{ return value[a]<value[b];}
	};



};

using namespace Objectdb_impl;





void Objectdb::momPerRadius::save(FILE*f)
{

	uint32 i;
	getMomSets();
	uint32 buf=meta.size();
	fwrite(&buf,sizeof(buf),1,f);	
	if(buf>0){
		outputflow::out()<<"writing "<<meta.size()<<" moment invariant sets with meta info "<<endl;
		fwrite(&momentSets[0],sizeof(momentSets[0]),momentSets.size(),f);
		fwrite(&meta[0],sizeof(meta[0]),meta.size(),f);
	}
		
	int nmultiObj= objectsnearby.size();
	fwrite(&nmultiObj,sizeof(nmultiObj),1,f);
	for(tobjectsnearby::const_iterator it=objectsnearby.begin();it!=objectsnearby.end();++it)
	{
		int id=it->first;
		fwrite(&id,sizeof(id),1,f);
		int no=it->second.size();
		fwrite(&no,sizeof(no),1,f);
		fwrite(& ((it->second)[0]) , sizeof((it->second)[0]),no,f);	
	}

}

void Objectdb::momPerRadius::load(FILE*f)
{
	uint32 buf;
	fread(&buf,sizeof(buf),1,f);
	if(buf==0)return;
	if(dbgout)(*dbgout)<<"loading "<<buf<<" moment sets:"<<endl;
	momentSets.resize(buf*tmom::n);
	fread(&momentSets[0],momentSets.size(),sizeof(momentSets[0]),f);
	if(dbgout)(*dbgout)<<"loading "<<buf<<" meta infos:"<<endl;
	meta.resize(buf);
	fread(&meta[0],sizeof(meta[0]),meta.size(),f);

	int nmultiObj;
	fread(&nmultiObj,sizeof(nmultiObj),1,f);
	for(int i=0;i<nmultiObj;++i)
	{
		int id;
		fread(&id,sizeof(id),1,f);
		tvpim &v=objectsnearby[id];
		int no;
		fread(&no,sizeof(no),1,f);
		v.resize(no);
		fread(&v[0],sizeof(v[0]),no,f);	
	}
}



void Objectdb::save(const char*name)
{
	FILE* f=fopen(name,"wb+");
	if(!f){if(dbgout)(*dbgout)<<"couldn't open"<<name<<endl;throw;}
	mapmom::iterator mit;
	
	uint32 buf = moments.size();
	fwrite(&buf,sizeof(buf),1,f);	
	for(mit=moments.begin();mit!=moments.end();++mit)
	{
		fwrite(&mit->first,sizeof(mit->first),1,f);	
		mit->second.save(f);	
	}

	buf=objects.size();
	fwrite(&buf,sizeof(buf),1,f);

	char charbuf[256];
	if(buf>0){		
		for(int i=0;i<objects.size();++i)
		{
			memset(charbuf,0,sizeof(charbuf));
			strncpy(charbuf,objects[i].name.c_str(),256);
			fwrite(charbuf,sizeof(charbuf),1,f);
		}
	}


	fclose(f);
}

void Objectdb::load(const char*name)
{
	FILE * f = fopen(name,"rb");
	if(!f){
		if(dbgout)(*dbgout)<<"couldn't open"<<name<<endl;
		throw std::runtime_error( (string("file")+name+"not found").c_str());
	}

	uint32 buf;
	mynum r;

	fread(&buf,sizeof(buf),1,f);
	for(uint32 i=0;i<buf;++i)
	{
		fread(&r,sizeof(r),1,f);
		moments[r].load(f);	
	}

	fread(&buf,sizeof(buf),1,f);
	if(dbgout)(*dbgout)<<"loading "<<buf<<" objects:"<<endl;	
	objects.clear();
	objects.resize(buf);

	char charbuf[256];
	for(int i=0;i<objects.size();++i)
	{
		memset(charbuf,0,sizeof(charbuf));
		fread(charbuf,sizeof(charbuf),1,f);
		objects[i].name=charbuf;
		if(dbgout)(*dbgout)<<"object"<<i<<": "<<charbuf<<endl;
	}

	fclose(f);

	/**calculate stats/object id lists*/	
	
	rmin=moments.begin()->first;
	rmax=moments.rbegin()->first;
	mapmom::iterator i;	
	vector<momsetMeta>::iterator it;
	for(i=moments.begin();i!=moments.end();++i)
	{
		vector<mynum>::iterator mit=i->second.momentSets.begin();			
		/*
		for(it=i->second.meta.begin();it!=i->second.meta.end();++it)
		{

			stats * ms= i->second.momstats, * msend=ms+tmom::n;
			for(;ms!=msend;++ms,++mit){
				ms->add(*mit);
			}	
		}
		*/

		for(;mit!=i->second.momentSets.end();mit+=tmom::n)
			i->second.vstats.add(&*mit);
	}

	updateObjects();
}


Objectdb::mynum Objectdb::momPerRadius::distance(const mynum*a,const mynum*b) const
{
	mynum ret=0;
	for(int i=0;i<treedim;++i)
	{
		const mynum* evi = convert[i];
		mynum c =  (a[0]-b[0]) * evi[0];
		for(int j=1;j<tmom::n;++j)
			c += (a[j]-b[j]) * evi[j];
		ret += c*c;
	}				

	return sqrt(ret);
}





void Objectdb::momPerRadius::tosearchable(mynum*out,const mynum*in) const
{
		for(int i=0;i<treedim&&i<tmom::n;++i)
		{
			out[i] = inner_product(in,in+tmom::n,convert[i],(mynum)0); 
		}				
}


void Objectdb::buildTrees()
{
	cerr<<"building all kdtrees..."<<flush;
	if(*dbgout)*dbgout<<"Building all kdtrees of db..."<<endl;
	mapmom::iterator it;
	for(it =	moments.begin();it!=moments.end();++it)
	{
		clock_t tstart=clock();
		if(*dbgout)*dbgout<<"Build kdtree for moments of radius: "<<it->first<<endl;
		it->second.getTree();
		if(*dbgout)*dbgout<<"ready. time used:"<<
			(clock()-tstart)*1000/CLOCKS_PER_SEC<<" ms"<<endl;
	}
	cerr<<"ready."<<endl;


	for(int i=0;i<tmom::n;++i){
		if(*dbgout)*dbgout
			<<"angles between principal dirs of "<<i<<"th moments :"<<endl;


		for(it =	++moments.begin();it!=moments.end();++it)
		{
			mynum c=0,asq=0,bsq=0;
			mynum minang=90;int minangid;

			const mynum* a=moments.begin()->second.convert[i];

			for(int k=0;k<tmom::n;++k)
			{
				const mynum* b=it->second.convert[k];
				for(int j=0;j<tmom::n;++j)
				{
					c+=a[j]*b[j];
					asq+=a[j]*a[j];		
					bsq+=b[j]*b[j];
				}
				mynum ang=acos(c/sqrt(asq*bsq)) *180/3.1415926;

				if(ang>90)ang=180-ang;

				if(ang<minang)minang=ang,minangid=k;
				//if(*dbgout)*dbgout<<ang<<"� "; 

			}
			if(*dbgout)*dbgout
				<<"   best match "<<minang<<" � at id"<<minangid<<endl;			
		}		
	}

}

void Objectdb::momPerRadius::initConvert(const vecstats<tmom::n,mynum>& vs)
{
	vs.initConvert(convert);
	treedim=tmom::n;
}

const Momkdtree &Objectdb::momPerRadius::getTree()
{
	if(momentSets.empty()||tree.getDim()!=0)
		return tree;




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
	vecstats<tmom::n,mynum> vs	= vstats.getstatscov();

	initConvert(vs);

	int nmoms=(momentSets.size()/tmom::n);
	vector<mynum> momsets2(treedim * nmoms ) ;	
	

	vecstatscollector<tmom::n,mynum> vc2;

	static const int maxtest=100;

	int numgt[tmom::n][2*maxtest+1];
	int numgtoi[tmom::n][2*maxtest+1];

	int maxmulti[maxtest+1];

	memset(numgt,0,sizeof(numgt));
	memset(numgtoi,0,sizeof(numgtoi));
	memset(maxmulti,0,sizeof(maxmulti));


	mynum avg2[tmom::n];
	tosearchable(avg2,vs.avg);

	//multiply moment vectors with matrix ev.
	vector<mynum>::iterator ai,oi;
	for(ai=momsets2.begin(),oi=momentSets.begin();
			oi!=momentSets.end();
			oi+=tmom::n , ai+=treedim )
	{
		//multiply with convert matrix
		tosearchable(&*ai,&*oi);
	
		//make statistics:not needed, it works now	
		vc2.add(&*ai);

		int iimin=maxtest;
		for(int i=0;i<tmom::n;++i)
		{
			int ii=int( floor(.5+ai[i]-avg2[i])  );		
			if(ii>maxtest)ii=maxtest;
			if(ii<-maxtest)ii=-maxtest;
			numgt[i][ii+maxtest]++;

			if(abs(ii)<iimin)iimin = abs(ii);
			
			ii= int( floor( .5+ (oi[i]-vs.avg[i])/vs.sigma[i]) );
			if(ii >  maxtest)ii=maxtest;
			if(ii < -maxtest)ii=-maxtest;
			numgtoi[i][ii+maxtest]++;

		}

		for(int i=0;i<iimin;++i)
		{
			maxmulti[i]++;
		}

	}

	vecstats<tmom::n,mynum> vvs = vc2;

	if(dbgout)
	{
		*dbgout<<"total no. of moment sets:"<<vs.num<<endl;

		*dbgout<<"for every number of std deviations,\n" 
			"the number of moment sets that are greater in all 28 dimenions"
			<<endl;

		for(int j=0;j<4 || maxmulti[j];++j)
		{
			*dbgout<<j<<":"<<maxmulti[j]<<endl;
		}

		for(int i=0;i<tmom::n;++i)
		{	

			*dbgout<<"original moment number"<<i<<":"<<endl;
			int no=0,noneg=0;
			for(int j=maxtest;j>=0 && no<100;--j)
			{
				if(numgtoi[i][maxtest+j])
				{
					no+=numgtoi[i][maxtest+j];
					*dbgout<<" # >="<<j<<":"<<no<<endl;					
				}

				if(numgtoi[i][maxtest-j])
				{
					noneg+=numgtoi[i][maxtest-j];
					*dbgout<<" # <=" << -j << ":"<<noneg<<endl;					
				}
			}
		}

		for(int i=0;i<tmom::n;++i)
		{
			*dbgout<<"avg"<<i<<":"<<avg2[i]<<" ==? " << vvs.avg[i] << endl;
			*dbgout<<"converted moment number"<<i<<":"<<endl;
			int no=0,noneg=0;
			for(int j=maxtest;j>=0 && no<100;--j)
			{
				if(numgt[i][j+maxtest])
				{
					no+=numgt[i][j+maxtest];
					*dbgout<<" # >="<<j<<":"<<no<<endl;					
				}

				if(numgt[i][maxtest-j])
				{
					noneg+=numgt[i][maxtest-j];
					*dbgout<<" # <= "<<-j<<":"<<noneg<<endl;					
				}
			}
		}
	}

	/*
	vector<int> momsleft;
	throwOutMoments2(momsets2,meta,0,0,momsleft);
	if(dbgout)(*dbgout)
		<< momsleft.size() * 100 / meta.size() 
		<<"% left after throwOutMoments2"
		<<endl;
	*/
	tree.init(momsets2,treedim);
	if(dbgout)(*dbgout)<<"tree has"<<treedim<<"dimensions and a height of"<<tree.getHeight()<<endl;



	/*
	//compute new(global) scaling factors
	for(int j=0;j<tmom::n;++j)
	{	
		mynum v =momstats[j].variance();
		if(v<1e-6)
			scale[j]=1e3;
		else
			scale[j]= 1.0/sqrt(v); 		
	}

	vector<mynum>::iterator mit;		

	//scale values by sc
	for(mit=momentSets.begin();mit!=momentSets.end();)
	{
		for(int j=0;j<tmom::n;++j)
			*(mit++) *= scale[j];
	}
	//insert them into tree
	tree.init(momentSets,tmom::n);
	*/


	return tree;
}

std::vector<Objectdb::mynum> & Objectdb::momPerRadius::getMomSets()
{
/*
	if(!momentSets.empty())
		return momentSets;

	tree.swapleaves(momentSets);

	mynum sc[tmom::n];
	
	//compute inverse of scalin factors
	for(int j=0;j<tmom::n;++j)
	{	
		sc[j]= 1.0/scale[j]; 		
	}

	vector<mynum>::iterator mit;		

	//restore values in original scale
	for(mit=momentSets.begin();mit!=momentSets.end();)
	{
		for(int j=0;j<tmom::n;++j)
			*(mit++) *= sc[j];
	}
*/
	return momentSets;
}


/*add an object to the database
*( it is assumed that its in the same scale/coordinates
*  as the other objects )
*\param name: name of the object. if an object with the 
* same name already exists	
*, add the moments of this pointcloud to this Object.
*\param points: the pointcloud describing that object	
*\param clusterIt: should we throw out moments too near to each other ?
*\param cleanIt: should we "clean" the database by throwin out moments too near to other ones ?
*/
const Objectdb::obj& Objectdb::addObject
	(  const std::string&name,
	const std::vector<mynum> & points,
	bool clusterit
	)
{


	if(dbgout)(*dbgout)<<
		"---------------------------------------\n"
		"adding object "<<name<<"\n"
		"---------------------------------------"<<endl;


	int objectid=0;
	while(objectid<objects.size() && objects[objectid].name!=name)
		++objectid;

	bool objectIsNew = (objectid>=objects.size());

	if(objectIsNew)
		objects.resize(objects.size()+1);

	obj & ob =objects[objectid];
	ob.name = name;


	TsTree t(points);
	if(dbgout)(*dbgout)<<"built tree."<<endl;

	int raex=0,raexmin=min(t.getDepth(),max(1,t.getDepth()-3)),
				raexmax=t.getDepth();

	mynum ra = ldexpl( t.getScale(),raexmin),
		rama=ldexpl( t.getScale(),raexmax)
		;
	if(objects.size()==1)
	{
		rmin=ra; rmax=rama;
	}
	else{
		if(rmin>ra)rmin=ra;
		if(rmax<rama)rmax=rama;
	}

	vector<mynum> m,centers;

	for(int raex=raexmin;raex<raexmax;++raex)
	{
		mynum r=ldexpl(t.getScale(),raex);
		obj::stlist & s = ob.momstats[r];
		momPerRadius &mpr=moments[r];

		int nmomsets = getMoms(m,centers,s,t,r, max(-1,raex-3) );

		if(dbgout)(*dbgout)<<"computed moments for radius "<<r<<endl;

		vector<int> added;


		vecstats<tmom::n,mynum> vs;
		vs = s; 

		if(clusterit)
		{

			mynum sc[tmom::n],sub[tmom::n];
			for(int j=0;j<tmom::n;++j)
			{	
				//compute (local) scaling factors
				mynum v = vs.sigma[j];
				//if(dbgout)(*dbgout)<<"var:"<<v<<endl;
				sc[j]= v<1e-6 ? 1e6 : 1.0/v; 
				sub[j]=vs.avg[j];
			}
			clusterMoms(m,sub,sc,added);
		}
		else //don't cluster
		{
			added.resize(centers.size()/3);
			for(int i=0;i<added.size();++i)
				added[i]=i;
		}

		//join statistics		
		mpr.vstats.add(s);		
		mpr.getMomSets();
		int omsiz=mpr.momentSets.size();
		mpr.momentSets.resize(omsiz+ tmom::n * added.size());
		_TYPENAME_ vector<mynum>::iterator mit;
		mit=mpr.momentSets.begin()+omsiz;	

		//add m at the end of the old array momentSets,
		//init momentSets of objects
		vector<int> &obmomsets=ob.momentSets[r];

		obmomsets.reserve(added.size());
		int ometasiz = mpr.meta.size();
		mpr.meta.resize(mpr.meta.size()+added.size());

		vector<momsetMeta>::iterator metaIt = mpr.meta.begin()+ometasiz;
		for(int i=0;i<added.size();++i,++metaIt)
		{
			int ai=added[i];
			obmomsets.push_back(mit-mpr.momentSets.begin());
			vector<mynum>::iterator it = m.begin()+ai * tmom::n;
			mit = copy(it,it+tmom::n,mit);

			it=centers.begin()+ ai * 3;

			metaIt->center[0]=it[0];
			metaIt->center[1]=it[1];
			metaIt->center[2]=it[2];

			metaIt->radius = r;
			metaIt->objectid = objectid;

		}		


		if(dbgout)(*dbgout)<<"--------------------------------\n"<<endl;
	}
	return ob;
}





namespace Objectdb_impl{namespace for_reduce{


	size_t hash_value(const Objectdb_impl::ivec& v)
	{
			size_t ret=0;
			for(int i=0;i < v.size();++i)
			{
				ret*=0x522b2a7b;ret+=v[i];
			}
			return ret;
	};


	struct mymapper
	{
		mynum delt,scale;

		Objectdb::momPerRadius &om;

		mutable ivec buf;

		mutable mynum mm[tmom::n];

		mymapper(Objectdb::momPerRadius& om_):om(om_){
			delt=0;
		}

		ivec operator() (int id) const
		{
			//transform to unit covariance
			om.tosearchable(mm, &om.momentSets[id * tmom::n] );

			for(int j=0;j<tmom::n;++j)
			{	

				//if(dbgout)(*dbgout)<<v<<' '<<flush;
				buf[j]=(int)floor(mm[j]*scale+delt);			
				//if(dbgout)(*dbgout)<<buf[j]<<' '<<flush;
			}
			return buf;
		}
	};


	/**
	* create downwardly - linked lists for every group, given an upwardly-linked list for ids of the representatives.
	* the last item in the downwardly linked list is always a representative.
	* mapper maps an int to the destination type Dest
	* Dest has to be comparable via < .
	*/
	template<class Dest,class mapper>
	void joinEquallyMapped(vector<int> & added, const mapper& m)
	{
	
		typedef _HASH_MAP_<Dest,int> mymap;

		mymap momi;

		int prev=0,id=0;

		while(id<added.size())
		{
			int next=added[id];

			pair< _TYPENAME_ mymap::iterator,bool> ret=momi.insert(_TYPENAME_ mymap::value_type(m(id),id));
			if(!ret.second)
			{
				//remove repr.from upwardly linked list
				added[prev]=next;
				//set representative for this moment
				added[id]=ret.first->second;
			}
			else
			{
				prev=id;
			}

			id=next;
		}
	
	}
	
	/** extract groups from the lists in added. destroys added[]*/
	void getGroups(vector<int>&offs,vector<int>&groups,vector<int>&added)
	{
		//count num. repr, replace repr. ids by their position in the off array.
		int numleft=0,id=0;
		while(id<added.size())
		{
			int aid=added[id];
			added[id] = -numleft-1;
			id=aid;
			++numleft;
		}



		offs.resize(numleft+1);
		groups.resize(added.size());
	
		//calc. group sizes by adding to off.
		for(int i=0;i<added.size();++i)
		{
			int offid=added[i];

			while(offid>=0) //if not a representative:
			{
				offid=added[offid];
			}
			offid=-offid-1;

			offs[offid]++;
		}

		for(int i=1;i<offs.size();++i)
		{
			offs[i]+=offs[i-1];
		}
		offs.back() = added.size();

		//insert into groups array:(s.t. repr is first, order is reversed)
		for(int i=added.size()-1;i>=0;--i)
		{
			int offid=added[i];

			while(offid>=0) //if not a representative:
			{
				offid=added[offid];
			}
			offid=-offid-1;

			groups[--offs[offid]]=i;
		}

		added.clear();
	}


	inline void test_JoinEquallyMapped()
	{
	
		int a[9]={1,2,3,4,5,6,7,8,9};

		vector<int> added(a,a+9);

		struct mapper
		{

			int operator()(int i) const
			{
				const int map[] = {5,2,2,7,7,7,2,2,3};
				return map[i];
			}
		};

		joinEquallyMapped<int,mapper>(added,mapper());

	
		vector<int> off;
		vector<int> groups;

		getGroups(off,groups,added);
	
	}


}}

//hash_value must be put in this namespace to be found by instanciation of hash_compare
namespace stdext
{
	size_t hash_value(const Objectdb_impl::ivec& v)
	{
		return  Objectdb_impl::for_reduce::hash_value(v);
	}
};


void Objectdb::momPerRadius::reduce
	(Objectdb::mynum clusterSize)
{
	OUTPUTFLOW();

	using namespace Objectdb_impl::for_reduce;


	initConvert(vstats.getstatscov());

	/** map an int vector to the id of a representative moment*/

	typedef map< ivec ,int > mymap;

	vector<mynum>& m=getMomSets();

	int nmomsets=m.size()/tmom::n;

	ivec buf;

	mynum mm[tmom::n];



	/** the remaining moment ids. in the first run, used for liked lists of equal moments.*/
	vector<int> added(nmomsets);

	//the representatives, if i>added[i], a represented if i<added[i]
	//is alinked list, always if i> added i.
	for(int i=0;i<added.size();++i)
		added[i]=i+1;

	if(dbgout)(*dbgout)<<added.size()<<" moment sets before."<<endl;

	vector<int>::const_iterator ait;
	//throw out moments in 29 steps.
	//scale up ,, so a distance of 0.2 variance still matters.
	

	mymapper mapper(*this);
	
	mapper.scale=1.0/clusterSize;
	for(int i=0;i<=tmom::n;++i)
	{
		mapper.delt=0.5+ i*1.0/(tmom::n+1);
		joinEquallyMapped<ivec,mymapper>(added,mapper);
	}


	vector<int> offs,groups;
	getGroups(offs,groups,added);

	tvpim objbuf;
	
	tobjectsnearby nhm;

	vector<int>::iterator offstart=offs.begin(),offit;
	for(offit=offstart+1;offit!=offs.end();offstart=offit,++offit)
	{

		int remainid = groups[*offstart];
		added.push_back(remainid);

		if( *offit - *offstart > 1)
		{

			int baseoid=meta[remainid].objectid;

			tobjectsnearby::iterator ridit = objectsnearby.find(remainid);
			if(ridit!=objectsnearby.end())
				objbuf.insert(objbuf.end(),ridit->second.begin(),ridit->second.end());

			for(int momids= *offstart + 1;momids<*offit;++momids)
			{
				int gm = groups[momids];
				Objectdb::momsetMeta &met= meta[gm];

				mynum dist=distance( &m[0]+ remainid * tmom::n, &m[0]+ gm * tmom::n );

				if(met.objectid!=baseoid)
					objbuf.push_back(pair<int,mynum>(met.objectid,dist));

				tobjectsnearby::iterator it = objectsnearby.find( gm);
				if(it!=objectsnearby.end())	{
					for(vector<pair<int,mynum> > ::iterator vit=it->second.begin();vit!=it->second.end();++vit)
					{
						if(vit->first!=baseoid)
							objbuf.push_back(pair<int,mynum>(vit->first, vit->second + dist));
					}
				}
			}

			sort(objbuf.begin(),objbuf.end());

			struct equal1st
			{
				bool operator()(const pair<int,mynum>&a,const pair<int,mynum> & b) const
				{
					return  a.first==b.first;
				}
			};

			//only keep lowest distances per object
			tvpim::iterator eit = unique(objbuf.begin(),objbuf.end(),equal1st());
			if(eit!=objbuf.begin())
			{
				tvpim &newv = nhm[added.size()-1];
				newv.insert(newv.begin(),objbuf.begin(),eit);
			}
			objbuf.clear();
		}
	}

	//throw out those which are not in added
	//from momentSets, objectofmoment
	vector<mynum>::iterator mit=m.begin(),mitbuf;
	vector<momsetMeta>::iterator oit=meta.begin();


	for(int i=0;i<added.size();++i,mit+=tmom::n,++oit)
	{
		int ai=added[i];

		if(ai==i)
			continue;

		*oit = meta[ai];

		mitbuf=m.begin()+ ai *tmom::n;
		copy(mitbuf,mitbuf+tmom::n,mit);
	}

	objectsnearby.swap(nhm);

	if(dbgout)(*dbgout)<<added.size()<<" moment sets afterwards."<<endl;

	if(dbgout)*dbgout <<"ratio of moments with multiple objects:"<< 1.* objectsnearby.size()/added.size()<<endl;

	m.resize(added.size()*tmom::n);
	meta.resize(added.size());		

}



void Objectdb::updateObjects()
{
	OUTPUTFLOW();
	mapmom::iterator i;	
	vector<momsetMeta>::iterator it;

	for(int i=0;i<objects.size();++i)
		objects[i].momstats.clear(),objects[i].momentSets.clear();

	for(i=moments.begin();i!=moments.end();++i)
	{
		vector<mynum>::iterator mit=i->second.momentSets.begin();			
		for(it=i->second.meta.begin();it!=i->second.meta.end();++it)
		{
			Objectdb::obj & o=objects[it->objectid];

			vector<int>&momids = o.momentSets[i->first];
			momids.push_back(mit - i->second.momentSets.begin());

			o.momstats[i->first].add(&*mit);
		}
	}

}


void Objectdb::reduce
	(Objectdb::mynum clusterSize)
{
	OUTPUTFLOW();
	for(mapmom::iterator i=moments.begin();i!=moments.end();++i){
		if(dbgout)*dbgout<<"reduce r="<<i->first<<":"<<endl;
		i->second.reduce(clusterSize);		
	}
	updateObjects();		
}

int Objectdb::nmoments()
{
				int ndbsiz=0;
			for(Objectdb::mapmom::const_iterator it=getMomsets().begin();
				it!=getMomsets().end();++it)
			{
				ndbsiz+= it->second.momentSets.size()/tmom::n;
			}
			return ndbsiz;
}




namespace Objectdb_impl{

	typedef Objectdb::match match;
	typedef Objectdb::mapmom mapmom;
	typedef Objectdb::momPerRadius momPerRadius;
	

	/**
	get and transform invariant set to be searched in momPerRadius
	*/ 
	int getInv( double ra,	const vector<mynum> & points, const TsTree&t, mynum &pnpr,
		mynum*m)
	{
		tensorSet3D<> ts,buf;
		tmom bufm;
	

		//get ball around arbitrarily chosen point			
		int ptid= int(double(points.size()/3)*rand()/RAND_MAX);			
		if(ptid>= points.size()/3)
			ptid=points.size()/3 -1;
		t.sumInCircle(ts,buf,
			(mynum*)&points[ptid*3],
			ra,ra<=t.getScale()*4 ? -1 : 0);
		pnpr +=ts.A[0];
		ts.transinv();
		Objectdb::kdtreeEntryFromTset(m,bufm,ts,ra);
		return ptid;
	}

	/**
		process found moment ids together with current invariant set
	*/
	void processMomids(
		const vector<int>&momids,
		vector<mynum>&dpo,
		momPerRadius& itsecond,
		
		const mynum*m,
		const vector<int>*classes,
		bool useEuclideanMetric,
		bool matches,
		vector<int>&bestmomperobject
		)
	{
		const Momkdtree & tr=itsecond.getTree();


		fill(dpo.begin(),dpo.end(),1e+300);



		for(int i=0;i<momids.size();++i)
		{
			//d is now the max. 
			//distance between the current moment 
			//and the moment of an object
			mynum d;
			if(!useEuclideanMetric)
				d= maxd(
				& (tr.getleaves()[momids[i]]),
				m , tr.getDim() );
			else
				d= eucd(
				& (tr.getleaves()[momids[i]]),
				m , tr.getDim() );

			int momid = momids[i]/tr.getDim();
			int objectid = itsecond.meta[momid].objectid;




			if(classes)
				objectid=(*classes)[objectid];

			mynum &dp =dpo[objectid];
			if(dp>d){
				dp=d;
				if(matches)
					bestmomperobject[objectid] = momid;
			}

			momPerRadius::tobjectsnearby::const_iterator itonb= itsecond.objectsnearby.find(momid);
			if(itonb!=itsecond.objectsnearby.end())
			{

				for(momPerRadius::tvpim::const_iterator vp=itonb->second.begin();vp!=itonb->second.end();++vp)
				{
					objectid=vp->first;
					mynum dd=d+vp->second;
					if(classes)
						objectid=(*classes)[objectid];

					mynum &dp =dpo[objectid];
					if(dp>dd){
						dp=dd;
						if(matches)
							bestmomperobject[objectid] = momid;
					}
				}
			}
		}

	}











	void getOids(vector<int>&oids,
		momPerRadius&itsecond,
		vector<int>&momids,
		mynum&maxdist,
		bool useEuclideanMetric,
		mynum*box,
		const mynum*m,
		vector<mynum>&dpo,
		int&numids,
		vector<int>&timespresent,
		vector<mynum>&timespresentfuzzy,
		vector<int>&bestmomperobject,
		vector<match> * matches,
		const vector<int>*classes
		)
	{



		oids.clear();


		const Momkdtree & tr=itsecond.getTree();

		const mynum maxdd=4;

		momids.clear();
		for( maxdist=ldexpl(1,-4);	maxdist<=maxdd &&oids.size()<5; maxdist*=2)
		{

			if(!useEuclideanMetric)
			{
				mynum*bxit=box;
				for(int j=0;j<tr.getDim();++j,bxit+=2)
				{ 
					bxit[0]=m[j]-maxdist;
					bxit[1]=m[j]+maxdist;
				}
				tr.momIdsInBox(momids,box);		
			}
			else
				tr.momIdsInBall(momids,m,maxdist);

			/*
			if(maxdist>=maxdd)
			{
				//use all moment ids!
				momids.resize(itsecond.meta.size());
				for(int i=0;i<momids.size();++i)
					momids[i]=i*tr.getDim();
			}
			else
			*/
			if(momids.empty()) continue;
			


			//processMom
			processMomids(
				momids,
				dpo,
				itsecond,

				m,
				classes,
				useEuclideanMetric,
				matches,
				bestmomperobject
				);


			numids += momids.size();


			oids.clear();

			//collect object ids which are near
			for(int i=0;i<dpo.size();++i)
				if(dpo[i]<maxdd)
					oids.push_back(i),++timespresent[i],timespresentfuzzy[i]+=1/(1+dpo[i]);

			//	allmomids.insert(momids.begin(),momids.end());
		}






	}

	void processOids2(
		vector<int>&oids,
		const vector<mynum>&dpo,
		vector<mynum>&ret,
		mynum&minprob,
		vector<int>&bestmomperobject,
		const vector<mynum> & points,
		int ptid,
		vector<match> * matches,
		double ra,
		const Objectdb::momPerRadius & itsecond
		)
	{
		if(oids.size()==0)
			return;

		mynum mindpo=dpo[0];
		int mindpoid=0;
		for(int i=0;i<oids.size();++i)
		{
			mynum d=dpo[oids[i]];
			if(d<mindpo)
			{
				mindpo=d;
				mindpoid=i;
			}
		}
		
		mynum isigma= 1.0/(mindpo+.1);

		mynum sump=0;
		for(int i=0;i<oids.size();++i)
		{
			mynum d = dpo[oids[i]]*isigma;
			sump+=bellcurve(d);
		}

		//big mindpo values should contribute less.
		sump/=bellcurve(mindpo);


		for(int i=0;i<oids.size();++i)
		{
			int oi=oids[i];
			mynum d = dpo[oi]*isigma;
			ret[oi] *= (1 - bellcurve(d)/sump );
			if(minprob>ret[oi])minprob=ret[oi];
		}

		if(matches)
		{
			matches->push_back(match());
			match&m=matches->back();
			m.dbmom = itsecond.meta[bestmomperobject[oids[mindpoid]]];

			copy( &points[ptid*3],&points[ptid*3]+3,m.testmom.center);
			m.testmom.objectid=oids[mindpoid];
			m.testmom.radius=ra;

			m.dist=dpo[mindpoid];

			m.seconddist=1e300;

			for(int i=0;i<oids.size();++i)
			{
				if(i==mindpoid)
					continue;

				mynum d=dpo[oids[i]];
				if(d<m.seconddist)
				{
					m.seconddist=d;
				}
			}
		}	
	
	}







	/** process vector of object ids which were found*/
	void processOids(
		vector<int>&oids,
		const vector<mynum>&dpo,
		vector<mynum>&ret,
		mynum&minprob,
		vector<int>&bestmomperobject,
		const vector<mynum> & points,
		int ptid,
		vector<match> * matches,
		double ra,
		const Objectdb::momPerRadius & itsecond
		){
			if(oids.size()>1){

				//sort object ids by their dpo values
				sort(oids.begin(),oids.end(),cmpindices(dpo));


				//look for a gap in the distance values
				int gapid=0;

				while(gapid+1<oids.size() && 
					dpo[oids[gapid]]*1.2 + 0.1 > dpo[oids[gapid+1]] )
					++gapid;


				mynum mult=1;
				if(gapid+1 == oids.size())
				{
					if(gapid==0)mult=.5;
					else mult=.99;

					if(dbgout)(*dbgout)<<"gapid "<<gapid<<';'<<flush;
				}
				else
				{
					mult= 1 - ( (dpo[oids[gapid+1]] -dpo[oids[gapid]] ));
					if(mult<0)mult=0;
					
					if(dbgout)(*dbgout)<<"gap r="<<ra<<':'<<dpo[oids[gapid+1]]<<' '<<dpo[oids[gapid]]<<';'<<flush;

				}

				//in case no gap was found:
				if(dpo[oids[gapid]]<0.5 && /*gapid<=objects.size()/2*/gapid==0 )
				{
					if(dbgout)(*dbgout)<<'<'<<ra<<':'<<oids[0]<<'>'<<flush;
					for(int i=0;i<=gapid;++i)
					{
						ret[oids[i]] *= mult;

						if(matches)
						{
							matches->push_back(match());
							match&m=matches->back();
							m.dbmom = itsecond.meta[bestmomperobject[oids[i]]];

							copy( &points[ptid*3],&points[ptid*3]+3,m.testmom.center);
							m.testmom.objectid=-1;
							m.testmom.radius=ra;

							m.dist=dpo[oids[i]];
							m.seconddist = i+1 < oids.size() ? dpo[oids[i+1]] :1 ;
						}

						if(minprob>ret[oids[i]])minprob=ret[oids[i]];
					}

				}

			}

			/*
			if(dbgout)(*dbgout)<<"# momids:"<<momids.size()<<" oids:"<<oids.size()<<" maxdist:"<<maxdist<<" radius:"<<ra<<endl;
			*/
			if(oids.size()==1)
			{
				if(dbgout)(*dbgout)<<'('<<ra<<':'<<oids[0]<<')'<<flush;

				if(matches)
				{
					matches->push_back(match());
					match&m=matches->back();
					m.dbmom = itsecond.meta[bestmomperobject[oids[0]]];

					copy( &points[ptid*3],&points[ptid*3]+3,m.testmom.center);
					m.testmom.objectid=-1;
					m.testmom.radius=ra;

					m.dist=dpo[oids[0]];
					m.seconddist = 1 ;
				}
				ret[oids[0]]*=dpo[oids[0]];
				if(minprob>ret[oids[0]])minprob=ret[oids[0]];
			}

	}


};




void Objectdb::findObjects(
		//the probabilities
		std::vector<mynum> &ret,
		const std::vector<mynum> & points
	
		,std::vector<match> * matches
		,const std::vector<int> * classes
		,bool ignoreRadii
	)
	{

		OUTPUTFLOW();

		//for evry radius, sum up no. of points in summing
		vector<mynum> nperRadius(moments.size(),0);

		bool useEuclideanMetric=true;

		if(dbgout)(*dbgout)<<(useEuclideanMetric? " euklidean":"max")<<"-metric"<<flush;
		clock_t t0=clock();
		TsTree t(points);
		if(dbgout)
		{
			*dbgout<<"time used for building tree:"
				<<1000*(clock()-t0)/CLOCKS_PER_SEC<<" ms "<<endl;
		}


		mynum box[ tmom::n * 2 ];

		//moments from kdTreeEntry,moments converted to the tree format
		mynum m[tmom::n],om[tmom::n];

		//maximum distance of moments we search for
		mynum maxdist;
		ret.resize(classes? classes->size():objects.size());

		//ret is now the probability that an object is not inside
		fill(ret.begin(),ret.end(),1);

		//buffer for distance per object
		vector<mynum> dpo(classes? classes->size():objects.size());
		//buffer for momids
		vector<int> momids;

		//buffer for object ids near an actual moment
		vector<int> oids;

		tmom bufm;

		tensorSet3D<> ts,buf;

		mynum minprob=1;


		//the number
		vector<int> timespresent(objects.size(),0);
		vector<mynum> timespresentfuzzy(objects.size(),0);

		int numids=0,forgetids=0;
		int numtsets=0;
		//if(dbgout)(*dbgout)<<"testing distance "<< maxdist<<'\r'<<flush;

		if(dbgout){
			*dbgout <<"radii under consideration:";

			for( mapmom::iterator it = moments.begin();	it!=moments.end();++it)
				if(it->first>t.getScale()*2)
					*dbgout <<' '<<it->first;
			*dbgout<<endl;
		}


		for(int ii=1;
			(ii<10|| ( minprob> (0.01))
			&&ii<1000);
		++ii)
		{		
			int ptid;	
			mynum ra;


			if(ignoreRadii)
			{
				ra=ldexp(t.getScale(),t.getDepth() *rand()/RAND_MAX );
				mynum num;
				ptid=getInv(ra,points,t,num,om);
			}


			//if(dbgout)(*dbgout)<<"testing radii "<< rmin <<" - "<<rmax<<'\r'<<flush;
			vector<mynum>::iterator npr=nperRadius.begin();
			for( mapmom::iterator it = moments.begin();	it!=moments.end();++it,++npr)
			{
				if(!ignoreRadii)
				{
					ra=it->first;
					if(ra<=t.getScale()*2) continue;

					if(ra>ldexp(t.getScale(),t.getDepth())) break;
					ptid=getInv(it->first,points,t,*npr,om);
				}
				//convert to scale for current inv.tree
				it->second.tosearchable(m,om);
				++numtsets;

				//for stats/vis:
				vector<int> bestmomperobject;
				if(matches)
					bestmomperobject.resize(objects.size());


				getOids(oids,
					it->second,
					momids,
					maxdist,
					useEuclideanMetric,
					box,
					m,
					dpo,
					numids,
					timespresent,
					timespresentfuzzy,
					bestmomperobject,
					matches,
					classes
					);



				processOids2(
					oids,
					dpo,
					ret,
					minprob,
					bestmomperobject,
					points,
					ptid,
					matches,
					ra,
					it->second);


			}



		}
		//			if(dbgout)(*dbgout)<<"                                              \r";

		for(int i=0;i<ret.size();++i)
			ret[i]=1-ret[i];

		

		if(dbgout)(*dbgout)<<"num. of moment ids found during getMomIdsInBox:"<<numids<<endl;

		if(dbgout)(*dbgout)
			<<"num. of tensorsets computed:"<<numtsets<<endl;

		if(dbgout)
		{
			(*dbgout)<<"presence per object:";
			vector<pair<mynum,int> > s;
			for(int i=0;i<timespresent.size();++i)
				if(timespresent[i])
					s.push_back(pair<mynum,int>(timespresent[i],i));

			sort(s.begin(),s.end());

			for(int i=0;i<s.size();++i)
					(*dbgout)<<s[i].second<<':'<<s[i].first<<' ';
			(*dbgout)<<endl;
			
			s.clear();
			for(int i=0;i<timespresentfuzzy.size();++i)
				if(timespresentfuzzy[i])
					s.push_back(pair<mynum,int>(timespresentfuzzy[i],i));

			sort(s.begin(),s.end());

			for(int i=0;i<s.size();++i)
					(*dbgout)<<s[i].second<<':'<<s[i].first<<' ';
			(*dbgout)<<endl;
		}
		






}

int Objectdb::getInv( double ra,	const vector<mynum> & points, const TsTree&t, mynum &pnpr,
		mynum*m)
{
	return Objectdb_impl::getInv(ra,points,t,pnpr,m);
}



int Objectdb::getMoms(
		//the invariant moments
		std::vector<mynum>& m,
		//the corresponding points
		std::vector<mynum>& pt,
		//statistik
		vecstatscollector<tmom::n,mynum> &s,
		//Baum, der die tensorsets enthaelt				
		const TsTree&t,
		//radius of the balls
		mynum ra,
		//minimum height to go down in the tree
		int minh
	)
{
	return Objectdb_impl::getMoms(m,pt,s,t,ra,minh);
}

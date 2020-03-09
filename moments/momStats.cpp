/*
	collections of tests /stats output using the moment invarints for object recognition

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

#include<time.h>
#include<cassert>
#include<fstream>
#include<sstream>
#include<vector>
#include<map>
#include<algorithm>
#define _USE_MATH_DEFINES
#include<cmath>
#include<stdlib.h>
#include<sys/stat.h>
#include<sys/types.h>
#include<fcntl.h>

#ifdef _DEBUG
#define TNT_BOUNDS_CHECK
#endif

#include<tnt/jama_eig.h>


#ifndef GNUCC

#include<io.h>
#include"dirent.h"
#include<direct.h>


#else

#include<dirent.h>
#include<errno.h>
#include<unistd.h>
#include<sys/times.h>

#endif

#include "pointCloud.h"
#include "tensorSet.h"
#include "TsTree.h"
#include "TsTree2.h"
#include "Objectdb.h"
#include "clustertree.h"
#include "vecstats.h"
#include "PointCloudIO.h"
//#include "mls.h"
#include "intervalnewton.h"
#include "Momkdtree_templ.h"
#include "cluster.h"
#include "infbruch.h"

#include "momStats.h"
#include "osgio.h"

#ifdef USE_PSB
#include "PSBCategories.h"
#endif

#ifdef USE_SVM
#include "Objectdbsvm.h"
#endif

#ifndef GNUCC
float myclock()
{
	return float(clock())/CLOCKS_PER_SEC;
}
#else

clock_t ticks_per_sec=sysconf(_SC_CLK_TCK);

float myclock()
{
	struct tms t;
	times(&t);
	return float(t.tms_utime)/ticks_per_sec;
}
#endif

using namespace std;



void testreadGAVABvrmlpts()
{

	vector<mytype> pts;

	PointCloudIO::read(pts,"C:\\u\\max\\models\\GAVAB\\allfaces\\cara1_abajo.wrl");

	cout<<pts.size()%3<<' '<<pts.size()/3<< "pts read"<<endl;

}


ostream& operator<<(ostream&out,Objectdb::stats&s)
{
	out<<s.average()<<' '<<sqrt(s.variance())<<' '<<s.min<<' '<<s.max;
	return out;
}









inline void covps(const mynum cov[28][28],ostream&out,int rh)
{

	for(int i=0;i<28;++i)for(int j=0;j<28;++j)
	{
		out<<rh*i<<' '<<rh*j<<" moveto"<<endl;
		out<<rh<<" 0 rlineto  0 "<<rh<<" rlineto "
				<< -rh << " 0 rlineto"<<endl;
		out<<1.0-fabs(cov[i][j])<<" setgray fill"<<endl;							
	}

	int h=rh*28;

	out<<"0 0 moveto "<<h<<" 0 rlineto  0 "<<h<<" rlineto "
			<< -h << " 0 rlineto closepath 0 setgray stroke"<<endl;
}


inline string nameps(string s)
{

	int lasti=0;
	for(int i=0;i<s.size();++i)
	{
		if(s[i]=='\\')
			s[i]='/';
		if(s[i]==')')
			s[i]=']';
		if(s[i]=='(')
			s[i]='[';
	}
	return s;
}













struct momStats2{

	typedef momentSet<mytype,28,0> tmom;
	static const int n=28;

	//calculate variance , gravitycenter, bbox of dataset

	vecstats<28> 
	stats, //direct statistis of moments
	stats2,//statistics of moments rooted by domain order
	stats3;//statistics of moments rooted by value order

	string name;

	/**
	 * create a stats from dump file fname
	 *\param valdomsi:
	 * should we norm the moments in value+domain ?
	 */
	momStats2(){};

	momStats2(const vector<mytype> &moms,bool valdomsi,string nam)
	{name=nam;init(moms,valdomsi);}

	momStats2(const char*fname,bool valdomsi=false)
	{
		name=fname;
		vector<mynum> moms;
		cout<<"reading in data from "<<fname<<"..."<<endl;
		PointCloudIO::read(moms,fname,".raw");

		init(moms,valdomsi);
	}

	void init(const vector<mytype> &moms,bool valdomsi)
	{
		cout
		<<(valdomsi?"normalize,":"")
		<<"calculating val. order root, "
		"dom.order root, min,max,avg,var ..."<<flush;

		vector<mytype> moms2(moms.size()),moms3(moms.size());

		vecstatscollector<n> st,st2,st3;

		for(int i=0;i<moms.size();i+=n)
		{		
			tmom &m=*((tmom*)&moms[i]);
			tmom &m2=*((tmom*)&moms2[i]);
			tmom &m3=*((tmom*)&moms3[i]);

			//scale moms
			if(valdomsi){
				double vscale=sqrt(m.M[0]/m.M[21]);
				m.scaleDomain(m.M[21]/(m.M[0]*m.M[0]));
				m.scaleValue(vscale);
			}

			for(int j=0;j<28;++j)
			{		
				mytype pwdord=1.0/tmom::metainfo.me[j].dord; 
				mytype pwvord=1.0/tmom::metainfo.me[j].vord; 
				m2.M[j]=(m.M[j]>0?1:-1)*pow(fabs(m.M[j]),pwvord);
				m3.M[j]=(m.M[j]>0?1:-1)*pow(fabs(m.M[j]),pwdord);			
			}

			st.add(m.M);st2.add(m2.M);st3.add(m3.M);
		}

		stats  = st;
		stats2 = st2;
		stats3 = st3;
	}



	void write(ostream&o) const
	{
		o<<"\nF\n"<<(name.size()>0?name:"unknown")
					<<"\nS\n"<<stats
					<<"\nT\n"<<stats2
					<<"\nU\n"<<stats3;
	}

	bool read(istream&in)
	{
		char c;
		in>>c;
		if(in.eof()||in.bad()||c!='F')
		{in.setstate(ios::badbit|ios::failbit);return false;}
		in>>name;
		in>>c;
		if(c!='S')return false;
		in>>stats;
		in>>c;
		if(c!='T')return false;
		in>>stats2;
		in>>c;
		if(c!='U')return false;
		in>>stats3;
		return true;
	}


	void writecorrps(ostream&out)
	{
		typedef momStats2::tmom tmom;

		int a4w = int(21 *72 / 2.54);
		int a4h = int(29 *72 / 2.54);



		cout<<"writing to ps file..."<<flush;

		out<<"/Times-Roman findfont 20 scalefont setfont"<<endl;
		out<< 30 << ' '<<a4h -60<<" moveto ("<<nameps(name)<<") show"<<endl;

		out<<"gsave 30 30 translate"<<endl;
		covps(stats.cov,out,a4h/4/28);
		out<<a4h/4<<" "<<a4h/8<<" moveto (correlation matrix) show"<<endl;
		out<<"grestore"<<endl;

		out<<"gsave 30 "<<30+a4h/4<<" translate"<<endl;
		out<<a4h/4<<" "<<a4h/8<<" moveto (correlation of  pow\\(x,1/value_order\\)) show"<<endl;
		covps(stats2.cov,out,a4h/4/28);
		out<<"grestore"<<endl;

		out<<"gsave 30 "<<30+a4h/2<<" translate"<<endl;

		out<<a4h/4<<" "<<a4h/8<<" moveto (correlation of pow\\(x,1/domain_order\\)) show"<<endl;
		covps(stats3.cov,out,a4h/4/28);
		out<<"grestore"<<endl;
		out<<"showpage"<<endl;


		float covf[28*28];
		string covn[28];
		for(int i=0;i<28;++i){
			string s = tmom::metainfo.me[i].des;

			for(int j=0;j<s.size()-1;++j)
			{
				if(s[j]==')')s[j]=']';
				if(s[j]=='(')s[j]='[';
				if(s[j]=='\\')s[j]=' ';
			}
			covn[i]=s;
		}


		std::copy(&stats.cov[0][0],&stats.cov[0][0]+28*28,covf);
		out<<"/Times-Roman findfont 3 scalefont setfont"<<endl;
		out<<"20 20 translate 3 4 scale"<<endl;
		out<<"100 110 moveto (direct moments) show"<<endl;
		drawCorrTreePS(covf,28,out,covn);
		out<<"showpage"<<endl;

		std::copy(&stats2.cov[0][0],&stats2.cov[0][0]+28*28,covf);
		out<<"/Times-Roman findfont 3 scalefont setfont"<<endl;
		out<<"20 20 translate 3 4 scale"<<endl;
		out<<"100 110 moveto (value orderth squareroot taken) show"<<endl;
		drawCorrTreePS(covf,28,out,covn);
		out<<"showpage"<<endl;

		std::copy(&stats3.cov[0][0],&stats3.cov[0][0]+28*28,covf);
		out<<"/Times-Roman findfont 3 scalefont setfont"<<endl;
		out<<"20 20 translate 3 4 scale"<<endl;
		out<<"100 110 moveto (domain orderth squareroot taken) show"<<endl;
		drawCorrTreePS(covf,28,out,covn);
		out<<"showpage"<<endl;		

		cout<<"ready."<<endl;

	}
};



struct momStats{

	typedef Objectdb::stats stat;
	//calculate variance , gravitycenter, bbox of dataset

	stat ptstats[3],
	momstats[28], //direct statistis of moments
	momstats2[28],//statistics of moments rooted by domain order
	momstats3[28];//statistics of moments rooted by value order

	mytype r;

	typedef momentSet<mytype,28,0> tmom;
	typedef tensorSet3D<4,mytype> ttset;

	int rexp; //ln2(radius) of balls for mom.comp.

	int npts;
	//the moments
	vector<tmom> moms;



	/**
	 * create a stats from dump file fname
	 *\param valdomsi:
	 * should we norm the moments in value+domain ?
	 */
	momStats(const char*fname,bool valdomsi=false,bool i2=false)
	{
		vector<mytype> pts;
		cout<<"reading in data from "<<fname<<"..."<<endl;
		const char* pos= strchr(fname,'.');

		memset(momstats,0,sizeof(momstats));
		memset(momstats2,0,sizeof(momstats));
		memset(momstats3,0,sizeof(momstats));

		PointCloudIO::read(pts,fname);

		if(i2)
			init2(pts,valdomsi);
		else
			init(pts,valdomsi);
	}

	void init(const vector<mytype> &pts,bool valdomsi)
	{
		static const int n=1000;
		moms.resize(n);
		npts=pts.size();
		if(!npts)
			return;

		vector<mytype>::const_iterator it;


		for(it=pts.begin();it<pts.end();)
		{
			for(int i=0;i<3;++i,++it){
				ptstats[i].add(*it);
			}
		}


		//radius should be 1/2 of the std deviation
		r = sqrt(ptstats[0].variance()
				+ptstats[1].variance()
				+ptstats[2].variance())*0.5;
		frexp(r*0.5,&rexp);

		//the radius, rounded down to next lower power of 2
		r = ldexpl(1,rexp);

		//radius,squared
		mytype rsq=r*r;


		//the centers
		vector<mytype> ct(n*3);


		for(int i=0;i<n*3;i+=3)
		{
			int id= 3*( rand()*(pts.size()/3)/RAND_MAX );
			ct[i]=pts[id];
			ct[i+1]=pts[id+1];
			ct[i+2]=pts[id+2];
		}


		tmom momset;




		//compute moments:
		mytype x[3];

		vector<mytype>::const_iterator cit;

		vector<ttset > tsets(n);
		vector<ttset>::iterator tsit;

		for(tsit=tsets.begin();tsit!=tsets.end();++tsit)
			for(int j=0;j<ttset::ncomps;++j)
				tsit->A[j]=0;

		cout<<"computin tsets..."<<endl;



		//for every point: add it to every tensor set
		//whose center is less than r away
		for(it=pts.begin();it!=pts.end();it+=3)
		{

			//for every ball: check if pt is inside
			for(
					cit=ct.begin(),tsit=tsets.begin();
					cit!=ct.end();
					cit+=3,++tsit)
			{
				x[0]=it[0]-cit[0];
				x[1]=it[1]-cit[1];
				x[2]=it[2]-cit[2];

				if(x[0]*x[0]+x[1]*x[1]+x[2]*x[2] < rsq)
				{
					tsit->add3dcoords(x /* &*it */);
				}
			}

			int i=(it-pts.begin())/3 ;
			int nn=pts.size();
			if((i&0xfff)==0)
				cout<<i<<'/'<<nn<<'('<<i*100/nn<<"%)\r"<<flush;
		}

		vector<tmom>::iterator mit;
		cout<<"computin moments..."<<endl;

		//compute the moments
		tensorSet3D<4> newA;
		for(tsit=tsets.begin(),mit=moms.begin();
				mit!=moms.end();++mit,++tsit)
		{

			//newA is the translated version of tsit->A
			tsit->getCenter(x);
			x[0]=-x[0];x[1]=-x[1];x[2]=-x[2];
			newA.translated(*tsit,x);
			mit->compute(newA.A);
			if(valdomsi)
			{
				mit->scaleValue(sqrt(mit->M[0]/mit->M[21]));
				mit->scaleDomain(mit->M[21]/(mit->M[0]*mit->M[0]));
			}
			else{
				mit->scaleValue(1.0/newA.A[0]);
				mit->scaleDomain(1.0/r);
			}

		}


		//now: make statistics with the moments.

		for(int i=0;i<28;++i)
		{
			stat &ms=momstats[i],ms2=momstats2[i];
			mytype pw=1.0/tmom::metainfo.me[i].dord;
			for(mit=moms.begin();mit!=moms.end();++mit)
			{
				mytype v=mit->M[i];
				ms.add(v);
				ms2.add((v>0?1:-1)*pow(fabs(v),pw));
			}
		}

	}

	void init2(const vector<mytype> &pts,bool valdomsi)
	{
		static const int n=100000;
		moms.resize(n);
		npts=pts.size();
		if(!npts)
			return;

		vector<mytype>::const_iterator it;


		for(it=pts.begin();it<pts.end();)
		{
			for(int i=0;i<3;++i,++it){
				ptstats[i].add(*it);
			}
		}


		//radius should be 1/2 of the std deviation
		r = sqrt(ptstats[0].variance()
				+ptstats[1].variance()
				+ptstats[2].variance())*0.5;
		frexp(r*0.5,&rexp);

		//the radius, rounded down to next lower power of 2
		r = ldexpl(1,rexp);

		//radius,squared
		mytype rsq=r*r;

		cout<<"radius="<<r<<endl;
		//the centers
		/*
	vector<mytype> ct(n*3);


	for(int i=0;i<n*3;i+=3)
	{
		int id= 3* ( ( rand()*(pts.size()/3)/RAND_MAX )%(pts.size()/3)  );
		ct[i]=pts[id];
		ct[i+1]=pts[id+1];
		ct[i+2]=pts[id+2];
	}
		 */

		tmom momset;


		TsTree t(pts);



		vector<tensorSet3D<> > in,bo;

		cout<<"getting tsets..."<<flush;
		t.numbersummed=0;
		t.ptssummed=0;
		t.getBallTsets(in,bo,0,(int)(r/t.getScale()));
		cout<<"ready: #tsets summed:"<<t.numbersummed;
		if(t.ptssummed=0)
			cout<<" #points summed:"<<t.ptssummed;

		cout<<endl;

		cout<<"computing "<<in.size()<<" moments ..."<<flush;
		mynum ptit[3];
		moms.resize(in.size());
		for(int i=0;i<in.size();++i)
		{
			in[i].getCenter(ptit);
			ptit[0]=-ptit[0];ptit[1]=-ptit[1];ptit[2]=-ptit[2];
			in[i].translate(ptit);

			moms[i].compute(in[i].A);
		}
		cout<<"ready. "<<flush;


		cout<<(valdomsi?"normalize...":"scale...")<<flush;

		vector<tmom>::iterator mit;
		for(mit=moms.begin();mit!=moms.end();++mit)
		{
			if(valdomsi){
				double vscale=sqrt(mit->M[0]/mit->M[21]);
				mit->scaleDomain(mit->M[21]/(mit->M[0]*mit->M[0]));
				mit->scaleValue(vscale);
			}
			else
			{

				//cout<< in[mit-moms.begin()].A[0]<<' '<<r<<endl;
				mit->scaleDomain(1.0/r);
				mit->scaleValue(1.0/in[mit-moms.begin()].A[0]);
			}
		}
		cout<<"ready. "<<flush;

		cout<<"compute stats..."<<flush;

		//now: make statistics with the moments.
		for(int i=0;i<28;++i)
		{
			stat &ms=momstats[i],&ms2=momstats2[i],&ms3=momstats3[i];

			mytype pwdord=1.0/tmom::metainfo.me[i].dord;
			mytype pwvord=1.0/tmom::metainfo.me[i].vord;

			for(mit=moms.begin();mit!=moms.end();++mit)
			{
				mytype v=mit->M[i];
				ms.add(v);
				ms2.add((v>0?1:-1)*pow(fabs(v),pwvord));
				ms3.add((v>0?1:-1)*pow(fabs(v),pwdord));
			}
		}
		cout<<"ready."<<endl;

	}



	inline void write(ostream &out)
	{
		out<<"//point statistics"<<endl;
		out<<"P"<<npts<<endl;
		out<<"//Middel std.dev min max"<<endl;

		for(int i=0;i<3;++i){
			stat&s=ptstats[i];
			out<<' '<<s;
		}

		out<<"R"<<r<<endl;
		out<<"//Moment Statistics:"<<endl;
		out<<"S"<<tmom::n<<endl;
		out<<"//Middel Std.dev Min Max"<<endl;
		for(int i=0;i<tmom::n;++i)
		{
			stat&s = momstats[i];
			out<< "//"<<tmom::metainfo.me[i].des<<":"<<endl;
			out<<' '<<s<<endl;
		}

		out<<"//Moment Statistics2:"<<endl;
		out<<"T"<<tmom::n<<endl;
		out<<"//Middel Std.dev Min Max"<<endl;
		for(int i=0;i<tmom::n;++i)
		{
			stat&s = momstats2[i];
			out
			<< "//"<<tmom::metainfo.me[i].des<<":"<<endl;
			out<<' '<<s<<endl;
		}

		out<<"//Moment Statistics3:"<<endl;
		out<<"U"<<tmom::n<<endl;
		out<<"//Middel Std.dev Min Max"<<endl;
		for(int i=0;i<tmom::n;++i)
		{
			stat&s = momstats3[i];
			out
			<< "//"<<tmom::metainfo.me[i].des<<":"<<endl;
			out<<' '<<s<<endl;
		}

		out<<"//number_of_moment_sets vals_per_moment"<<endl;
		out<<"M"<<moms.size()<<"  " <<tmom::n<<endl;
		for(int i=0;i<moms.size();++i)
		{
			for(int j=0;j<tmom::n;++j)
				out<<' '<<moms[i].M[j];
			out<<endl;
		}
	}

};



//draw covariance matrix
inline void drawCorrTreePS_old(const char*infile,ostream&out)
{
	cout<<"trying to write covariance to ps file..."<<flush;
	ifstream in(infile);
	typedef momStats::tmom tmom;

	int a4w = int(21 *72 / 2.54);
	int a4h = int(29 *72 / 2.54);

	char buf[1024];
	char name[1024];
	int dummy=0;

	mynum avg[28],avg2[28],avg3[28];
	mynum var[28],var2[28],var3[28];
	mynum cov[28][28],cov2[28][28],cov3[28][28];


	while(!in.eof()&&!in.fail()){

		do{	in.getline(buf,1023);	} while(buf[0]!='F' && ! in.eof() );		
		if(in.eof()) break;
		sscanf(buf,"F%256s",name);

		cout<<"reading moments of "<<name<<"..."<<flush;
		//read in moments
		do{ in.getline(buf,1023);} while(buf[0]!='M' && ! in.eof() );
		int nmoms;
		sscanf(buf,"M%d%d",&nmoms,&dummy);
		vector<mynum> moms(nmoms*dummy),
				moms2(nmoms*dummy),
				moms3(nmoms*dummy);
		for(int i=0;i<moms.size();++i)
			in>>moms[i];

		cout<<"ready."<<flush;
		cout<<"calulating stats..."<<flush;
		//calculate values and averages:			
		for(int i=0;i<28;++i)
			avg[i]=avg2[i]=avg3[i]=var[i]=var2[i]=var3[i]=0;

		for(int i=0;i<moms.size();)
		{
			for(int J=0;J<28;++J,++i){
				mytype pwdord=1.0/tmom::metainfo.me[J].dord; 
				mytype pwvord=1.0/tmom::metainfo.me[J].vord; 
				mynum v=moms[i]; avg[J]+=moms[i];
				moms2[i]=(v>0?1:-1)*pow(fabs(v),pwvord);
				avg2[J]+=moms2[i];
				moms3[i]=((v>0?1:-1)*pow(fabs(v),pwdord));					
				avg3[J]+=moms3[i];
			}		
		}

		for(int i=0;i<28;++i)
		{avg[i]/=nmoms;avg2[i]/=nmoms;avg3[i]/=nmoms;}


		//calculate variance and subtract average
		for(int i=0;i<moms.size();)
		{
			for(int J=0;J<28;++J,++i){
				moms[i]-=avg[J];
				var[J]+= moms[i]*moms[i];
				moms2[i]-=avg2[J];
				var2[J]+= moms2[i]*moms2[i];
				moms3[i]-=avg3[J];
				var3[J]+= moms3[i]*moms3[i];
			}
		}

		for(int i=0;i<28;++i)
		{
			var[i]/=nmoms-1;var2[i]/=nmoms-1;var3[i]/=nmoms-1;
			//cout<<"average["<<i<<"]="<<var[i]<<endl;
		}

		memset(cov,0,sizeof(cov));
		memset(cov2,0,sizeof(cov));
		memset(cov3,0,sizeof(cov));

		for(int i=0;i<moms.size();i+=tmom::n)
		{
			for(int J=0;J<28;++J)for(int k=0;k<=J;++k)
			{
				cov[J][k]+= moms[i+J]*moms[i+k];			
				cov2[J][k]+= moms2[i+J]*moms2[i+k];			
				cov3[J][k]+= moms3[i+J]*moms3[i+k];			
			}
		}

		for(int J=0;J<28;++J)for(int k=0;k<=J;++k)
		{
			cov[J][k] /= (nmoms-1)*sqrt(var[J]*var[k]);
			cov[k][J] = cov[J][k];
			cov2[J][k] /= (nmoms-1)*sqrt(var2[J]*var2[k]);
			cov2[k][J] = cov2[J][k];
			cov3[J][k] /= (nmoms-1)*sqrt(var3[J]*var3[k]);
			cov3[k][J] = cov3[J][k];

			//cout<<"cov["<<k<<"]["<<J<<"]="<<cov[k][J];
		}


		cout<<"ready. writing to ps file..."<<flush;

		out<<"/Times-Roman findfont 20 scalefont setfont"<<endl;
		out<< 30 << ' '<<a4h -60<<" moveto ("<<name<<") show"<<endl;

		out<<"gsave 30 30 translate"<<endl;
		covps(cov,out,a4h/4/28);
		out<<"grestore"<<endl;

		out<<"gsave 30 "<<30+a4h/4<<" translate"<<endl;
		out<<a4h/4<<" 0 moveto (with pow\\(x,1/value_order\\)) show"<<endl;
		covps(cov2,out,a4h/4/28);
		out<<"grestore"<<endl;

		out<<"gsave 30 "<<30+a4h/2<<" translate"<<endl;

		out<<a4h/4<<" 0 moveto (with pow\\(x,1/domain_order\\)) show"<<endl;
		covps(cov3,out,a4h/4/28);
		out<<"grestore"<<endl;
		out<<"showpage"<<endl;


		float covf[28*28];
		string covn[28];
		for(int i=0;i<28;++i){
			string s = tmom::metainfo.me[i].des;

			for(int j=0;j<s.size()-1;++j)
			{
				if(s[j]==')')s[j]=']';
				if(s[j]=='(')s[j]='[';
				if(s[j]=='\\')s[j]=' ';
			}
			covn[i]=s;
		}


		std::copy(&cov[0][0],&cov[0][0]+28*28,covf);
		out<<"/Times-Roman findfont 3 scalefont setfont"<<endl;
		out<<"20 20 translate 3 4 scale"<<endl;
		out<<"100 110 moveto (direct moments) show"<<endl;
		drawCorrTreePS(covf,28,out,covn);
		out<<"showpage"<<endl;

		std::copy(&cov2[0][0],&cov2[0][0]+28*28,covf);
		out<<"/Times-Roman findfont 3 scalefont setfont"<<endl;
		out<<"20 20 translate 3 4 scale"<<endl;
		out<<"100 110 moveto (value orderth squareroot taken) show"<<endl;
		drawCorrTreePS(covf,28,out,covn);
		out<<"showpage"<<endl;

		std::copy(&cov3[0][0],&cov3[0][0]+28*28,covf);
		out<<"/Times-Roman findfont 3 scalefont setfont"<<endl;
		out<<"20 20 translate 3 4 scale"<<endl;
		out<<"100 110 moveto (domain orderth squareroot taken) show"<<endl;
		drawCorrTreePS(covf,28,out,covn);
		out<<"showpage"<<endl;		

		cout<<"ready."<<endl;
	}

}




inline void writeCorrelationPS(const char*infile,ostream&out)
{
	ifstream in(infile);

	momStats2 st;
	while(!in.eof())
	{
		if(st.read(in))
			st.writecorrps(out);	
		else 
			break;
	}
}











struct statist
{
	char nam[80];
	mytype avg,sigma,min,max;
};

struct statisti{
	char nam[80];
	statist s[28];
	statist s2[28];
	statist s3[28];

	double r;
};


inline void writestatsps(const char*infile,ostream&out)
{
	ifstream in(infile);

	if(!in.is_open())
		return;

	char buf[1024];
	int dummy;

	vector<statisti> s;



	cout<<"reading in stats"<<endl;
	cout.flush();

	while(!in.eof()){
		s.push_back(statisti());
		statisti &ss=s.back();


		do{in.getline(buf,1024,'\n');}
		while(buf[0]!='F' && !in.eof());
		if(in.eof())break;

		sscanf(buf,"F%s",ss.nam);

		do{in.getline(buf,1024,'\n');}
		while(buf[0]!='R' && !in.eof());
		if(in.eof())break;

		sscanf(buf,"R%lf",&ss.r);

		char stchar[]={'S','T','U'};
		statist*st[]={ss.s,ss.s2,ss.s3};

		for(int k=0;k<3;++k){

			do{ in.getline(buf,1024,'\n'); }
			while(buf[0]!=stchar[k] && !in.eof());
			if(in.eof())break;		
			for(int i=0;i<28;++i){
				statist &s=st[k][i];
				//overjump comments
				do{in.getline(buf,1024,'\n');}while(buf[0]=='/');		
				sscanf(buf,"%lf%lf%lf%lf",&s.avg,&s.sigma,&s.min,&s.max);
				cout<<s.avg<<' '<<s.sigma<<' '<<s.min<<' '<<s.max<<endl;
			}
		}

	}

	cout<<"read "<<s.size()<<"stats"<<endl;
	if(s.size()>1)s.pop_back();

	mytype momi[3][28],moma[3][28];

	for(int j=0;j<s.size();++j)
	{
		statist * xx[]= {s[j].s,s[j].s2,s[j].s3};
		for(int k=0;k<3;++k)for(int i=0;i<28;++i)
		{
			mytype &mi = momi[k][i],&ma=moma[k][i];
			statist &x = xx[k][i];
			if(j==0){
				mi=x.min;ma=x.max;
			}
			else{
				if(mi>x.min) mi=x.min;
				if(ma<x.max) ma=x.max;
			}					
		}
	}




	int a4w = int(21 *72 / 2.54);
	int a4h = int(29 *72 / 2.54);
	int lh=a4h/(s.size()+1);

	out<<"/Times-Roman findfont 20 scalefont setfont"<<endl;

	char * names[]={"normal","value order'th root","domain order'th root"};
	for(int i=0;i<28;++i)
	{
		for(int k=0;k<3;++k)
		{
			mytype sc = a4w/(moma[k][i]-momi[k][i]),
					st = momi[k][i]*sc;

			vector<statisti>::iterator j;
			for(j=s.begin();j!=s.end();++j)
			{

				statist & x = (k==0 ? j->s[i] : ( k==1 ? j->s2[i] :j->s3[i] ));


				int mi=int(sc*x.min-st),
						ma=int(sc*x.max-st),
						avg=int(sc*x.avg-st),
						si=int(sc*x.sigma);

				out<< avg-si<<" 2 moveto "<<avg+si<<" 2 lineto "
						<< avg+si<<" "<<lh-2<< " lineto "
						<< avg-si<<" "<<lh-2<< " lineto 0.5 setgray fill "
						<<endl;

				//draw minmax:
				out<< mi <<" 2 moveto "
						<< ma << " 2 lineto "
						<< ma<<" "<< lh-2 <<" lineto "
						<< mi<<" "<<lh-2<<" lineto closepath 0 setgray stroke"
						<<endl;

				//draw med:
				out<< avg<<" 0 moveto "<<avg<<" "<<lh<<" lineto stroke "<<endl;



				out<<2<<" 2 moveto ("<<j->nam<<") show"
						<< "( r: "<<j->r
						<<" range:["<<x.min<<","<<x.max<<"]) show "<<endl;

				out<<"0 "<<lh<< " translate"<<endl; 		
			}

			out<<"0 0 moveto ("<< faceMomSet::metainfo.me[i].des 
					<<" "<<names[k]<<") show"<<endl;


			out<<"showpage"<<endl;
		}
	}

}



inline void writestatsps2(const char*infile,ostream&out)
{
	ifstream in(infile);

	if(!in.is_open())
		return;

	char buf[1024];
	int dummy;


	vector<momStats2> s;



	cout<<"reading in stats"<<endl;
	cout.flush();

	momStats2 ms;

	while(!in.eof()){
		if(ms.read(in)&&!in.fail())
			s.push_back(ms);
		else break;
	}

	cout<<"read "<<s.size()<<"stats"<<endl;

	mytype momi[3][28],moma[3][28];

	for(int j=0;j<s.size();++j)
	{
		vecstats<28> * xx[]= {&s[j].stats,&s[j].stats2,&s[j].stats3};
		for(int k=0;k<3;++k)for(int i=0;i<28;++i)
		{
			mytype &mi = momi[k][i],&ma=moma[k][i];
			mytype xmin=xx[k]->min[i],xmax=xx[k]->max[i];
			if(j==0){
				mi=xmin;ma=xmax;
			}
			else{
				if(mi>xmin) mi=xmin;
				if(ma<xmax) ma=xmax;
			}					
		}
	}




	int a4w = int(21 *72 / 2.54)-40;
	int a4h = int(29 *72 / 2.54)-40;
	int lh=a4h/(s.size()+1);


	out<<"/Times-Roman findfont 20 scalefont setfont"<<endl;
	char * names[]={"normal","value order'th root","domain order'th root"};
	for(int i=0;i<28;++i)
	{
		for(int k=0;k<3;++k)
		{

			out<<"20 20 translate"<<endl;

			mytype sc = a4w/(moma[k][i]-momi[k][i]),
					st = momi[k][i]*sc;

			vector<momStats2>::iterator j;
			for(j=s.begin();j!=s.end();++j)
			{

				vecstats<28> & x = (k==0 ? j->stats : 
						( k==1 ? j->stats2 :j->stats3 ));


				int mi=int(sc*x.min[i]-st),
						ma=int(sc*x.max[i]-st),
						avg=int(sc*x.avg[i]-st),
						si=int(sc*x.sigma[i]);

				out<< avg-si<<" 2 moveto "<<avg+si<<" 2 lineto "
						<< avg+si<<" "<<lh-2<< " lineto "
						<< avg-si<<" "<<lh-2<< " lineto 0.5 setgray fill "
						<<endl;

				//draw minmax:
				out<< mi <<" 2 moveto "
						<< ma << " 2 lineto "
						<< ma<<" "<< lh-2 <<" lineto "
						<< mi<<" "<<lh-2<<" lineto closepath 0 setgray stroke"
						<<endl;

				//draw med:
				out<< avg<<" 0 moveto "<<avg<<" "<<lh<<" lineto stroke "<<endl;

				out<<2<<" 2 moveto ("<<nameps(j->name)<<") show"
						<<"( range:["<<x.min[i]<<","<<x.max[i]<<"] )show "<<endl;

				out<<"0 "<<lh<< " translate"<<endl; 		
			}

			out<<"0 10 moveto ("<< faceMomSet::metainfo.me[i].des 
					<<" "<<names[k]<<") show"<<endl;


			out<<"showpage"<<endl;
		}
	}
}


inline void createstats(int argc, const char*argv[])
{
	ofstream stf("stats.txt",ofstream::app);
	ofstream stf2("stats2.txt",ofstream::app);
	const char* defaultname
	= "c:\\u\\max\\models\\da\\dino.raw";


	for(int i=0;i<argc;++i)		
	{
		printf("%S",argv[i]);
	}
	if(argc<2){
		stf<<"F"<<defaultname<<endl;
		momStats(defaultname,false,true).write(stf);
	}
	else
		for(int i=1;i<argc;++i)
		{
			char buf[80];
			sprintf(buf,"%s",argv[i]);
			cout<<buf<<"("<<i<<'/'<<argc<<")"<<endl;
			stf<<"F"<<buf<<endl;
			momStats x(buf,false,true);
			x.write(stf);			
			vector<mynum> xx(x.moms.size()*momStats::tmom::n);
			vector<mynum>::iterator it2=xx.begin();
			for(int i=0;i<x.moms.size();++i)
				it2=std::copy(&x.moms[i].M[0],&x.moms[i].M[0]+28 ,it2);
			momStats2 y(xx,false,buf);
			y.write(stf2);

		}

	stf.close();
	stf2.close();
}


inline void testTsTreeBuild(const vector<long double>&points)
{
	clock_t start,sum;
	int it;
	cout<<"building tree of "<<points.size()/3<<" points"<<endl;

	for(int ptperleaf=8;ptperleaf<128;ptperleaf=(ptperleaf *5 +3 ) / 4)
	{
		cout<<ptperleaf<<" desired points per leaf"<<endl;
		start=clock();it=0;
		while((sum=clock()-start)<1000)
		{
			TsTree2 *t=new TsTree2(points,0,ptperleaf);++it;
			//if(it==1)cout<<"built tree of height"<<t->getDepth()<<'\n';
			delete t;
		}
		cout<<" used "<<1e3 * sum/it/CLOCKS_PER_SEC<<" ms per building"<<endl;

		TsTree2 t2(points,0,ptperleaf);

		cout<<"testing speed TsTree2::sumInCircle"<<endl;
		start=clock();

		TsTree2::tset ts2,bord2;
		int sumnum=0;
		int trynum=1000;
		long double ra = ldexp(t2.getScale(),t2.getDepth()-2);
		do{
			for(int i=0;i<trynum;++i)
				t2.sumInCircle(ts2,bord2,(long double*)&points[0],
						(long double)ra,-1);
			sumnum+=trynum;
			sum=clock()-start;
		}while(sum<1000);
		cout<<"used "<<1e3 * sum/CLOCKS_PER_SEC/sumnum<<" ms per summing\n"<<endl;

	}



	cout<<"test speed building TsTree2"<<endl;
	start=clock();
	it=0;
	while((sum=clock()-start)<2000)
	{
		TsTree2 *t=new TsTree2(points,0.01,0);++it;
		if(it==1)cout<<"built tree of max. height"<<t->getDepth()<<'\n';
		delete t;
	}
	cout<<" used "<<1e3 * sum/it/CLOCKS_PER_SEC<<" ms per building"<<endl;

	cout<<"test speed building TsTree method 1"<<endl;
	start=clock();
	it=0;
	while((sum=clock()-start)<2000)
	{
		TsTree *t=new TsTree(points,0.01,0);++it;
		if(it==1)cout<<"built tree of height"<<t->getDepth()<<'\n';
		delete t;
	}
	cout<<" used "<<1e3 * sum/it/CLOCKS_PER_SEC<<" ms per building"<<endl;

	cout<<"...with method2"<<endl;
	start=clock();it=0;
	while((sum=clock()-start)<2000)
	{
		TsTree *t=new TsTree(points,0.01,1);++it;
		if(it==1)cout<<"built tree of height"<<t->getDepth()<<'\n';
		delete t;
	}
	cout<<" used "<<1e3 * sum/it/CLOCKS_PER_SEC<<" ms per building"<<endl;






}

inline void testTsTreesumCircle(const vector<long double>&points)
{			
	clock_t start,sum;

	cout<<"building tree of "<<points.size()/3<<" points"<<endl;

	start=clock();
	TsTree2 t2(points,0.02);
	sum=clock()-start;
	cout<<" used "<<1e3 * sum/CLOCKS_PER_SEC<<" ms for building TsTree2"<<endl;
	cout<<"tree height:"<<t2.getDepth()<<endl;

	start=clock();
	TsTree t(points,0.02,1);
	sum=clock()-start;
	cout<<" used "<<1e3 * sum/CLOCKS_PER_SEC<<" ms for building TsTree"<<endl;
	cout<<"tree height:"<<t.getDepth()<<endl;

	long double ra = ldexp(t2.getScale(),t2.getDepth()-2);

	tensorSet3D<> tstest,tstest2,tsouter,ts,bord,ts2,bord2;

	t2.sumInCircle(ts,bord,(long double*)&points[0],(long double)ra,-1);

	t2.sumInCircle(ts2,bord2,(long double*)&points[0],(long double)ra,-1);

	tsouter=ts;tsouter+=bord;
	vector<mynum>::const_iterator it,
	ptit=points.begin();


	long double center[3];
	const long double 
	rb0=t.getScale(),
	*mi=t.getOrig();
	for(int i=0;i<3;++i)
		center[i] = int( (points[i] - mi[i])*2.0/rb0 ) 
		* rb0/2.0 +mi[i];

	tstest.reset();
	tstest2.reset();
	mynum buf,lsq;

	for(it=points.begin();it!=points.end();it+=3)
	{
		lsq=ra*ra;
		buf=it[0]-center[0];lsq-=buf*buf;
		buf=it[1]-center[1];lsq-=buf*buf;
		buf=it[2]-center[2];lsq-=buf*buf;
		if(lsq>0)
			tstest.add3dcoords(&*it);															

		lsq=ra*ra;
		buf=it[0]-points[0];lsq-=buf*buf;
		buf=it[1]-points[1];lsq-=buf*buf;
		buf=it[2]-points[2];lsq-=buf*buf;
		if(lsq>0)
			tstest2.add3dcoords(&*it);															
	}


	mynum maxdist =0,maxdist2=0,len2=0;
	mynum len= 0;
	for(int i=0;i<tstest.ncomps;++i){
		len+=fabs(tstest.A[i])+fabs(ts.A[i]);
		maxdist+=fabs(tstest.A[i]-ts.A[i]);
		len2+=fabs(tstest2.A[i]+ts2.A[i]);
		maxdist2+=fabs(tstest2.A[i]-ts2.A[i]);
	}
	cout<<"ERROR TsTree2 integer procedure: "<<" : "<<100.*maxdist/len<<"% of "<<len<<endl;
	cout<<"ERROR TsTree2 float procedure: "<<" : "<<100.*maxdist2/len2<<"% of "<<len2<<endl;


	int trynum=1000;

	int sumnum;


	cout<<"testing speed float procedure TsTree2:"<<endl;
	start=clock();

	sumnum=0;
	do{
		for(int i=0;i<trynum;++i)
			t2.sumInCircle(ts2,bord2,(long double*)&points[0],(long double)ra,0);
		sumnum+=trynum;
		sum=clock()-start;
	}while(sum<500);
	cout<<" used "<<1e3 * sum/CLOCKS_PER_SEC/sumnum<<" ms per summing"<<endl;



	cout<<"testing speed float procedure TsTree:"<<endl;
	start=clock();
	sumnum=0;
	do{
		for(int i=0;i<trynum;++i)
			t.sumInCircle(ts2,bord2,(long double*)&points[0],(long double)ra,0);
		sumnum+=trynum;
		sum=clock()-start;
	}while(sum<500);
	cout<<" used "<<1e3 * sum/CLOCKS_PER_SEC/sumnum<<" ms per summing"<<endl;


	trynum/=10;

	cout<<"testing speed float procedure TsTree with points summing:"<<endl;
	start=clock();
	sumnum=0;
	do{
		for(int i=0;i<trynum;++i)
			t.sumInCircle(ts2,bord2,(long double*)&points[0],(long double)ra,-1);
		sumnum+=trynum;
		sum=clock()-start;
	}while(sum<500);
	cout<<" used "<<1e3 * sum/CLOCKS_PER_SEC/sumnum<<" ms per summing"<<endl;


	cout<<"testing speed float procedure TsTree2 with points summing:"<<endl;
	start=clock();
	sumnum=0;
	do{
		for(int i=0;i<trynum;++i)
			t2.sumInCircle(ts2,bord2,(long double*)&points[0],(long double)ra,-1);
		sumnum+=trynum;
		sum=clock()-start;
	}while(sum<500);
	cout<<" used "<<1e3 * sum/CLOCKS_PER_SEC/sumnum<<" ms per summing"<<endl;


	/*
	cout<<"building kdtree"<<endl;
	vector<int> bla(points.size());
	for(int i=0;i<points.size();++i)
		bla[i]=(int)points[i]*1000;
	Momkdtree_templ<int> k(bla,3);
	cout<<"testing speed kdtree procedure with no summing:"<<endl;
	vector<int> ids;
	start=clock();
	sumnum=0;
	int cnter[]={points[0]*1000,points[1]*1000,points[2]*1000};
	do{
	for(int i=0;i<trynum;++i)
		k.momIdsInBall(ids,cnter,(int)(ra*1000));
		sumnum+=trynum;
		sum=clock()-start;
	}while(sum<500);
	cout<<" used "<<1e3 * sum/CLOCKS_PER_SEC/sumnum<<" ms per summing"<<endl;

	 */



	cout<<"testing multiple radii summing TsTree2:"<<endl;

	static const int nr=1000;
	mynum radii[nr]; TsTree::tset tsets[nr], tsets2[nr];
	for(int i=0;i<nr;++i)
		radii[i] = pow(2.0f, - (nr-i-1) * 3.0f/nr);


	t2.sumInCircle(tsets,nr,(long double*)&points[0],radii);

	//for reference: the original procedure
	for(int i=0;i<nr;++i)
		t.sumInCircle(tsets2[i],bord2,(long double*)&points[0],radii[i],-1);


	//compare:
	for(int j=0;j<nr;++j){
		//cout<<" r="<<radii[j]<<flush;
		len=0;maxdist=0;	
		for(int i=0;i<tstest.ncomps;++i){

			//cout<<"tsets A["<<i<<"]="<<tsets[j].A[i]<<endl;
			//cout<<"tsets2 A["<<i<<"]="<<tsets2[j].A[i]<<endl;
			len+=fabs(tsets[j].A[i])+fabs(tsets2[j].A[i]);
			maxdist+=fabs(tsets[j].A[i]-tsets2[j].A[i]);
		}

		if(maxdist>1e-9*len)
			cout<<"error: tset No. "<<j<<",r="<<radii[j]<<":"<<100*maxdist/len<<"% "<<endl;
		else cout<<"."<<flush;
	}

	cout<<"..ready."<<endl;

	cout<<"testing speed mult-radii-summing TsTree2 with "<<nr<<" radii"<<endl;
	start=clock();
	sumnum=0;
	trynum = trynum/nr +1;
	do{
		for(int i=0;i<trynum;++i)
			t2.sumInCircle(tsets,nr,(long double*)&points[0],radii);
		sumnum+=trynum;
		sum=clock()-start;
	}while(sum<1000);
	cout<<" used "<<1e3 * sum/CLOCKS_PER_SEC/sumnum<<" ms per summing"<<endl;

	cout<<"testing speed of the same task with single-radius-summing"<<endl;
	start=clock();
	sumnum=0;
	do{
		for(int i=0;i<trynum;++i)
			//for reference: the original procedure
			for(int i=0;i<nr;++i)
				t2.sumInCircle(tsets2[i],bord2,(long double*)&points[0],radii[i],-1);
		sumnum+=trynum;
		sum=clock()-start;
	}while(sum<1000);
	cout<<" used "<<1e3 * sum/CLOCKS_PER_SEC/sumnum<<" ms per summing"<<endl;


}


inline void outputranks(vector<int> &ranks,ostream&out)
{
	for(int i=1;i<ranks.size();++i)
		ranks[i]+=ranks[i-1];

	if(ranks.empty() || ranks.back()==0)return;

	int i;
	for(i=0;i<ranks.size()&& i<10 && ranks[i]<ranks.back();++i)
		out<< ranks[i]*100/ranks.back()<< "% rank<= "<<i+1<<", ";

	if(i<ranks.size())
		out<< ranks[i]*100/ranks.back()<< "% rank<= "<<i+1<<endl;
}


inline void testTsTree(const char*nam
		="c:\\u\\max\\models\\denker\\da\\dino.raw")
{


	vector<long double> pts;
	cout<<"reading file "<<nam<<" ..."<<endl;

	PointCloudIO::read(pts,nam);


	testTsTreeBuild(pts);

	for(int i=0;i<300&&i<pts.size();++i)
		cout<<pts[i]<<' '<<flush;

	//norm points so they are all in the range -1,1

	long double min[3],max[3];
	vector<long double>::iterator it=pts.begin();

	for(int i=0;i<3;++i,++it)
		min[i]=max[i]=*it;

	while(it!=pts.end()){
		for(int i=0;i<3;++i,++it)
		{
			if(min[i]>*it)min[i]=*it;
			if(max[i]<*it)max[i]=*it;
		}
	}

	long double mul[3],add[3];
	for(int i=0;i<3;++i)
	{
		mul[i]=2.0/(max[i]-min[i]);
		add[i]=-(max[i]+min[i])*0.5;
	}

	it=pts.begin();
	while(it!=pts.end())
		for(int i=0;i<3;++i,++it)
			*it= (*it+add[i])*mul[i];



	testTsTreesumCircle(pts);


	cout<<"building tree..."<<endl;
	TsTree *t=new TsTree(pts);
	cout<<"fractal dim=" << t->fractdim()<<endl;


	vector<TsTree::tset> inner,border;

	t->getBallTsets(inner,border,0,1);
	cout<<inner.size()<<"tsets computed"<<endl;



	delete t;
	cout<<"after del"<<endl;

	cout<<"now testing for memory leaks"<<endl;
	for(int i=0;i<100000;++i)
	{
		cout<<"create..."<<endl;
		t=new TsTree(pts,0.002,i&1);
		cout<<"delete..."<<endl;
		delete t;
	}

}




struct matcheslt{
	bool operator () (
			const Objectdb::match&a,
			const Objectdb::match&b)
	{
		return a.dist<b.dist;
	}
};

struct toremove {
	int id;
	toremove(int i):id(i){}
	bool operator () (const Objectdb::match&a)
	{
		return a.dbmom.objectid!=id || a.dbmom.radius>64;
	}
};

/** get filenames in dir ending with ext*/
inline void getnames(vector<string>&names,const std::string &dir , const std::string& ext)
{
	DIR* d = opendir(dir.c_str());	
	assert(d!=0);
	struct dirent* ent = readdir(d);	
	while(ent)
	{
		const char* s=strstr(ent->d_name,ext.c_str());
		if(s && strlen(s)==ext.length() )
			names.push_back(dir+ent->d_name);
		ent=readdir(d);
	}
	closedir(d);
}


inline int testfind(Objectdb&db,const vector<mynum> &pts,int oid,
		ostream&out,
		mynum *likelihoodoid=0,
		int *rankoid=0,
		const vector<mynum> * origpts=0,
		const vector<int>*classes=0
)
{
	vector<Objectdb::match> matches;
	vector<Objectdb::mynum> r;
	clock_t tstart=clock();
	db.findObjects( r,pts,origpts ? &matches : 0,classes, classes!=0 );
	out<<"time for analyzing a pointcloud with"
			<<pts.size()<<" points:"
			<<1000*(clock()-tstart)/CLOCKS_PER_SEC<<" ms"<<endl;

	if(origpts){

		mynum delt[3]={400,0,0};


		ofstream osg((db.getObjects()[oid].name+".osg").c_str());
		osg<<"Group{\n";


		osgpointcloud(osg,*origpts);

		osg<<
				"MatrixTransform {\n"
				;
		osgplaceandscalemat(osg,delt,1);
		osgpointcloud(osg,pts);
		osg<<
				"}\n"
				;


		matches.resize(
				std::remove_if(matches.begin(),matches.end(),
						toremove(oid)
				)-	matches.begin());

		sort(matches.begin(),matches.end(),matcheslt());


		//just that it's there and can be referenced
		osg<<"MatrixTransform { \n"
				"Matrix { 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 } \n";
		osgcircball(osg);
		osg<<"}\n";


		//for 2 best matches to orig do sth.:
		//draw balls and connector
		for(int i=0;i<matches.size() && i <2;++i)
		{
			Objectdb::match&m=matches[i];
			mynum center2[3];
			for(int i=0;i<3;++i)
				center2[i] = m.testmom.center[i]+delt[i];

			osgconnector(osg,m.dbmom.center,center2);

			osgplaceandscale(osg,m.dbmom.center,m.dbmom.radius);
			osgplaceandscale(osg,center2,m.dbmom.radius);
		}

		osg<<
				"}\n"
				;

	}


	vector<pair<Objectdb::mynum,int> > v;


	for(int i=0;i<r.size();++i)
		if(r[i]!=0)
			v.push_back(pair<Objectdb::mynum,int>(r[i],i));

	std::sort(v.begin(),v.end());

	for(int k=0;k<v.size();++k)
		out<<(v[k].second==oid?" * ":"   ")<<v[k].second<< ' '<<v[k].first<<endl;

	if(likelihoodoid)
	{
		*likelihoodoid = r[oid];
	}

	*rankoid=r.size();
	for(int i=0;i<v.size();++i)
		if(v[i].second==oid)
			*rankoid = v.size()-i;
	if(v.size()>0)
		return v.back().second;
	else
		return -1;
}







inline void writeHistPS(ostream&out,const vector<mynum> &vals,
		mynum w,mynum h,mynum sigmafract=4)
{

	Objectdb::stats s;
	for(int i=0;i<vals.size();++i)
		s.add(vals[i]);

	mynum avg=s.average();
	mynum sigma=sqrt(s.variance());
	mynum bs=sigmafract/sigma;
	int nr = int(ceil(bs*(s.max-avg)));
	int nl = int(ceil(bs*(avg-s.min)));
	int n=nl+nr;
	mynum off=avg-nl/bs;


	vector<mynum> npb(n,0);
	for(int i=0;i<vals.size();++i)
	{
		int ind=(int)floor( (vals[i]-off)*bs );
		if(ind>=n) {cout<<ind<<" > "<<n<<endl;ind=n-1;}
		if(ind<0){cout<<ind<<"<0"<<endl;ind=0;}
		++npb[ind];
	}

	mynum maxnp=0;
	for(int i=0;i<n;++i)
		if(npb[i]>maxnp)maxnp=npb[i];


	out
	<<"0 0 moveto "
	<<"0 "<<npb[0]*h/maxnp<<" lineto \n";

	off = (s.min-off)*bs;

	for(int i=1;i<n;++i)
	{
		if(npb[i]!=npb[i-1])
		{
			out
			<<int(i*w/n -off)<<' '<<int(npb[i-1]*h/maxnp)<<" lineto "
			<<int(i*w/n -off)<<' '<<int(npb[i]*h/maxnp)<<" lineto\n";
		}		
	}

	out<<w<<' '<<int(npb[n-1]*h/maxnp)<<" lineto "
			<<w<<" 0 lineto "
			<<int(nl*w/n-off)<<" 0 lineto "
			<<" 0 "<<h<<" rlineto "
			<<" 0 "<<-h<<" rlineto "
			<<w<<" 0 lineto\n";

}

inline void writeDistPS(ostream&out,vector<mynum> vals,
		mynum w,mynum h)
{

	sort(vals.begin(),vals.end());

	out<<"0 0 moveto ";
	mynum sc=w/(vals.back()-vals.front());

	int ox=0;
	int oy=0;

	out<<"/r { 1 rlineto } def\n";
	out<<"matrix currentmatrix .1 .1 scale\n";
	for(int i=0;i<vals.size();++i)
	{
		int x=int( 10*(vals[i]-vals.front())*sc),
				y=int( 10*(i+1)*h/(vals.size()));
		if(x!=ox && y!=oy){
			out<<(x-ox)<<" r\n";
			ox=x;oy=y;
		}
	}
	out<<"setmatrix\n";
	out<<w<<" 0 lineto 0 0 lineto\n";
}


inline void writeMomHistPS(ostream&out,const vector<mynum>& moms)
{

	vector<mynum> vals(moms.size()/tmom::n);

	int a4w = int(pow(2.0,-2.25)*72/0.0254);
	int a4h = int(pow(2.0,-1.75)*72/0.0254);

	//border offset
	int off = 40;

	//inter-line distance
	int linedist = 6; 



	int w=a4w-2*off,h=a4h-2*off;

	int dh=(h-(tmom::n-1)*linedist)/tmom::n;


	//save current transformation matrix
	out<<"gsave\n";
	out<<off<<" dup translate\n";
	//12 pt times font
	out<<"/Times-Roman findfont 12 scalefont setfont\n";

	for(int i=0;i<tmom::n;++i)
	{

		int k=i;
		for(int j=0;j<vals.size();++j,k+=tmom::n)
			vals[j]=moms[k];

		writeHistPS(out,vals,w,dh);
		out<<"0.6 setgray fill\n";
		writeHistPS(out,vals,w,dh);
		out<<"0 setgray stroke\n";
		writeDistPS(out,vals,w,dh);
		out<<"0 setgray stroke\n";


		out<<"2 2 moveto .4 setgray ("<<tmom::metainfo.me[i].des<<") show\n";


		out<<"0 "<<dh+linedist<<" translate\n";
	}


	//restore transformation matrix
	out<<"grestore\n"<<endl;;

}

template<class x>
std::string tostr(const x&y)
{
	std::ostringstream os;
	os<<y;
	return os.str();
}




inline void testWriteMomHistPS()
{

	static const int n=40000;

	vector<mynum> xx(n*tmom::n);

	for(int i=0;i<xx.size();++i)
		xx[i]=rand()+rand()+rand()+rand(), xx[i]*=xx[i];

	ofstream of("test.ps");

	writeMomHistPS(of,xx);

}






inline void orig2domord(vector<mynum> & origmoms)
{
	vector<mynum>::iterator it=origmoms.begin();
	while(it!=origmoms.end())
	{
		for(int i=0;i<tmom::n;++i,++it)
		{
			*it= (*it>0?1:-1) * pow(fabs(*it),1.0l/tmom::metainfo.me[i].dord);
		}
	}
}

inline void valord2orig(vector<mynum> & origmoms)
{
	vector<mynum>::iterator it=origmoms.begin();
	while(it!=origmoms.end())
	{
		for(int i=0;i<tmom::n;++i,++it)
		{
			*it=pow(*it,tmom::metainfo.me[i].vord);
		}
	}
}

inline void valord2tree(vector<mynum> & origmoms,const Objectdb::momPerRadius&m)
{
	vector<mynum>::iterator it=origmoms.begin();
	vector<mynum> buf(Objectdb::tmom::n);
	while(it!=origmoms.end())
	{
		m.tosearchable(&buf[0],&*it);
		it=copy(buf.begin(),buf.end(),it);
	}
}


inline void writeMomHistObjectdb(Objectdb&db)
{

	Objectdb::mapmom::const_iterator mit;

	for(mit=db.getMomsets().begin();mit!=db.getMomsets().end();++mit)
	{
		cout<<"Histograms for radius "<<mit->first<<endl;
		ofstream momhist;

		momhist.open(("histogram_r="+tostr(mit->first)+".ps").c_str());
		writeMomHistPS(momhist,mit->second.momentSets);
		momhist.close();

		momhist.open(("histogram_pca_r="+tostr(mit->first)+".ps").c_str());
		writeMomHistPS(momhist,mit->second.tree.getleaves());
		momhist.close();

		momhist.open(("histogram_original_r="+tostr(mit->first)+".ps").c_str());
		vector<mynum> origmoms(mit->second.momentSets);

		valord2orig(origmoms);
		writeMomHistPS(momhist,origmoms);
		momhist.close();

		orig2domord(origmoms);
		momhist.open(("histogram_domordersqr_r="+tostr(mit->first)+".ps").c_str());

		writeMomHistPS(momhist,origmoms);
		momhist.close();

	}




	//note that here the db should be fresh (loaded db's do not contain momentids for object)
	const vector<Objectdb::obj> &ob=db.getObjects();
	vector<mynum> objmoms,orig;

	for(int i=0;i<ob.size();++i)
	{
		cout<<"Histograms for object "<<i<<endl;
		Objectdb::obj::mapmoms::const_iterator mito;
		for(mito=ob[i].momentSets.begin();mito!=ob[i].momentSets.end();++mito)
		{

			const Objectdb::momPerRadius &mp
			=db.getMomsets().find( (mynum)mito->first )->second;

			const vector<int> & offs=mito->second;
			objmoms.resize(offs.size()*Objectdb::tmom::n);			
			const vector<mynum> & momper 
			= mp.momentSets;

			vector<mynum>::iterator oit=objmoms.begin();			
			//copy object's moments into contigous space
			for(int j=0;j<offs.size();++j)
			{
				oit = std::copy(momper.begin()+offs[j],
						momper.begin()+offs[j]+tmom::n, oit );
			}

			ofstream momhist;			
			momhist.open(("histogram_id="+tostr(i)+"r="+tostr(mito->first)+".ps").c_str());
			writeMomHistPS(momhist,objmoms);
			momhist.close();

			orig=objmoms;
			valord2orig(orig);

			momhist.open(("histogram_orig_id="+tostr(i)+"r="+tostr(mito->first)+".ps").c_str());
			writeMomHistPS(momhist,orig);
			momhist.close();

			orig2domord(orig);

			momhist.open(("histogram_domordsqr_id="+tostr(i)+"r="+tostr(mito->first)+".ps").c_str());
			writeMomHistPS(momhist,orig);
			momhist.close();

			orig=objmoms;			
			valord2tree(orig,mp);

			momhist.open(("histogram_pca_id="+tostr(i)+"r="+tostr(mito->first)+".ps").c_str());
			writeMomHistPS(momhist,orig);
			momhist.close();




		}	
	}
}	



enum datasetid
{DENKER , GAVAB, _3DRMA ,PRINCETON};

/** 
 *get the file names for the datasets and version:
 names[v][o] is the filename of object o with version v
 */
void getFileNamesOfPointClouds(vector<vector<string> > & names
		,datasetid type)
{
	switch(type)
	{

	case _3DRMA:
	{

		names.resize(6);

		char*vers[]={"1v1","1v2","1v3","2v1","2v2","2v3",0};

		for(int v=0;v<6;++v)
			getnames(names[v],"C:\\u\\max\\models\\3D_RMA\\all_auto\\",
					string(vers[v])+".xyz");

	}break;


	case GAVAB:
	{
		char*blavers[9]={"abajo","arriba","izquierda","frontal1","frontal2","derecha","gesto","risa","sonrisa"};
		char**vers=blavers;

		names.resize(9);

		for(int v=0;v<9;++v)
		{
			names[v].resize(61);

			for(int i=0;i<61;++i)
			{
				ostringstream os;
				os<<"C:\\u\\max\\models\\GAVAB\\allfaces\\"
						<<"cara"<<i+1<<'_'<<vers[v]<<".wrl"<<ends;
				names[v][i]=os.str();
			}
		}


	}break;

	case DENKER:default:
	{


		string prefix="c:\\u\\max\\models\\denker\\da\\";


		const char* nams[]={ 
				"dino",
				"huhn",
				"handy",
				"kaestchen2",
				"kaestchen",
				"steckdose",
				"test",
				"pig",
				0
		};

		names.resize(1);
		for(int i=0;nams[i];++i)
		{
			names[0].push_back(prefix+nams[i]+".dump");
		}

	}break;
	}

	/*
	//for test:
	names.resize(1);
	names[0].resize(1);
	 */
}




//1 row of object recognition statistics table
struct tablerow_creation
{
	int db_id;
	int version_in_db;
	int object_id;
	string objects_file_name;
	int number_of_points_original;
	float time; //time for computation in seconds

	void write(ostream&o)
	{
		o
		<<db_id
		<<' '<<version_in_db
		<<' '<<object_id
		<<' '<< ( objects_file_name.length()==0 ? "x" : objects_file_name)
		<<' '<<time
		<<' '<<number_of_points_original
		<<' ';
	}

	void read(istream&i)
	{
		i
		>>db_id
		>>version_in_db
		>>object_id
		>>objects_file_name
		>>time
		>>number_of_points_original
		;
	}


};


struct tablerow_recognition:public tablerow_creation
{

	int object_version_id;
	float noiselevel; //in ratio of std deviation radius
	float ratio_points_left; //in ratio of points left
	int rank_correct_object; //rank of corret object
	int first_ranked_id; //id of object which was first-ranked
	int number_of_points_test;

	void write(ostream&o)
	{
		tablerow_creation::write(o);

		o
		<<' '<<object_version_id
		<<' '<< noiselevel
		<<' '<< ratio_points_left
		<<' '<<rank_correct_object
		<<' '<<first_ranked_id
		<<' '<<number_of_points_test
		<<' ';
	}

	void read(istream&i)
	{
		tablerow_creation::read(i);

		i
		>>object_version_id
		>>noiselevel
		>>ratio_points_left
		>>rank_correct_object
		>>first_ranked_id
		>>number_of_points_test
		;
	}






};





//create statistics from a table containing tablerow_recognition data
struct sft
{


	/*
	map
	 db_id
	 x
	 version_in_db 
	 x 
	 object_version_id
	 x
	 noiselevel
	 x
	 ratio_points_left

	 to

		number_entries,
		num_ranked[],
		sum_time 
	 */


	struct int0
	{
		int i;
		int0(){i=0;}
	};

	struct val{
		int number_entries;
		float sum_time;

		map<int,int0> num_ranked;

		val(){
			number_entries=0;
			sum_time=0;
		}
	};


	typedef map<float, val > tratio;
	typedef map<float, tratio > tnoise;
	//to avoid compiler warning of too-long names
	struct trat { tnoise tr; };
	typedef map<int ,trat > tversion;
	typedef map<int,tversion > tversiondb;
	typedef map<int,tversiondb > tvalperdb;

	tvalperdb	data0;



	void addrow(tablerow_recognition&r)
	{
		val & v
		=data0[r.db_id][r.version_in_db][r.object_version_id]
										 .tr[r.noiselevel][r.ratio_points_left];

		v.number_entries ++;
		v.num_ranked[r.rank_correct_object].i ++;
		v.sum_time += r.time;	
	}


	void addFile(const string& f)
	{
		cout<<"adding file "<<f<<endl;
		ifstream i(f.c_str());
		tablerow_recognition r;		
		int it=0;
		for(r.read(i);!i.bad()&&!i.fail();r.read(i),++it){
			if((it%1000)==0)cout<<"\r"<<it<<flush;
			addrow(r);
		}
	}


	void addFiles()
	{
		vector<string>names;
		getnames(names,".\\","find_db_table.txt");
		for(int i=0;i<names.size();++i)
			addFile(names[i]);
	}




	void createforgnuplot()
	{

		tvalperdb::iterator i;
		tversiondb::iterator j;
		tversion::iterator k;
		tnoise::iterator l; 
		tratio::iterator m;


		ofstream cmds("createplots.gnuplot");
		for(i=data0.begin();i!=data0.end();++i)
			for(j=i->second.begin();j!=i->second.end();++j)
				for(k=j->second.begin();k!=j->second.end();++k)
				{
					ostringstream ijk;
					ijk<<i->first<<'_'<<j->first<<'_'<<k->first
							<<"_data.dat";

					cout<<"writing to "<<ijk.str()<<endl;
					ofstream out(ijk.str().c_str());

					for(l=k->second.tr.begin();l!=k->second.tr.end();++l)
					{
						for(m=l->second.begin();m!=l->second.end();++m)
						{					


							out<<l->first<<' '	
									<<m->first<<' ';


							out<< m->second.sum_time/m->second.number_entries<<' ';

							int tot=0;
							for(int ii=0;ii<4;++ii)
							{
								tot+=m->second.num_ranked[ii].i;
								out<< float(tot) / m->second.number_entries<<' ';
							}		
							out<<'\n';
						}	
						out<<'\n'<<'\n';
					}
				}
	}


	sft(){
		addFiles();
		createforgnuplot();
	}





};






struct pt{ 
	mynum x[3];	
	bool operator <(const pt& b) const
	{
		return x[0]<b.x[0];
	}
};

void sortPointCloud(vector<mynum> &ret)
{
	if(ret.size()<6)return;
	pt *a = (pt*)&ret[0];
	sort( a,a+ret.size()/3 );
}


//create a random rotation matrix
void createRandRot(mynum m[3][3])
{
	for(int j=0;j<2;++j)
	{
		for(int i=0;i<3;++i)
			m[j][i]=rand();
	}


	//make m0,m1 perpendicular
	//m1 -=  m0  (m1 m0) / (m0 m0)

	mynum m0m1=0,m0m0=0,invd;
	for(int i=0;i<3;++i) m0m1+= m[0][i]*m[1][i],m0m0+=m[0][i]*m[0][i];
	m0m1/=m0m0;
	for(int i=0;i<3;++i) m[1][i]-= m[0][i]*m0m1;

	//norm m0
	invd=1.0/sqrt(m0m0);
	for(int i=0;i<3;++i) m[0][i]*=invd;

	//norm m1
	invd=0;
	for(int i=0;i<3;++i) invd += m[1][i]*m[1][i];
	invd=1.0/sqrt(invd);
	for(int i=0;i<3;++i) m[1][i]*=invd;

	//m2 = m0 % m1
	for(int i=0;i<3;++i)
	{
		int j=(i+1)%3,k=(i+2)%3;
		m[2][k]=m[0][i]*m[1][j] -m[1][i]*m[0][j];
	}
}


void noisePointCloud(vector<mynum>& ret,const vector<mynum>& orig,mynum ratio)
{
	Objectdb::stats s[3];

	for(int i=0;i<orig.size();)
	{
		s[0].add(orig[i]);++i;
		s[1].add(orig[i]);++i;
		s[2].add(orig[i]);++i;	
	}

	ratio *= sqrt(s[0].variance()+s[1].variance()+s[2].variance());

	ret.resize(orig.size());
	for(int i=0;i<orig.size();++i)
	{
		ret[i] =  orig[i] + ratio * rand() / RAND_MAX;
	}

}







void rotatePointCloud(vector<mynum>& ret,const vector<mynum>& orig,const mynum m[][3])
{
	ret.resize(orig.size());
	vector<mynum>::iterator a;
	vector<mynum>::const_iterator b;
	//a is a rot.version of b.
	for(a=ret.begin(),b=orig.begin();a<ret.end();a+=3,b+=3)
	{
		for(int i=0;i<3;++i)
		{
			a[i]=0;
			for(int j=0;j<3;++j)
				a[i]+=m[i][j] * b[j];
		}
	}
}

#ifndef min
#define min(a,b) (a<b?a:b)
#endif

void createTable()
{

	vector<mynum> pts;
	tablerow_recognition r;



	for(r.db_id=1;r.db_id<3;++r.db_id)
	{
		vector< vector<string> > names;		
		getFileNamesOfPointClouds(names,(datasetid)r.db_id);

		for(r.version_in_db=3;r.version_in_db<min(5,names.size());++r.version_in_db)
		{

			ostringstream a;
			a<<r.db_id<<'_'<<r.version_in_db<<'_';
			ofstream creation_table((a.str()+"create_db_table.txt").c_str());
			ofstream finding_table((a.str()+"find_db_table.txt").c_str());
			ofstream logfile((a.str()+"create_big_table.log").c_str());

			Objectdb db;
			db.dbgout=&logfile;

			for(r.object_id=0;r.object_id<names[r.version_in_db].size();++r.object_id)
			{				
				r.objects_file_name=names[r.version_in_db][r.object_id];				
				PointCloudIO::read(pts,r.objects_file_name.c_str());				
				r.number_of_points_original = pts.size()/3;				
				clock_t start=clock();				
				db.addObject(r.objects_file_name,pts);				
				r.time = 1.0 * (clock()-start) /CLOCKS_PER_SEC;

				((tablerow_creation&)r).write(creation_table);creation_table<<endl;			
				((tablerow_creation&)r).write(cout);cout<<endl;			

			}

			for(r.object_version_id=3;r.object_version_id<min(5,names.size());++r.object_version_id)
			{				
				for(r.object_id=0;r.object_id<names[r.object_version_id].size();++r.object_id)
				{
					r.objects_file_name=names[r.object_version_id][r.object_id];				
					PointCloudIO::read(pts,r.objects_file_name.c_str());
					r.number_of_points_original=pts.size()/3;

					vector<mynum> testpts;

					for(r.noiselevel=0;r.noiselevel<=0.1;r.noiselevel=r.noiselevel*2+0.01)
					{
						mynum mat[3][3];
						createRandRot(mat);

						rotatePointCloud(testpts,pts,mat);

						noisePointCloud(testpts,testpts,r.noiselevel);
						sortPointCloud(testpts);										

						for(r.ratio_points_left=1;r.ratio_points_left>0.4;r.ratio_points_left-=0.15f)
						{
							r.number_of_points_test 
							=(int) (r.number_of_points_original*r.ratio_points_left);

							testpts.resize(r.number_of_points_test*3);




							vector<mynum> probs;

							clock_t start=clock();			
							db.findObjects(probs,testpts);
							r.time = 1.0 * (clock()-start) /CLOCKS_PER_SEC;

							vector<pair<mynum,int> > probs2;
							for(int i=0;i<probs.size();++i)
								if(probs[i])probs2.push_back(pair<mynum,int>(probs[i],i));

							sort(probs2.begin(),probs2.end());

							r.first_ranked_id = probs2.empty() ? -1 : probs2.back().second;
							r.rank_correct_object=probs.size()-1;
							for(int i=probs2.size()-1;i>=0;--i)
								if(probs2[i].second == r.object_id)
								{
									r.rank_correct_object = probs2.size()-1-i;
									break;
								}										
							r.write(cout);cout<<endl;
							r.write(finding_table);
							finding_table<<endl;
						}
					}
				}
			}
		}
	}	
}



void createascs(float ptleft,float noise,const string& fname)
{

	vector<mynum> pts;
	ostringstream outfile;outfile<<fname<<'_'<<ptleft<<'_'<<noise<<".asc";

	PointCloudIO::read(pts,fname.c_str());
	int n=pts.size()/3;

	vector<mynum> testpts;
	mynum mat[3][3];
	createRandRot(mat);						
	rotatePointCloud(testpts,pts,mat);						
	noisePointCloud(testpts,testpts,noise);
	sortPointCloud(testpts);										

	int n2=(int)( ptleft * n);
	testpts.resize(n2*3);

	PointCloudIO::write(testpts,outfile.str().c_str());


}


inline void writeObjMomsCsv(Objectdb&db)
{

	Objectdb::mapmom::const_iterator it;
	const Objectdb::mapmom &m=db.getMomsets();

	int i=0;
	for(it=m.begin();it!=m.end();++i,++it)
	{
		ostringstream oss;
		oss<<"facemoms\\facemoms_r="<<it->first<<".csv";
		PointCloudIO::write(it->second.momentSets,oss.str().c_str(),"csv28");	
	}




}




inline string intermediatedbname (const string&dbname,int nameid)
{
	ostringstream os;
	os<<"Objectdb"<<dbname<<"_intermediate"<<nameid<<".raw"<<ends;
	return os.str();
}



inline void testObjectdb(ostream &out)
{
	outputflow::dbgout=&out;

	outputflow ofl("testObjectdb");

	vector<mynum> pts;

	Objectdb db;

	char buf[256];



	datasetid datb=DENKER;
	/** for db reduction: which distance (in std.dev units) to be used as "equal" */
	mynum REDUCE_DEV=0;

#ifdef USE_SVM
	static const bool svm=true;
#else
	static const bool svm=false;
#endif


	char**vers=0;
	int onlyversion;
	vector<string> names;
	bool testmode=false;

	bool resampleVrml=true;
#ifdef USE_PSB
	PSBCategories dbclasses,testclasses;
#endif

	switch(datb)
	{
	case GAVAB:
	{
		char*blavers[]={"abajo","arriba","izquierda","frontal1","frontal2","derecha","gesto","risa","sonrisa",0};
		vers=blavers;

		//only two objects to test program correctness ?
		onlyversion=3;

		for(int i=0;i<(testmode? 2:61);++i)
			for(int j= (onlyversion<0? 0 :onlyversion) ;
					onlyversion>=0 ? j==onlyversion : vers[j]!=0 ;++j)
			{
				sprintf(buf,
						"C:\\u\\max\\models\\GAVAB\\allfaces\\"
						"cara%d_%.60s",i+1,vers[j]);
				names.push_back(buf);
				out<<names.back()<<endl;
			}
	}break;
	case _3DRMA:
	{

		onlyversion=1;

		char*blavers[]={"1v1","1v2","1v3","2v1","2v2","2v3",0};
		vers=blavers;

		getnames(names,"C:\\u\\max\\models\\3D_RMA\\all_auto\\",
				onlyversion>=0 ? string(vers[onlyversion]) : "" 
		);

		if(testmode&&names.size()>6)names.resize(6);

		cout<<"db candidates:"<<endl;
		for(int i=0;i<names.size();++i)
			cout<<names[i]<<endl;

	}break;
#ifdef USE_PSB
	case PRINCETON:
	{
		getPrincetonClasses(dbclasses,testclasses);

		const vector<int>&ids=dbclasses.mods();
		for(vector<int>::const_iterator cit=ids.begin();cit!=ids.end();++cit)
		{
			names.push_back(PSBCategories::getModelName(*cit));

			assert(*cit==PSBCategories::idForModelName(names.back()));
		}
		break;

	}break;
#endif

	case DENKER:
	{


		string prefix="c:\\u\\max\\models\\denker\\da\\";


		const char* nams[]={ 
				"dino",
				"huhn",
				"handy",
				"kaestchen2",
				"kaestchen",
				"steckdose",
				"test",
				"pig",
				0
		};


		for(const char**nam=nams;*nam!=0;++nam)
		{
			names.push_back(prefix+*nam);
		}

	}break;




	}

	string dbpostfix;

	switch(datb)
	{
	case GAVAB:dbpostfix="_gavab";
	break;
	case _3DRMA:dbpostfix="_3drma";
	break;
	case PRINCETON:dbpostfix="_prince";
	break;
	case DENKER:dbpostfix="_denker";
	break;
	default: dbpostfix="";
	}

	string dbname="Objectdb"+dbpostfix+".raw";

	char yesno;
	cout<<"Create database(y/n)"<<flush;
	cin>>yesno;

	if(yesno=='y')
	{

		int odbsiz=0;

		int nameid=0;
		string intermeddb;
		for(nameid=names.size()-1;nameid>=0;--nameid)
		{
			intermeddb= intermediatedbname(dbpostfix,nameid);
			struct ::stat stt;
			errno=0;
			int retval= ::stat(intermeddb.c_str(),&stt);
			int ee=errno;
			if( retval==0 )
				break;
		}

		if(nameid>=0)
		{
			cout<<"load intermediate db"<<intermeddb<<" (y/n)"<<flush;
			cin>>yesno;

			if(yesno!='y')
			{
				nameid=0;
			}
			else
			{
				db.load(intermeddb.c_str());
				++nameid;
			}
		}
		else
		{
			nameid=0;
		}

		for(;nameid<names.size();++nameid)
		{		
			outputflow xxxx("loadAddObject",(string("name=\"")+names[nameid]+"\"").c_str());

			cout<<nameid+1<<'/'<<names.size()<<':'<<names[nameid]<<"...\r"<<flush;
			pts.resize(0);
			out<<"trying to read"<<names[nameid]<<endl;
			try{
				PointCloudIO::read(pts,(names[nameid]+".rw2").c_str());
			}catch(...){
				out<<"no rw2 for "+names[nameid]<<" present yet"<<endl;
			}

			if(pts.size()==0)
			{
				switch(datb)
				{
				case GAVAB:
					PointCloudIO::read(pts,(names[nameid]+".wrl").c_str(),resampleVrml? "equaldist":0);
					break;

				case PRINCETON:
					PointCloudIO::read(pts,(names[nameid]+".off").c_str());
					break;
				case _3DRMA:
					PointCloudIO::read(pts,(names[nameid]+".xyz").c_str());
					break;
				case DENKER:
					PointCloudIO::read(pts,(names[nameid]+".dump").c_str());
					break;
				default:
					break;
				}

				PointCloudIO::write(pts,(names[nameid]+".rw2").c_str());
			}

			clock_t tstart=clock();

			string commonname=names[nameid];

			switch(datb){
			case GAVAB:
				commonname=commonname.substr(0,commonname.find('_'));break;
			case _3DRMA:
				commonname=commonname.substr(0,commonname.length()-7);
				break;
			default:
				break;
			}


			out<<"commonname:"<<commonname;


			const Objectdb::obj &ob=db.addObject(commonname,pts,false);		

			/*

			Objectdb::obj::mapmoms::const_iterator mit;

			for(mit=ob.momentSets.begin();mit!=ob.momentSets.end();++mit)
			{

				vector<int>::const_iterator 
					imom=mit->second.begin(),
					imomend=mit->second.end();

				vector<mynum> centers(mit->second.size()*3);
				vector<mynum>::iterator cit=centers.begin();
				for(;imom!=imomend;++imom)
				{
					const Objectdb::mynum* c=
						db.getMomsets().find(mit->first)->second.meta[*imom/Objectdb::tmom::n].center;
					cit=copy(c,c+3,cit);
				}

				sprintf(buf,"%.60s%.10s_centers_r=%lg.osg",
					commonname.c_str(),
					(onlyversion>=0&&vers!=0 ? vers[onlyversion]:""),
					mit->first);
				PointCloudIO::write(centers,buf);


			}
			 */

			out<<" time for adding a "<<pts.size()/3000<<" k object to the database:"<<
					float(clock()-tstart)/CLOCKS_PER_SEC<<" seconds"<<endl;


			int ndbsiz=db.nmoments();
			out<<"current db size:"
					<<ndbsiz<<" moments"<<endl;

			if(ndbsiz > 1.2*odbsiz+50000)
			{
				tstart=clock();

				if(REDUCE_DEV)
					db.reduce(REDUCE_DEV);

				odbsiz=db.nmoments();
				out<<"time for db reduction:"<<float(clock()-tstart)/CLOCKS_PER_SEC<<" seconds"<<endl;
				out<<"new db size:"
						<<odbsiz<<" moments"<<endl;

			}


			if( (nameid %10)==0 )
			{	
				db.save(intermediatedbname(dbpostfix,nameid).c_str());
			}

		}

		if(REDUCE_DEV)
			db.reduce(REDUCE_DEV);

		db.save(dbname.c_str());
	}
	else
		db.load(dbname.c_str());


	//writeObjMomsCsv(db);return;


	//out<<db.getMeta().size()<<" moments before reduce"<<endl;
	//db.reduce(6);
	//out<<db.getMeta().size()<<" moments after reduce"<<endl;

	/*
	cout<<"building kdtree of the moments..."<<flush;

	db.getTree();
	cout<<"ready."<<endl;
	 */
	mynum rmin,rmax;
	db.getradii(rmin,rmax);

	out<<"radii in tree:"<<rmin<<" - "<<rmax<<endl;

#ifdef USE_SVM
	Objectdbsvm odsvm;
	ostringstream ost;
	ost<<"Objectdb"<<dbpostfix<<".svm"<<ends;
	bool loaded=odsvm.load(ost.str().c_str());

	if(loaded)
	{
		cout<<"successfully loaded svm model from"<<ost.str()
								<<".\n still re-create svm model(y/n) ?"<<endl;
		char c;
		cin>>c;
		if(c=='y'){
			loaded=false;
		}
	}

	if(!loaded)
	{
		odsvm.trainModel(db);
		odsvm.save(ost.str().c_str());
	}
#else
	db.buildTrees();
#endif

	//writeMomHistObjectdb(db);
	//	out<<"ready. number of moments in the database:"<<db.objectofmoment.size()<<endl;

	//now try(with all versions)
	vector<int> rankatcut(names.size(),0),rankatrot(names.size(),0),
			rankatnoise[6];

	for(int i=0;i<6;++i)
		rankatnoise[i].resize(names.size(),0);

	if(datb!=DENKER)
		names.clear();


	if(vers)
	{
		cout<<"which version of the faces should be used for testing:"<<endl;
		int j=0;
		for(;vers[j];++j)
			cout<<j<<':'<<vers[j]<<endl;
		cout<<"select one version[default="<<onlyversion<<"]:"<<flush;
		do{
			cin>>buf;
			if(strlen(buf)!=0)
				sscanf(buf,"%d",&onlyversion);
		}while(onlyversion<0 || onlyversion>=j);
	}


	vector<int>*classes=0;

	switch(datb)
	{

	case GAVAB:
		for(int i=0;i<61;++i){
			for(int j=onlyversion;j==onlyversion&&vers[j];++j)
			{
				sprintf(buf,
						"C:\\u\\max\\models\\GAVAB\\allfaces\\"
						"cara%d_%.60s",i+1,vers[j]);
				names.push_back(buf);
				out<<"will test:"<<names.back()<<endl;
			}
		}
		break;

#ifdef USE_PSB
	case PRINCETON:
	{
		const vector<int>&ids=testclasses.mods();
		for(int i=0;i<100;++i)
		{
			names.push_back(PSBCategories::getModelName(ids[ (ids.size()-1) *rand()/RAND_MAX]));
		}

		const vector<Objectdb::obj> &obj=db.getObjects();
		classes=new vector<int>(obj.size(),-1);

		for(int i=0;i<obj.size();++i)
			(*classes)[i] =  dbclasses.catIdForModelName(obj[i].name);
	}
	break;
#endif
	case _3DRMA:


		getnames(names,"C:\\u\\max\\models\\3D_RMA\\all_auto\\",
				onlyversion>=0 ? string(vers[onlyversion]) : ""
		);

		if(testmode&&names.size()>6)
			names.resize(6);
		break;

	}




	for(int nameid=0;nameid<names.size();++nameid)
	{	
		cerr<<
				"\r#################################################################"
				"\r testing "<<nameid+1<<" / " <<names.size()<<" : "<<flush;
		out
		<<
		"\n-----------------------------------------------------------------------\n"
		<<names[nameid]<<
		"\n\n"
		<<endl;

#ifdef USE_PSB
		if(classes)
		{
			out<<"expected category "<< testclasses.catForModelName(names[nameid]).fullName<<endl;
		}
#endif

		try{
			PointCloudIO::read(pts,(names[nameid]+".rw2").c_str());
		}
		catch(...){
			switch(datb){
			case GAVAB:
				PointCloudIO::read(pts,(names[nameid]+".wrl").c_str(),resampleVrml?"equaldist":0);
				break;
#ifdef USE_PSB
			case PRINCETON:
				PointCloudIO::read(pts,(names[nameid]+".off").c_str());
				break;
#endif
			case _3DRMA:
				PointCloudIO::read(pts,(names[nameid]+".xyz").c_str());
				break;
			case DENKER:
				PointCloudIO::read(pts,(names[nameid]+".dump").c_str());
				break;
			default:
				PointCloudIO::read(pts,names[nameid].c_str());
			}	
			PointCloudIO::write(pts,(names[nameid]+".rw2").c_str());
		}


		vector<mynum> opt=pts;


		int oid=-1;

		const vector<Objectdb::obj> & v=db.getObjects();

		mynum maxrad;
		if(classes &&datb==PRINCETON)
		{
#ifdef USE_PSB
			oid = 
					dbclasses.classWithFullName(testclasses.catForModelName(names[nameid]).fullName);
#endif
		}
		else
		{
			for(oid=0;oid<v.size();++oid)
				if(v[oid].name==names[nameid].substr(0,names[nameid].find('_')))break;

			if(oid==v.size())
			{
				for(oid=0;oid<v.size();++oid)
					if(names[nameid].substr(0,v[oid].name.size())==v[oid].name)
						break;
			}

			if(oid==v.size())
				continue;

			if(v[oid].momentSets.empty())
				continue;
		}

		cerr<<" noise:"<<flush;

		int nid=0;



		if(classes==0)
		{
			maxrad=(-- v[oid].momentSets.end())->first;
		}
		else
		{

			Objectdb::stats ms[3];
			for(int i=0;i<pts.size();i+=3)
			{
				ms[0].add(pts[i]);
				ms[1].add(pts[i+1]);
				ms[2].add(pts[i+2]);
			}

			maxrad=sqrt(ms[0].variance()+ms[1].variance()+ms[2].variance());
		}



		for(mynum noise=0;noise <=0 /* maxrad/16*/;noise=noise*2+maxrad/64,++nid){




			for(int j=0;j<pts.size();++j)
				pts[j]+= noise*((2.0*rand()/RAND_MAX)-1);

			cerr<<noise<<' '<<flush;
			out<<"recognition at noise level "<<noise<<":"<<endl;
			//sprintf(buf,"%.60s_noiselevel%d.osg",names[nameid].c_str(),(int)noise);

			//PointCloudIO::write(pts,buf);



			int rank=1;

			vector<double> probs;

#ifdef USE_SVM
			int found=odsvm.predict(probs,pts);
#else
			int found=testfind(db,pts,oid,out,0,&rank,0,classes);
#endif
			if(classes )
			{
				#ifdef USE_PSB
				out<<"expected:"<< (oid>=0?dbclasses.cats()[oid].fullName:"nothing")<<endl;
				out<<"found:"<< (found>=0? dbclasses.cats()[found].fullName: "nothing")<<endl;
				#else
				out<<"expected:"<< (oid>=0?db.getObjects()[oid].name:"nothing")<<endl;
				out<<"found:"<< (found>=0? db.getObjects()[found].name: "nothing")<<endl;
				#endif
			}
			else
			{
				out<<"expected:"<< (oid>=0?db.getObjects()[oid].name:"nothing")<<endl;
				out<<"found:"<< (found>=0? db.getObjects()[found].name: "nothing")<<endl;			
			}

			if(nid<6){
				rankatnoise[nid][rank-1]++;
			}

			out<<endl;
		}

		//taint the cloud by arbitrary cuts:
#ifdef ROTATION

		cerr<< " rotation " << flush;
		out<<"\nnow testing rotations"<<endl;
		for(int ii=0;ii<1;++ii)
		{

			mynum m[3][3];

			createRandRot(m);

			out<<"rotation matrix:"<<endl;
			{for(int i=0;i<3;++i)
				for(int j=0;j<3;++j)
					out<<m[i][j]<<(j==2?'\n':' ')<<flush;}


			rotatePointCloud(pts,opt,m);


			out<<"ii="<<ii<<':'<<endl;
			sprintf(buf,"%s_rot%d.osg",names[nameid].c_str(),ii);
			PointCloudIO::write(pts,buf);

			//test recognition

			int rank;
			int found=testfind(db,pts,oid,out,0,&rank,&opt,classes);
			rankatrot[rank-1]++;


		}
#endif
#ifdef CUTS
		cerr<<"cuts:"<<flush;
		out<<"\nnow testing cuts:"<<endl;
		for(int ii=0;ii<3;++ii)
		{
			cerr<<ii<<flush;
			pts=opt;
			int npts=opt.size()/3;

			mynum dir[]={rand()-RAND_MAX/2,
					rand()-RAND_MAX/2,rand()-RAND_MAX/2};

			vector<mynum>::iterator it;
			mynum mid[3]={0,0,0};
			for(it=pts.begin();it!=pts.end();it+=3)
				mid[0]+=it[0],mid[1]+=it[1],mid[2]+=it[2];


			mid[0]/=npts,mid[1]/=npts,mid[2]/=npts;

			int npt=npts;
			for(it=pts.begin();it!=pts.end();it+=3)
			{
				mynum dec=0;
				for(int j=0;j<3;++j)
					dec+=(it[j]-mid[j])*dir[j];
				//set all pts left dec to 0
				if(dec<0)
				{
					for(int j=0;j<3;++j)
						it[j]=0;
					--npt;
				}


			}

			out<<npt <<" of "<<opt.size()/3<<"points left:"<<endl;

			out<<"ii="<<ii<<':'<<endl;
			sprintf(buf,"%s_cut%d.osg",names[nameid].c_str(),ii);
			PointCloudIO::write(pts,buf);

			mynum likelihoodoid;int rankoid;
			int found=testfind(db,pts,oid,out,0,&rankoid,0,classes);
			rankatcut[rankoid-1]++;


		}

#endif

	}


	int nl=0;
	for(int i=0;i<4;++i,nl = nl*2+1)
	{
		out<<"noise level "<<nl<<":\n";
		outputranks(rankatnoise[i],out);
	}
	/*
	out<<"cuts:\n";
	outputranks(rankatcut,out);
	out<<"rotations:\n";
	outputranks(rankatrot,out);
	 */


}




//create a lipschitz-continuous surface which is repetitive
//(can be used as tile)
inline void createRandFractSurf(int pn, mynum *  surf)
{
	int n= (1<<pn), andn = n-1;		

	mynum*sur=surf;
	for(int y=0;y<n;++y)
		for(int x=0;x<n;++x,sur+=3)
		{
			sur[0]=x;
			sur[1]=y;
		}

	surf[2]=0;

	for(int d=1<<pn ; d>1;d>>=1)
	{
		int d2=d>>1;


		//compute z values with respect to left and right neighbors

		sur=surf+2; 	
		for(int y=0;y<n;y+=d,sur+=d*3*n)
		{
			mynum zr = sur[ 0 ],zl;
			for(int x=d2;x<n;x+=d)
			{
				zl=zr; zr = sur[ ((x + d2)&andn)*3 ];
				mynum r = (d - fabs(zr-zl)) ;	
				sur[ x*3 ]
					 = 0.5*( zr + zl - r )  + r * rand() / RAND_MAX; 				
			}
		}

		//compute z values with respect to top and bottom neighbors

		sur = surf+2;
		for(int x=0;x<n;x+=d2,sur+=d2*3)
		{
			mynum zr = sur[ 0 ],zl;
			for(int y=d2;y<n;y+=d)
			{
				zl = zr; zr = sur[ ((y + d2)&andn)*3*n];
				mynum r = (d - fabs(zr-zl)) ;	
				sur[ y*3*n ]
					 = 0.5*( zr + zl - r )  + r * rand() / RAND_MAX; 				
			}
		}
	}
}


//non-repetitive, grid is (1<<pn)+1 x (1<<pn)+1
inline void createRandFractSurf2(int pn, mynum *  surf)
{
	int n= (1<<pn)+1 ;

	mynum*sur=surf;
	for(int y=0;y<n;++y)
		for(int x=0;x<n;++x,sur+=3)
		{
			sur[0]=x;
			sur[1]=y;
		}

	surf[2] =rand()*n / RAND_MAX;

	surf[ 3*n -1 ]= rand()*n / RAND_MAX;

	surf[ 3*n*(n-1) +2 ]= rand()*n / RAND_MAX;

	surf[ 3*n*n -1 ]= rand()*n / RAND_MAX;

	for(int d=1<<pn ; d>1;d>>=1)
	{
		int d2=d>>1;


		//compute z values with respect to left and right neighbors

		sur=surf+2; 	
		for(int y=0;y<n;y+=d,sur+=d*3*n)
		{
			mynum zr = sur[ 0 ],zl;
			for(int x=d2;x<n;x+=d)
			{
				zl=zr; zr = sur[ (x + d2)*3 ];
				mynum r = (d - fabs(zr-zl)) ;	
				sur[ x*3 ]
					 = 0.5*( zr + zl - r )  + r * rand() / RAND_MAX; 				
			}
		}

		//compute z values with respect to top and bottom neighbors

		sur = surf+2;
		for(int x=0;x<n;x+=d2,sur+=d2*3)
		{
			mynum zr = sur[ 0 ],zl;
			for(int y=d2;y<n;y+=d)
			{
				zl = zr; zr = sur[ (y + d2)*3*n];
				mynum r = (d - fabs(zr-zl)) ;	
				sur[ y*3*n ]
					 = 0.5*( zr + zl - r )  + r * rand() / RAND_MAX; 				
			}
		}
	}
}


inline void kdtreeEntryFromTset(
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
				mynum(1)/faceMomSet::metainfo.me[j].vord);

		*(ret++)=v;
	}
}




inline void makecovrandsurf()
{
	mynum f[9*9*3];

	tensorSet3D<> ts;

	static const int nn=faceMomSet::n;

	vecstatscollector<nn> coll;

	faceMomSet bufm;
	mynum moms[nn];

	TNT::Array2D<mynum> ev;
	TNT::Array1D<mynum> ew;

	printf("%15s %15s %15s %15s\n",
			"ew[nn-1]","ev[nn-1][nn-1]","ew[10]","ev[10][10]"
	);

	ofstream out("bla.asc");

	for(int it=1;it<1e7;++it)
	{
		createRandFractSurf2(3,f);

		if(it<20)
			for(int i=0;i<9*9*3;++i)
			{
				out<<f[i]
					   +(i%3==1?it*10:0)
					   +(i%3==0? (it&1)*3 : 0 )
					   <<(i%3==2?'\n':' ');
			}
		else if(it==20)
			out.close();

		/*
		for(int y=0;y<8;++y){
			for(int x=0;x<8;++x)
				cout<<char('a' + char(f[(y*8+x)*3+2]+8) );
			cout<<endl;
		}
		 */


		for(int i=0;i<9*9;++i)
			ts.add3dcoords(f+i*3);

		ts.transinv();
		kdtreeEntryFromTset(moms,bufm,ts,4);

		coll.add(moms);


		if(it%100000==10)
		{

			vecstats<nn> vs = coll.getstatscov();

			TNT::Array2D<mynum> a2d(nn,nn,&vs.cov[0][0]);
			JAMA::Eigenvalue<mynum> eigenv(a2d);

			eigenv.getRealEigenvalues(ew);
			eigenv.getV(ev);


			printf("%+15.6e %+15.6e %+15.6e %+15.6e\n",
					ew[nn-1],ev[nn-1][nn-1],ew[10],ev[10][10]);


		}

	}

	cout<<"\neigenvalues of co.mat."<<endl;// of "<<it<<" surfaces:"<<endl;
	for(int i=0;i<ew.dim();++i)
		cout<<i<<':'<<ew[i]<<' '<<flush;
	cout<<endl;

	cout<<"eigenvector marix:\n";
	for(int i=0;i<nn;++i){
		for(int j=0;j<nn;++j)
			cout<<ev[i][j] <<' ';
		cout<<'\n';
	}
	cout<<endl;

	mynum convert[nn][nn];
	//mynum m=sqrt(1.0/ew[tmom::n-1]);
	for(int i=0;i<nn;++i){
		mynum m = sqrt(1.0/fabs(ew[i]));
		for(int j=0;j<nn;++j)
			convert[i][j] = ev[j][i]*m;
	}

	cout<<"conversion marix:\n";
	for(int i=0;i<nn;++i){
		for(int j=0;j<nn;++j)
			cout<<convert[i][j] <<' ';
		cout<<'\n';
	}
	cout<<endl;

}












inline void convert2raw(int argc, const char*argv[])
{

	if(argc>2)
	{
		//convert file to another format
		char buf[256];
		vector<mynum> pts;

		sprintf(buf,"%s",argv[1]);
		cout<<"reading in "<<buf<<" ... "<<flush;
		PointCloudIO::read(pts,buf);

		sprintf(buf,"%s",argv[2]);
		cout<<"ready.\n"
				"writing binary long double values to "
				<<buf<<"... "<<flush;
		PointCloudIO::write(pts,buf);
		cout<<"ready"<<endl;
	}
}


inline void testObjectdbAndLog()
{
	char buf[80];
	sprintf(buf,"testobjectdb.log");
	struct stat st;
	for(int i=0;stat(buf,&st) != -1; ++i)
		sprintf(buf,"testobjectdb_%d.log",i);		

	ofstream out(buf);

	testObjectdb(out);

	out.close();


}






inline mynum testf(const mynum*p)
{
	return p[0]*p[0] +p[1]*p[2] -p[2]*p[2] ;
	//+p[1]*p[1]*p[2]
	//+p[2]*p[2]*p[0];
}

struct tf
{
	template<class T>
	T operator()(const T &x)const
	{
		//return (x-4)*(x-3)*(x-100000)*x*x*x ;

		T xpi= x * M_PI;
		return xpi*sin(xpi)+xpi;

		//T xx=x*x;
		//return (((xx+1)*xx+2)*xx +3)*xx -1;
	}

};

struct tdf
{
	template<class T>
	T operator()(const T &x)const
	{

		//return x*x*(((x*6-500035)*x+2800048)*x-3600000) ;
		T xpi = x*M_PI;
		return (xpi*cos(xpi) +sin(xpi) + 1)*M_PI;
		//T xx=x*x;
		//		return (((xx*8 +6)*xx +8)*xx +6)*x; 
	}
};



inline void testsolve()
{

	cout.precision(30);
	typedef mynum T;

	//if using gnu g++, replace it with 
	//compiler switch -mfp-rounding-mode=d
#ifndef GNUCC
	_controlfp(_RC_DOWN,_MCW_RC);
#endif

	vector<myinterval<T> >ret;
	myinterval<T> s(-60,60);

	/*

	solve(ret,s,tf(),tdf());

  cout<<"output solve:"<<endl;
	for(int i=0;i<ret.size();++i)
		cout<<ret[i]<<' ';
	cout<<endl;

	solve2(ret,s,tf(),tdf());

  cout<<"output solve2:"<<endl;
	for(int i=0;i<ret.size();++i)
		cout<<ret[i]<<' ';
	cout<<endl;

	solve3< T,tf, tdf > solv(ret,s,1e-16);

  cout<<"output solve3:"<<endl;
	for(int i=0;i<ret.size();++i)
		cout<<ret[i]<<' '<<ret[i].width()<<' '<<tf()(ret[i])<<' '
		<<tf()(ret[i]).width()<<endl;
	cout<<endl;
	 */
	solve4< T,tf, tdf > solv4(ret,s,-1e-10);

	cout<<"output solve4:"<<endl;
	for(int i=0;i<ret.size();++i)
		cout<<ret[i]<<' '<<ret[i].width()<<' '<<tf()(ret[i])<<' '
		<<tf()(ret[i]).width()<<endl;
	cout<<endl;
	/*

	float start,diff;int num;int trynum=1000;

	start=myclock();num=0;
	while( (diff=myclock()-start)< 2 )
	{
		for(int i=0;i<trynum;++i)
			solve<T>(ret,s,tf(),tdf(),1e-6);

		num+=trynum;
	}
	cout<<" time per computation solve: "<<diff*1e6/1/num
		<<" mus"<<endl;

	start=myclock();num=0;
	while( (diff=myclock()-start)< 2 )
	{
		for(int i=0;i<trynum;++i)
			solve2<T>(ret,s,tf(),tdf(),1e-6);

		num+=trynum;
	}
	cout<<" time per computation solve2: "<<diff*1e6/1/num
		<<" mus"<<endl;

	start=myclock();num=0;
	while( (diff=myclock()-start)< 2 )
	{
		for(int i=0;i<trynum;++i)
			solv.solve(ret,s,1e-6);

		num+=trynum;
	}
	cout<<" time per computation solve3: "<<diff*1e6/1/num
		<<" mus"<<endl;

	start=myclock();num=0;
	while( (diff=myclock()-start)< 2 )
	{
		for(int i=0;i<trynum;++i)
			solv4.solve(ret,s,1e-6);

		num+=trynum;
	}
	cout<<" time per computation solve4: "<<diff*1e6/1/num
		<<" mus"<<endl;
	 */
}






//write the functions needed for a local-support symmetric weighting functions
//with smooth cut-off of continuity order n





inline void tounitcov
(vector<mynum>&momsets2,const vector<mynum>&momentSets,int&treedim)
{

	vector<mynum>::const_iterator momit;
	vecstatscollector<tmom::n,mynum> vstats;

	for(momit=momentSets.begin();momit!=momentSets.end();momit+=tmom::n)
		vstats.add(&*momit);

	vecstats<tmom::n,mynum> vs	= vstats.getstatscov();

	TNT::Array2D<mynum> a2d(tmom::n,tmom::n,&vs.cov[0][0]);
	JAMA::Eigenvalue<mynum> eigenv(a2d);

	TNT::Array2D<mynum> ev;
	TNT::Array1D<mynum> ew;

	eigenv.getRealEigenvalues(ew);
	eigenv.getV(ev);

	treedim=tmom::n;
	while( treedim<tmom::n && 
			ew[tmom::n-1- treedim ] > ew[tmom::n-1]*1e-3 )
		++treedim;

	mynum convert[tmom::n][tmom::n];

	//mynum m=sqrt(1.0/ew[tmom::n-1]);
	for(int i=0;i<treedim;++i){
		int revi = tmom::n-1-i;
		mynum m = sqrt(1.0/ew[revi]);
		for(int j=0;j<tmom::n;++j)
			convert[i][j] = ev[j][revi]*m;
	}


	int nmoms=(momentSets.size()/tmom::n);
	momsets2.resize(treedim * nmoms ) ;	

	//multiply with conversion matrix

	vector<mynum>::iterator momit2=momsets2.begin();
	for(momit=momentSets.begin(); momit!=momentSets.end(); 
			momit+=tmom::n)
	{
		for(int i=0;i<treedim;++i,++momit2)
		{
			*momit2=0;

			for(int j=0;j<tmom::n;++j)
				*momit2 += convert[i][j] * momit[j];
		}

	}

}



inline void terrainmomskmeans()
{
	cout<<"reading in points..."<<flush;
	vector<mynum> pts;
	PointCloudIO::read(pts,"C:\\u\\max\\models\\Inga\\subset_all.txt","xyzi");
	cout<<pts.size()/3 <<"points read."<<endl;

	cout<<"creating tree ..."<<flush;
	int npts=pts.size()/3;
	TsTree2 t(pts);
	cout<<"ready"<<endl;

	mynum r=16;//for(mynum r=1;r<128;r*=2)
	{


		cout<<"computing moments, radius="<<r<<"..."<<endl;

		int nmoms = npts;
		vector<mynum> moms(nmoms*tmom::n);
		vector<mynum> momctr(nmoms*3);

		vector<mynum>::iterator it=moms.begin();
		vector<mynum>::iterator cit=momctr.begin();
		tensorSet3D<> ts,tsbuf;
		tmom bufm;
		for(it=moms.begin();it!=moms.end();it+=tmom::n,cit+=3)
		{
			int ptid;
			/*ptid= int(float(npts)*rand()/RAND_MAX);
			if(ptid>=npts)ptid=npts-1;
			ptid*=3;
			 */
			ptid=cit-momctr.begin();
			std::copy(pts.begin()+ptid,pts.begin()+ptid+3,cit);
			t.sumInCircle(ts,tsbuf,(TsTree2::mynum*)&*cit,r,-1);
			ts.transinv();
			kdtreeEntryFromTset(&*it,bufm,ts,r);

			if(ptid%3000==0)
			{
				cout<< 100.*ptid/momctr.size()<<"%\r"<<flush;
			}
		}
		cout<<"ready"<<endl;


		PointCloudIO::write(moms,"momsfordanielengel.csv","csv28");

		cout<<"converting moments to unit covariance..."<<flush;
		vector<mynum> newmoms;
		int dim;
		tounitcov(newmoms,moms,dim);
		cout<<"dim. is:"<<dim<< "..ready"<<endl;


		//mink is the desired number of clusters
		int mink=3;//for(int mink=4;mink<=64;mink*=2)
		{

			cout<<"clustering"<<endl;
			vector<int> cellpermom;
			int nbx = cluster::cluster_boxes(cellpermom,newmoms,dim,mink);
			cout<<"kmeans..."<<endl;

			vector<mynum> kmeans_centers(dim*nbx);

			it=newmoms.begin();
			vector<int> num(nbx,0);
			for(int i=0;i<nmoms;++i,it+=dim)
			{
				vector<mynum>::iterator it2 = kmeans_centers.begin() 
							+ cellpermom[i] *dim;
				num[cellpermom[i]]++;

				for(int j=0;j<dim;++j)
					it2[j]+=it[j];
			}

			it=kmeans_centers.begin();
			for(int i=0;i<nbx;++i,it+=dim)
			{		
				mynum m=1.0/num[i];
				for(int j=0;j<dim;++j)
					it[j] *=m;			
			}

			cluster::kmeans_simple<mynum,cluster::eukliddist<mynum> >(kmeans_centers,cellpermom,newmoms,dim,0,20,0.1);
			cout<<"ready"<<endl;

			ostringstream os;
			os<<"C:\\u\\max\\models\\Inga\\subset_all_r="<<r<<"k="<<nbx<<".osg"<<ends;
			ofstream out(os.str().c_str());
			cout<<"writing osg file..."<<flush;
			//	points_colorcoded_osg(momctr,cellpermom,out);
			points_colored_momdist_osg(momctr,newmoms,kmeans_centers,dim,out);
			cout<<"ready"<<endl;
		}
	}
}

//create maple basis functions
//depending on the radius
//for a unit sphere
template<class T>
void bla(ostream&out)
{

	vector<vector<T> > p(10);

	//the squares of p's norm;
	vector<T> pn(10);

	for(int i=0;i<10;++i)
	{
		p[i].resize(i+1,T(0));
		p[i][i]=1;

		//orthogonalize
		for(int j=0;j<i;++j)
		{
			T sprod = 0;
			for(int k=0;k<=i;++k)
				for(int l=0;l<=j;++l)
					sprod += p[i][k] * p[j][l] / T( k+l+2);

			sprod/=pn[j];
			for(int l=0;l<=j;++l)
				p[i][l]-=p[j][l]*sprod;

		}

		//calculate norm,squared
		pn[i] = 0;
		for(int k=0;k<=i;++k)
			for(int l=0;l<=i;++l)
				pn[i] += p[i][k] * p[i][l] / T( k+l+2);

		out<<"p["<<i<<"]:=(";

		for(int l=0;l<=i;++l)
		{
			/*			if(p[i][l].n<0)
			{p[i][l].n=-p[i][l].n;p[i][l].z=-p[i][l].z;}*/

			out << (l>0 ? " + ":"")<<'('<<p[i][l] << ")*r^"<<l;
		}
		/*if(pn[i].n<0)
		{pn[i].n=-pn[i].n;pn[i].z=-pn[i].z;}*/
		out<<")/sqrt(2*Pi*("<<pn[i]<<"));\n";


	}
}







//uncomment whatever you want to.
void momStatsmain(int argc, const char*argv[])
{

	//ofstream blam("c:\\u\\max\\bla.maple");
	//bla<infbruch>(blam);

	for(int i=0;i<argc;++i)
		cout<<i<<" "<<argv[i]<<endl;

	//terrainmomskmeans();

#ifdef GNUCC
	cout<<"Clock ticks per second:"<<ticks_per_sec<<endl;
#endif
	//convert2raw(argc,argv);

	//	testreadGAVABvrmlpts();

	/*
	{
	ofstream o("test.osg");
	osgcircball(o);
	o.close();  
	}
	 */
	//makecovrandsurf();
	/*
#ifdef GNUCC
	testTsTree("dino.dump");
#else
	testTsTree();
#endif
	//testspeedtranslate();

	 */
	/*



	//createstats(argc,argv);
	ofstream stp;

	const char*nam[]={0,"face.moments","noise.moments","klaus.moments","debug.moments","bunny.moments",0};

	stp.open("stats_face.txt");
	for(int i=1;nam[i]!=0;++i)
	{
		momStats2 xx(nam[i]);
		xx.write(stp);
	}
	stp.close();


	stp.open("stats_face.ps");
	writestatsps2("stats_face.txt",stp);
	writeCorrelationPS("stats_face.txt",stp);
	stp.close();

	 */
	//testsolve();


	//testObjectdb(cout);	


	//testObjectdbAndLog();	

	//testWriteMomHistPS();


	//createTable();
	//sft();


	/*
	createascs(1,0,"C:\\u\\max\\models\\3D_RMA\\all_auto\\ysghl2v1.xyz");
	createascs(1,0.07,"C:\\u\\max\\models\\3D_RMA\\all_auto\\ysghl2v1.xyz");
	createascs(0.4,0.07,"C:\\u\\max\\models\\3D_RMA\\all_auto\\ysghl2v1.xyz");
	 */

	//writesmoothfun(10);

}

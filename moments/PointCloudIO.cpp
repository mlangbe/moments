/*
	I/O of pointclouds to/from different formats

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
#include<iostream>
#include<fstream>
#include<cassert>
#include<stdlib.h>
#include<sys/stat.h>
#include<sys/types.h>
#include<fcntl.h>
#ifndef GNUCC
#include<io.h>
#else
#include<errno.h>
#include<unistd.h>
#endif

#include "PointCloudIO.h"
#include "equaldist.h"
#include "osgio.h"

#ifndef OLD_VC
#include<string.h>
#include<stdexcept>
#define _EXCEPTION_ std::runtime_error
#else
#define _EXCEPTION_ exception
#endif


using namespace std;
PointCloudIO::PointCloudIO(void)
{
}

PointCloudIO::~PointCloudIO(void)
{
}

bool PointCloudIO::accept(const char *name, const char *format)
{

	if(format)
		return strcmp(format,this->format())==0;

	const char*s=strstr(name,ext());
	if(!s)return false;
	if(strlen(s)!=strlen(ext()))
		return false;

	return true;

}


PointCloudIO::Factory* PointCloudIO::instance=0;

PointCloudIO* PointCloudIO::Factory::getIO(const char *name, const char *format)
{
	for(int i=0;i<instances.size();++i)
		if(instances[i]->accept(name,format))
			return instances[i];

	cerr<<"no appropriate reader or writer for "<<name<<" found";
	throw;
}

//--------------------------------------------------------------------------
//now instances of reader/writers:


struct pcio_raw : public PointCloudIO
{


	
	void read_(vector<mynum> & points ,const char*fname)
	{
		struct stat s;
		if(stat(fname,&s)<0)
		{throw _EXCEPTION_("couldn't read");}

		points.resize(s.st_size/sizeof(mynum));
		FILE*f=fopen(fname,"rb");
		fread(&points[0],sizeof(points[0]),points.size(),f);
		fclose(f);
	}



void write_(const vector<mynum> & points ,const char*fname)
{
/*	struct stat s;
	errno=0;
	int fd=open(fname,O_CREAT|O_RDWR,S_IWRITE|S_IREAD);
	if(fd<0)
	{	cerr<<"couldn't open "<<fname<<" for writing:errno="<<errno<<endl;throw exception("couldn't write");}
	else
	{
		::write(fd,&points[0],sizeof(points[0])*points.size());
		close(fd);
	}
	*/
	FILE *f=fopen(fname,"wb");
	if(f==0)
	{	cerr<<"couldn't open "<<fname<<" for writing:errno="<<errno<<endl;throw _EXCEPTION_("couldn't write");}
	else
	{
		//Achtung: Dateigroesse bei 64-bit-Version ist groesser als die geschriebenen Daten
		size_t num=fwrite(&points[0],sizeof(points[0]),points.size(),f);
		fclose(f);
	}

}

 const char *ext(){return ".raw";}


	pcio_raw()
	{
		reg();
	}
} pcio_1;


struct pcio_rw2 : public PointCloudIO
{

	//typedef unsigned __int64 tlen;
	typedef u_int64_t tlen;

	
	void read_(vector<mynum> & points ,const char*fname)
	{
		FILE*f=fopen(fname,"rb");

		if(f==0)
		{
			cerr<<"could not open "<<fname<<endl;
			throw _EXCEPTION_("could not read");
		}

		tlen len;
		tlen x;
		fread(&x,sizeof(x),1,f);
		fread(&len,sizeof(len),1,f);
		if(x!=sizeof(mynum)){
			printf("\n%s: number size=%d bytes != read size=%d bytes\n",fname,(int)x,(int)sizeof(mynum));
			throw _EXCEPTION_("number format does not match");
		}
		points.resize(len);
		size_t nread=fread(&points[0],x,len,f);
		if(nread!=len){
			printf("\n%s: length=%d != read=%d \n",fname,(int)len,(int)nread);
			throw _EXCEPTION_("could not read full length");
		}
		fclose(f);
	}



void write_(const vector<mynum> & points ,const char*fname)
{
	FILE *f=fopen(fname,"wb");
	if(f==0)
	{	cerr<<"couldn't open "<<fname<<" for writing:errno="<<errno<<endl;throw _EXCEPTION_("couldn't write");}

	tlen x=sizeof(mynum);
	fwrite(&x,sizeof(x),1,f);
	tlen len=points.size();
	fwrite(&len,sizeof(len),1,f);
	//Achtung: Dateigroesse bei 64-bit-Version ist groesser als die geschriebenen Daten
	size_t num=fwrite(&points[0],sizeof(points[0]),len,f);
	
	fclose(f);


}

 const char *ext(){return ".rw2";}


	pcio_rw2()
	{
		reg();
	}
} pcio_1_1;



struct pcio_asc:public PointCloudIO
{

	void write_(const vector<mynum> & points ,
									const char*fname)
{
	ofstream out(fname);
	vector<mynum>::const_iterator i;
	for(i=points.begin();i<points.end();)
		out <<*(i++)<<' '
				<<*(i++)<<' '
				<<*(i++)<<endl;
}

	const char*ext(){ return ".asc"; }

	pcio_asc(){reg();}

} pcio_2;



struct pcio_gavabvrml:public PointCloudIO
{
void read_(vector<mynum> & points ,
									const char*fname)
{
	char buf[80];
	ifstream in(fname);
	
	for(;;){
		in.getline(buf,79);
		if(strstr(buf,"point [")!=0)
			break;
		if(in.eof())return;
	}

	mynum x,y,z;
	char cc;
	in.clear();
	points.clear();
	for(;;){
		in>>x>>y>>z;
		if(in.fail())break;
		in>>cc;
		points.push_back(x);
		points.push_back(y);
		points.push_back(z);
	}

	in.close();
}
	const char*ext(){ return ".wrl"; }

	pcio_gavabvrml(){reg();}

} pcio_3;

//create the points eually distributed over the surface.
struct pcio_gavabvrml_equaldist:public PointCloudIO
{
	void read_(vector<mynum> & points ,
		const char*fname)
	{

		equaldist eq;

		char buf[80];
		ifstream in(fname);

		for(;;){
			in.getline(buf,79);
			if(strstr(buf,"point [")!=0)
				break;
			if(in.eof())return;
		}

		mynum x,y,z;
		char cc;
		in.clear();
		eq.points.clear();
		for(;;){
			in>>x>>y>>z;
			if(in.fail())break;
			in>>cc;
			eq.points.push_back(x);
			eq.points.push_back(y);
			eq.points.push_back(z);
		}



		vector<int> bufidx;
		in.clear();

		for(;;){
			in.getline(buf,79);
			if(strstr(buf,"coordIndex [")!=0)
				break;
			if(in.eof())return;
		}

		int vidx;
		for(;;){
			in>>vidx;
			if(in.fail())break;
			if(vidx>=0)
				bufidx.push_back(vidx);
			else
			{
				if(bufidx.size()==3)
				{
					eq.tris.resize(eq.tris.size()+3);
					copy(bufidx.begin(),bufidx.end(),eq.tris.end()-3);
					bufidx.clear();
				}
				else if(bufidx.size()==4)
				{
					eq.quads.resize(eq.quads.size()+4);
					copy(bufidx.begin(),bufidx.end(),eq.quads.end()-4);
					bufidx.clear();
				}
				else
				{
					eq.polygons.push_back(bufidx);
					bufidx.clear();
				}

			}

			in>>cc;
		}





		in.close();


		eq.createEqualDist(points,eq.points.size()/3);


	}
	const char*ext(){ return ".wrl"; }
	const char*format(){ return "equaldist"; }

	pcio_gavabvrml_equaldist(){reg();}

} pcio_3_1;


struct pcio_3drma:public PointCloudIO
{

void read_(vector<mynum> & points ,
									const char*fname)
{
	
	ifstream in(fname,ios::binary);

	points.clear();
	short int n;
	assert(sizeof(short int)==2);
	assert(strstr(fname,".xyz")!=0);
	assert(sizeof(mynum)>=2);

	while(!in.eof())
	{
		n=0;
		in.read((char*)&n,2);
		
		if(n==0)continue;

		int onpt=points.size();		
		points.resize(points.size()+n*3);

		//use output vector as read buffer
		short int*buf=(short int*)&points[onpt];
		in.read((char*)buf,n*3*sizeof(short int));

		//assign backwards, so you don't overwrite values that are still to be converted
		for(int i=n*3-1;i>=0;--i)
		{
				points[onpt+i]=buf[i];
		}
	}
}


void write_(const vector<mynum> & points ,
						const char*fname)
{
	ofstream out(fname,ios::binary);
	short int n;

	assert(sizeof(short int)==2);
	static const int nmax=((1<<15)-1);
	short int buf[1+3*nmax];
	buf[0]=(short int)nmax;

	vector<mynum>::const_iterator it;
	
	for(it=points.begin();
	    points.end()-it >=nmax*3;
			)
	{
		short int*pts=buf+1, *ptsend=pts+nmax*3;
		for(;pts!=ptsend;++pts,++it)
			*pts=(short int)*it;
	
		out.write((const char*)buf,sizeof(buf));		
	}

	buf[0] = (short int) ((points.end()-it)/3); 

	for(short int*pts=buf+1;it!=points.end();++pts,++it)
		*pts=(short int)*it;
	
	out.write((const char*)&n,sizeof(n));
	out.write((const char*)buf,(int(buf[0])*3+1)*sizeof(buf[0]));

}


	
	const char*ext(){ return ".xyz"; }
	pcio_3drma(){reg();}


} pcio_4;


struct pcio_scandump_denker:public PointCloudIO
{
void read_(vector<mynum> & points ,
									const char*fname)
{


	ifstream in(fname);

	if(!in.is_open())
	{cout<<"couldn't open"<<fname<<endl;
	return;}

	char buf[80];

	points.clear();

	while(!in.eof())
	{
		buf[0]=0;
		in.getline(buf,80,'\n');
		if('A'<=buf[0] && buf[0] <='Z' || buf[0]==0)
			continue;

		size_t ps=points.size();
		points.resize(points.size()+3);
		
#ifdef GNUCC
  		const char*formatstr="%Lf,%Lf,%Lf";
#else
		const char*formatstr="%lf,%lf,%lf";
#endif
		
		if( sscanf(buf,formatstr,
			&points[ps],&points[ps+1],&points[ps+2])!=3 )
		{cout<<"sth wrong at pos:"<<ps<<":"<<buf<<endl;}
	}
	in.close();
}

const char*ext(){return ".dump";}
  pcio_scandump_denker()
	{reg();}


} pcio_5;


/*
point cloud from zaman with 4 values
per row:
x,y,z,intensity
*/
struct pcio_zaman:public PointCloudIO
{

	void read_(vector<mynum> & points ,
									const char*fname)
{
	ifstream in(fname);

	if(!in.is_open())
	{cout<<"couldn't open"<<fname<<endl;
	return;}

	char buf[80];

	points.clear();

	mynum x,y,z,v;
	while(!in.eof())
	{
		in>>x>>y>>z>>v;
		if(in.fail())
			break;

		points.push_back(x);
		points.push_back(y);
		points.push_back(z);

	}
	in.close();
}

	const char*ext(){ return ".txt"; }
	const char*format()  {return "xyzi";}
	pcio_zaman(){reg();}

} pcio_6;

/*
point cloud from zaman with 3 values
per row:
x,y,z
*/
struct pcio_csv:public PointCloudIO
{

	void read_(vector<mynum> & points ,
									const char*fname)
{
	ifstream in(fname);

	if(!in.is_open())
	{cout<<"couldn't open"<<fname<<endl;
	return;}

	char buf[80];

	points.clear();

	mynum x,y,z;
	while(!in.eof())
	{
		in>>x;
		in.ignore(1,',');
		in>>y;
		in.ignore(1,',');
		in>>z;
		if(in.fail())
			break;

		points.push_back(x);
		points.push_back(y);
		points.push_back(z);

	}
	in.close();
}

	void write_(const vector<mynum> & points ,
						const char*fname)
{
	ofstream out(fname);
	
	vector<mynum>::const_iterator it;	
	for(it=points.begin();
	    it!=points.end();
			it+=3
			)
	{
		out<<it[0]<<','<<it[1]<<','<<it[2]<<'\n';
	}
}

	const char*ext(){ return ".csv"; }
	pcio_csv(){reg();}

} pcio_7;

/*
point cloud from zaman with 4 values
per row:
x,y,z,intensity
*/
struct pcio_csv28:public PointCloudIO
{

	void read_(vector<mynum> & points ,
									const char*fname)
{
	ifstream in(fname);

	if(!in.is_open())
	{cout<<"couldn't open"<<fname<<endl;
	return;}

	char buf[80];

	points.clear();

	mynum x;
	while(!in.eof())
	{
		in>>x;
		if(in.fail())
			break;
		in.ignore(1,',');

		points.push_back(x);
	}
	in.close();
}


	void write_(const vector<mynum> & points ,
						const char*fname,int n)
{
	ofstream out(fname);
	
	assert(points.size()%n==0);

	out.precision(sizeof(mynum)*2);
	vector<mynum>::const_iterator it;	
	for(it=points.begin();
	    it!=points.end();
			it+=n
			)
	{
		out<<it[0];
		for(int i=1;i<n;++i)
			out<<','<<it[i];
		out<<"\r\n";
	}
}

	void write_(const vector<mynum> & points ,
						const char*fname)
	{
		write_(points,fname,28);
	}


	const char*ext(){ return ".csv"; }
	const char*format(){ return "csv28"; }
	pcio_csv28(){reg();}

}pcio_8;


//create the points eually distributed over the surface.
struct pcio_off_equaldist:public PointCloudIO
{
	void read_(vector<mynum> & points ,
		const char*fname)
	{

		equaldist eq;

		char buf[80];
		ifstream in(fname);

		in>>buf;

		int nverts,nfaces,nedges;

		in>>nverts>>nfaces>>nedges;

		mynum x,y,z;
		char cc;
		in.clear();
		eq.points.clear();
		for(int i=0;i<nverts;++i){
			in>>x>>y>>z;
			if(in.fail())break;
			eq.points.push_back(x);
			eq.points.push_back(y);
			eq.points.push_back(z);
		}



		in.clear();

		int vertsperface;
		for(int i=0;i<nfaces;++i){
			in>>vertsperface;

			if(vertsperface==3)
			{
				int a,b,c;
				in>>a>>b>>c;
				eq.tris.push_back(a);
				eq.tris.push_back(b);
				eq.tris.push_back(c);
			}
			else if(vertsperface==4)
			{
				int a,b,c,d;
				in>>a>>b>>c>>d;
				eq.tris.push_back(a);
				eq.tris.push_back(b);
				eq.tris.push_back(c);
				eq.tris.push_back(d);
			}
			else
			{

					eq.polygons.resize(eq.polygons.size()+1);
					vector<int> &last=eq.polygons.back();
					last.resize(vertsperface);
					for(int i=0;i<vertsperface;++i)
					{
						in>>last[i];
					}
			}

		}





		in.close();

		int npoints=(eq.points.size()/3) *10;
		if(npoints <1000)npoints= 1000;
		if(npoints>50000)npoints=50000;
		eq.createEqualDist(points,npoints);

	}
	const char*ext(){ return ".off"; }
	const char*format(){ return "offequaldist"; }

	pcio_off_equaldist(){reg();}

} pcio_9_1;

struct pcio_osg:public PointCloudIO
{
	void read_(vector<mynum> & points ,
		const char*fname)
	{
		std::ifstream in(fname);
		osgreadpoints(points,in);
		in.close();
	}

	void write_(const vector<mynum> & points ,
		const char*fname)
	{
		std::ofstream of(fname);
		osgpointcloud(of ,points);
		of.close();
	}

	const char*ext(){ return ".osg"; }

	pcio_osg(){reg();}

} pcio_10;


#undef _EXCEPTION_

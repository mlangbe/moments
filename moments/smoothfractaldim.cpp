/*
	"smooth fractal dimension"- calculation for point clouds

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

#include "PointCloudIO.h"
#include "TsTree2_impl.h"
#include "vecstats.h"
#include "stdio.h"
#include<vector>
#include<algorithm>
#include<cmath>
#include<fstream>
#include <time.h>
using namespace std;

typedef unsigned long long tzindex;




void outputPS(std::vector<tzindex> &zindex,ostream&out,int num,bool nonum)
{
	if(num>zindex.size())num=zindex.size();

	//calculate "interesting" area
	int minchange=sizeof(tzindex)*8;
	int maxchange=0;
	for(int i=1;i<num;++i)
	{
		tzindex diff = zindex[i-1] ^ zindex[i]; 
		int nch=0;
		while(diff!=0) { diff>>=3;nch++; }
		if(nch<minchange)minchange=nch;
		if(nch>maxchange)maxchange=nch;	
	}


	if(minchange>0)minchange-=1;
	out<<"%!PS-Adobe-3.0 EPSF-3.0"<<endl;
	out<<"%%BoundingBox: 0 0 "<<num*10 <<" "<<(maxchange-minchange)*30<<endl;
	out<<"/Times-Roman findfont 10 scalefont setfont"<<endl;
	
	for(int i=0;i<num;++i)
	{
		if(i>0){
			tzindex diff = (zindex[i-1] ^ zindex[i])>>minchange*3; 
			int nch=0;
			while(diff!=0) { diff>>=3;nch++; }
			out<<i*10<<" 0 moveto "<<i*10<<' '<<nch*3*10 <<" lineto stroke"<<endl;
		}	

		if(!nonum){
		tzindex z = zindex[i];
		z>>=(minchange*3);
		for(int j=0;j<(maxchange-minchange)*3;++j)
		{
			int zz=((z>>j)&1);

			if(zz)
				out<<" 0.8 setgray ";
			else
				out<<" 0.9 setgray ";			
			out<<i*10+1<<' '<<j*10<<" moveto "
				"8 0 rlineto 0 10 rlineto -8 0 rlineto fill ";

			out<<"0 setgray ";
			out<<i*10+2<<' '<<j*10+2<<" moveto ("<< zz <<") show ";
		}
		}
		out<<endl;


	}
	out<<"0.7 setgray"<<endl;
	for(int j=0;j<(maxchange-minchange);++j)
	{
		out<<"0 "<<j*30<<" moveto "<<num*10<<' '<<j*30<<" lineto ";
	}
	out<<" stroke "<<endl;

}


void smoothfractaldim()
{
	
	for(int i=0;i<1<<30;++i)
	{
		int x,y,z;
		invzord(x,y,z,i);
		if(zord2(x,y,z)!=i){
			cout<<"error at i="<<i<<endl;
			return;
		}
		if((i&0xffff)==0)
			cout<<'\r'<<i * 100. / (1<<30)<<"%"<<flush;
	}
	

	for(int i=0;i<64;++i)
	{
		int x,y,z;
		invzord(x,y,z,i^(i>>1));
		//cout<<oct<<z<<' ';

		cout<<oct<<"i"<<i<<"->("<<x<<','<<y<<','<<z<<')'<<endl;	
	}



	static const int octaves=sizeof(tzindex)*8/3;
	static const int numsub=12;
	

	cout<<"Calculation of smooth fractal dimension:"<<endl;

	cout<<"reading pointcloud:"<<endl;

	std:vector<mynum> points;

	const char*name="C:\\u\\max\\models\\denker\\da\\dino.raw";
	PointCloudIO::read(points,name);
	cout<<"done."<<endl;

	char oname[256],oname2[256],oname3[256];
	sprintf(oname,"%.*s_fractaldim.dat",strlen(name)-4,name);
	sprintf(oname2,"%.*s_zindex.eps",strlen(name)-4,name);
	sprintf(oname3,"%.*s_zindex2.eps",strlen(name)-4,name);
	
	ofstream out2(oname2);
	ofstream out3(oname3);
	cout<<"output will be written to:"<<oname<<endl;

	cout<<"collecting stats:"<<endl;
	vecstatscollector<3,mynum> vc;
	std::vector<mynum>::iterator it=points.begin(),itend=points.end();
	for(;it<itend;it+=3)
		vc.add(&*it);

	vecstats<3,mynum> vs = vc;
	cout<<"done."<<endl;

	double maxlen=0;
	for(int i=0;i<3;++i)
	{
		double d=vs.max[i]-vs.min[i];
		if(maxlen<d)maxlen=d;
	}

	
	double mulsub = pow(0.5 , 1.0/numsub);



	//the output array
	size_t num[octaves][numsub];
	size_t num2[octaves][numsub];

	std::vector<tzindex> zindex;

	zindex.resize(points.size()/3);

	double mul = ldexp(1.0/maxlen,octaves);
	for(int s=0;s<numsub;++s,mul*=mulsub)
	{
		vector<tzindex>::iterator zit = zindex.begin(),zitend = zindex.end();
		cout<<"calculating zindices box length="<<1.0/mul<<endl;

		clock_t tstart=clock();


		for(it=points.begin();it<itend;it+=3,++zit)
		{
			*zit = zord2<tzindex>	(
				(tzindex)(  ( it[0] - vs.min[0] )*mul ),
				(tzindex)(  ( it[1] - vs.min[1] )*mul ),
				(tzindex)(  ( it[2] - vs.min[2] )*mul )								
				)		;
		}
		
		float tim=(clock()-tstart)*1.0f/CLOCKS_PER_SEC;
		cout<<" used "<< tim*1000 <<" ms for calculating all "<< sizeof(tzindex)*8<<" bit  zindices = "<<tim/zindex.size()*1e9<<"ns per index"<<endl;
		tstart=clock();

		sort(zindex.begin(),zindex.end());

		cout<<" used "<<(clock()-tstart)*1e3/CLOCKS_PER_SEC<<" ms for sorting zindices "<<endl;
	

		if(s==0){
			outputPS(zindex,out2,50,false);
			outputPS(zindex,out3,50,true);
		}

		cout<<"extracting numbers:"<<endl;

		for(int o=0;o<octaves;++o)
			num[o][s]=1;		

		zit=zindex.begin();
		tzindex ozind = *zit;
		for( ++zit;zit!=zitend;++zit)
		{		
			tzindex diff = *zit ^ ozind; 
			ozind = *zit;		


			for(int o=0; o < octaves && (diff>>(3*o))!=0 ;++o)
				++num[o][s];				
		}				
	}
	


	
	
	ofstream out(oname);
	//out<<"'box size' 'box number'"<<endl;

	mul=ldexp(1.0/maxlen,octaves);
	for(int o=0;o<octaves;++o)
	{
		out<<"# octave number "<<o<<endl;
		for(int s=0;s<numsub;++s,mul*=mulsub)
		{
			out<< 1.0/mul<< ' '<< num[o][s]<<'\n';
		}
	}
	out.close();

	cout<<"ready with fractal dim for " <<points.size()/3<<" point pointcloud "<<endl;
	cout<<"sizeof(void*)="<<sizeof(void*)<<endl;
	
	char bufc;
	cin>>bufc;

}

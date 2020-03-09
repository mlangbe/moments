/*
	This file is part of example code using interval arithmetics for implicit surface triangulation

	Copyright (C) 2020 TU Kaiserslautern, Prof.Dr. Hans Hagen (AG Computergraphik und HCI)
	
	Author: Max Langbein	

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
#include <stdio.h>
#include<vector>
using namespace std;
#include "myinterval.h"
#include "intervalnewton.h"

#include "Implicit.h"
#include "d3op.h"
#include<assert.h>

#ifndef OLD_VC
#include<string.h>
#endif

class extractsurface
{
public:

	inline extractsurface()
	{
		memset(mctable,0,sizeof(mctable));

		initmctable();
	}
	int caseCounts[256];
	int caseCounts2[256];


	inline void toStl(const Implicit &x, const char*name, float epsilon=0.1,myinterval<float>* box=0,bool binary=true)
	{
		binaer=binary;
		eps=epsilon;
	
		printf("\n%s\n",name);
		memset(caseCounts,0,sizeof(caseCounts));
		memset(caseCounts2,0,sizeof(caseCounts));
		memset(wrongs,0,sizeof(wrongs));
		memset(rights,0,sizeof(rights));
	
		out=fopen(name,binaer?"wb":"wt");

		myinterval<float> bigbox[]		
		={
			myinterval<float>(-1<<10,1<<10),
			myinterval<float>(-1<<10,1<<10),
			myinterval<float>(-1<<10,1<<10)
		};

		if(box==0)
			box=bigbox;


		ntris=0;ontris=0;

		if(binary)
		{
			const char* header=
				"binary stl : hallo  "
				"binary stl : hallo  "
				"binary stl : hallo  "
				"binary stl : hallo ";
			fwrite(header,1,80,out);
			fwrite(&ntris,4,1,out);
		}
		else
		{
			fprintf(out,"solid hallo\n");
		}

		extract(x,box);

		if(binary)
		{
			fseek(out,80,SEEK_SET);
			//write number of triangles
			fwrite(&ntris,4,1,out);
		}
		else
			fprintf(out,"endsolid hallo\n");

		fclose(out);

		printf("\rntris: %d\n",ntris);


		/*
		printf("wrong orientated:");
		for(int i=0;i<256;++i)
		{
			if(caseCounts[i]>0)
			{
				printf("%d : %d/%d\n",i,caseCounts2[i], caseCounts[i]);
			}
		
		}
		*/


	}





private:

	float eps;
	FILE*out;

	bool binaer;
	int ntris;


	myinterval<float> extract(const Implicit &x, myinterval<float>* box)
	{

		myinterval<float> val=x(box);

		if(!val.contains(0))
			return val;;

		myinterval<float> * maxi=box+0;
		float maxw=box[0].width();
		float bw=box[1].width();
		if(maxw<bw){maxw=bw;maxi=box+1;}
		bw=box[2].width();
		if(maxw<bw){maxw=bw;maxi=box+2;}

		if(maxw<eps){

			//if case found:return
			if(mctris(x,box))
				return val;

			return val;
			

			if(bw<eps*0.1)
				return val;
		}
		myinterval<float> oi=*maxi;
		float mid=oi.mid();
		maxi->b=mid;
		myinterval<float> val1=extract(x,box);
		maxi->b=oi.b;
		maxi->a=mid;
		myinterval<float> val2=extract(x,box);
		maxi->a=oi.a;


		return val;
	}




	struct mirrorOp
	{
		int pos;
		int value;
		int newc[3];

		mirrorOp(int opcode)
		{
			init(opcode);
		}


		/**

		decode opcode.

		mirrorOp  = 	((((((((newx)<<1) | (newy-newx-1) )<<1) | mirrorvalues)<<1)| mirrorz)<<1)| mirrory)<<1)| mirrorx;  

		coordxch equiv_3  newx*3 + (1+newy-newx)*2   
		*/

		inline void init(int mirror)
		{
			pos=(mirror&7);
			mirror>>=3;
			value= (-(mirror&1))&0xff ;
			mirror>>=1;
			newc[0] = mirror>>1;
			newc[1] = (newc[0]+ 1+ (mirror&1) ) % 3 ;
			newc[2] = (newc[0]+ 1+ ((mirror&1)^1) ) % 3 ;
		}


		inline int mirrorVertex(int i)
		{
			int j =i^ pos;
			return ((j&1)<<newc[0]) | (((j&2)>>1)<<newc[1]) | (( (j&4)>>2)<<newc[2]);
		}


		inline int mirrorCase(int oi)
		{
			int ret=0;
			for(int i=0;i<8;++i)
			{
				int j=mirrorVertex(i);
				ret|= ((oi>>i)&1)<<j;
			}
			return ret^value;
		}



		inline int mirrorEdge(int oedge)
		{
			int nvtx = mirrorVertex(oedge&7);
			int ndir=newc[oedge>>3];
			nvtx &= ~ (1<<ndir);
			return (ndir<<3)|nvtx;
		}

		/** have tris orientation to be mirrored?*/
		inline int mirrortris()
		{
			return ((value^pos^(pos>>1)^(pos>>2))&1) ^  ( (newc[1]-newc[0]+3)%3 !=1)  ;
		}

	
	
	
	};




	struct mcentry
	{
		//triangles using edge ids as edge[ j*8 + i ] = edge i in direction j, i being the lower vertex index of the edge 
		int ntris;
		int tris[3*4];
		int nedges;
		int edges[12];


		int code;

		int origCode;

		int opcode;

		bool checkEdges(int code)
		{

			bool et[24];
			memset(et,0,sizeof(et));

			//calculate required vertices on edges 
			for(int i=0;i<ntris*3;++i)
			{
				int ei=tris[i];
				et[ei]=true;
			}


			for(int ei=0;ei<24;++ei)
			{			
				int v1=ei&7;
				int dir=ei>>3;
				int v2=v1| (1<< dir);
			
				if(v1==v2)
					continue;


				if( et[ei] != ( ((code>>v1)&1) != ((code>>v2)&1) ) )
				{
					printf("(%d-%d:%sthere)",v1,v2,et[ei]?"":"not "); 
					return false;
				}
			}

			return true;
		}

		void initEdges()
		{
			bool et[24];
			
			memset(et,0,sizeof(et));

			for(int i=0;i<ntris*3;++i)
			{
				et[tris[i]]=true;	
			}

			nedges=0;
			for(int i=0;i<24;++i)
			{
				if(et[i])
				{
					edges[nedges++]=i;
				}
			}
		}


		void mirroredOf(mcentry &o,mirrorOp&mo)
		{
			ntris=o.ntris;
			int m=mo.mirrortris();
	
			for(int i=0;i<ntris;++i)
			{
				 int t=i*3;
				 tris[t] = mo.mirrorEdge(o.tris[t]);
				 tris[t+1+ m] = mo.mirrorEdge(o.tris[t+1]);
				 tris[t+1+!m] = mo.mirrorEdge(o.tris[t+2]);
			}
						
			nedges=o.nedges;
			for(int i=0;i<nedges;++i)
			{
				edges[i]=mo.mirrorEdge(o.edges[i]);
			}


		}
	};


	mcentry mctable[256];



	int normalizeCase(int  c)
	{
		int ret=255;
		for(int mirrorCode=0;mirrorCode< 3*32;++mirrorCode)
		{

			mirrorOp mo(mirrorCode);
			int code=mo.mirrorCase(c);
			if(ret>code)
				ret=code;
		}

		return ret;

	}


			
	int addMirroredCases(mcentry&x,int c,bool noValMirror=false)
	{
		if(!x.checkEdges(c))
		{
			printf(" aua%d ",c);
		}

		x.initEdges();

		int ret=0;
		for(int mirrorCode=0;mirrorCode< 3*32;++mirrorCode)
		{

			mirrorOp mo(mirrorCode);

	
			if(noValMirror&&mo.value || mo.mirrortris())
			{
				continue;
			}



			int code=mo.mirrorCase(c);
			mcentry&ent=mctable[code];
			if(ent.ntris==0)
			{
				ent.mirroredOf(x,mo);
				if(!ent.checkEdges(code))
				{
					printf(" \" ");
				}

				ent.opcode=mirrorCode;
				ent.origCode=c;
				ent.code=code;
				++ret;
			}
		}
		return ret;
	}


	void initmctable(){
	
		//edge directions
		int X=0,Y=8,Z=16;

		//vertex ids
		int x=1,y=2,z=4;


		int cases=0;



		//#1
		//case  1:  vertex 0 < 0           
		mcentry c1={  1, { X , Y,  Z }};
		

		cases+=addMirroredCases(c1,1);

		//#2
		//case 3  ( vertex 0,1 <0 ) 
		mcentry c3={ 2,  
		 { Y, Z, Y|x,
			Y|x, Z, Z|x,     
			 }  };
		cases+=addMirroredCases(c3,3);

		//#3
		//vertex 0, x|z <0
		//=case 1 +  (1<<5) =33
		mcentry c33={
			4,{X,Y,Z, X|z,Y|x|z,Z|x,  
			}
		};
		//no "value" mirror for this cases:
		//cases+=addMirroredCases(c33,33,true);

		//#3
		//vertex 0, x|z <0
		//=case 1 +  (1<<5) =33
		mcentry c33_1={
			4,{X,Y,Z, X|z,Y|x|z,Z|x,  
			//added to close surface if two of these cases are at one face:
			X,Z,X|z,
			X,X|z,Z|x
			}
		};
		//the additional tris only for the value-mirrored case:
		cases+=addMirroredCases(c33_1,33);



		//#4
		//( vertex x,y|x,y <0 ) -> 1<<1 | 1<<3 | 1<<2 = 2+8+4 = 14
		mcentry c_4={ 3,  { X,Z|x,Y,   
						   Y,Z|x,Z|y,
						   Z|x,Z|x|y,Z|y
		}  
		};
		cases+=addMirroredCases(c_4,14);

		//#5
		//case 15: ( plane orthogonal to z: vertices 0,1,2,3  < 0 )
		mcentry c15={ 2,  { Z,  Z|x,Z|x|y,    
			Z, Z|x|y,Z|y  } };
		cases+=addMirroredCases(c15,15);

		//#6 : vertices  x,x|y,y,z 
		// = 1<<1 | 1<<2 | 1<<3 | 1<<4  =2|4|8|16 =30
		mcentry c30={
			4,{  Y,X,Z|x,
				 Y,Z|x,Z|y,
				 Z|x,Z|x|y,Z|y,

				 Z,Y|z,X|z
				 },
		};
		cases+=addMirroredCases(c30,30);


		//#7: vertices 0, x|y, x|z, y|z, <0
		//=  1 | 1<<3 | 1<<5 | 1<<6 = 1|8|32|64  = 105
		mcentry c105={
			4,{ X,Y,Z,
				Y|x,Z|x|y,X|y,
				Z|x,X|z,Y|x|z,
				Z|y,Y|z,X|z|y
			},
		};
		cases+=addMirroredCases(c105,105);




		//#8
		//case   ( vertex 0, y,x|y,z|y  <0 )
		//-> 1 + 4 + 8 +  64  =  77
		mcentry c_8={ 4,  {Z,X|y|z, Y|z,
							Z,X,X|y|z,
							X,Z|y|x,X|y|z,
							X,Y|x,Z|y|x },
		};
		cases+=addMirroredCases(c_8,77);
		
		//#9
		//vertex x=1, x|y =3, y =2, y|z =6
		//= 2+ 4+8+ 64  = 78
		mcentry c_9={ 4,  { Y,X,   Y|z,
			                X,Z|x|y,Y|z,
							X,Z|x,Z|x|y,
							Y|z,Z|x|y,X|y|z
						   },
		};
		cases+=addMirroredCases(c_9,78);


		//#10
		//vertex 0, x|y|z <0
		//case 1+ 1<<7 = 129 
		mcentry c_10={
			2,{X,Y,Z,  
			   X|y|z,Z|y|x,Y|z|x },
		};
		cases+=addMirroredCases(c_10,129);

		//#11
		//vertex 0, x, x|y|z <0
		//case 1+ 2+ 1<<7 = 131 
		mcentry c_11={
			3,{Y,Z,Y|x,
				Y|x,Z,Z|x,
			   X|y|z,Z|y|x,Y|z|x, }
		};
		cases+=addMirroredCases(c_11,131);


		//#12
		//vertex x, z, x|y|z <0
		//case 1<<1+ 1<<4+ 1<<7 = 146 
		mcentry c_12={
			3,{X,Z|x,Y|x,
				Z,Y|z,X|z,
			   X|y|z,Z|y|x,Y|z|x },
		};
		cases+=addMirroredCases(c_12,146);

		//#13
		//vertex 0, 0|z,  x|y, x|y|z <0
		//case 1<<0 +  1<<4  +  1<<3  +  1<<7  =  153 

		mcentry c_13={
			4,{X,  Y,    Y|z,
			   X,  Y|z,  X|z,
			   Y|x,  X|y|z,  X|y, 
			   Y|x, Y|x|z,X|y|z,     },
		};
		cases+=addMirroredCases(c_13,153);



	//printf("%d of %d cases covered",cases,256-2);

	
	}






	float edgebuf[8*3][3];


	template<class T>
	inline static void getVertex(T*buf,const myinterval<float>*box,int i)
	{
			buf[0]= i&1 ? box[0].b:box[0].a;
			buf[1]= i&2 ? box[1].b:box[1].a;
			buf[2]= i&4 ? box[2].b:box[2].a;
	}


struct singledir
{

	int coordid;
	mutable myinterval<float> val[3];
	const Implicit&imp;

	inline singledir(const Implicit&imp,int id,const float*vals)
		:imp(imp),coordid(id)
	{
		val[0]=vals[0];val[1]=vals[1];val[2]=vals[2];
	}
	
	inline myinterval<float> operator()(const myinterval<float>& x) const
	{
		val[coordid]=x;
		return imp(val);
	}
	inline float operator()(float x) const
	{
		val[coordid]=x;
		return imp(val).mid();
	}
};

struct singledirDeriv
{

	const int coordid;

	mutable myinterval<float> val[3],g[3];

	const Implicit&imp;

	inline singledirDeriv(const Implicit&imp,int id,const float*vals)
		:imp(imp),coordid(id)
	{
		val[0]=vals[0];val[1]=vals[1];val[2]=vals[2];
	}

	inline myinterval<float> operator()(const myinterval<float>& x) const
	{
		val[coordid]=x;
		imp.grad(g,val);
		return g[coordid];
	}
};


	int wrongs[256],rights[256];

	int ontris;

	/** output marching-cubes-triangles into stl-file  for this cube*/
	int mctris(const Implicit&x,myinterval<float>*box)
	{
		float vertVals[8];
		int code=0;

		myinterval<float> buf[3];
		for(int i=0;i<8;++i)
		{
			getVertex(buf,box,i);
			vertVals[i]=x(buf).mid();

			code|=  (vertVals[i]<0)<<i;

			if(!(vertVals[i]<0) &&!(vertVals[i]>=0))
			{
				printf("au");
			}
		}



		const mcentry&mce=mctable[code];



		caseCounts[ mce.origCode] ++;

		//calculate required vertices on edges 
		for(int i=0;i<mce.nedges;++i)
		{
			int ei=mce.edges[i];

			int v1=ei&7;
			int dir=ei>>3;
			int v2=v1| (1<< dir);

			float*fbuf=edgebuf[ei];

			getVertex(fbuf,box,v1);


			vector< myinterval<float> > ret;


			if(false){
				singledir a(x,dir,fbuf);
				singledirDeriv b(x,dir,fbuf);

				solve4<float,singledir,singledirDeriv>(ret,box[dir],eps*1e-2,a,b);
			}
			if(ret.size()==1)
			{
				fbuf[dir] =ret[0].mid();
			}
			else
			{
				float delt=(vertVals[v1]-vertVals[v2]);

				if(delt==0)
					printf("del0!");
				/*
				if(vertVals[v1]==0)
					printf("v10");
				if(vertVals[v2]==0)
					printf("v10");
					*/
				float interp = delt!=0? vertVals[v1] / delt :0.5;

				
				if(!(interp>=0 &&interp<=1))
					printf("ou");
					
			
				fbuf[dir] += box[dir].width()*interp;

			}
		}

		

		float normal[3];
		float bufa[3],bufb[3];
		//write triangles

		myinterval<float> g[3];
		bool hasz=false;
		for(int i=0;i<mce.ntris;++i)
		{
			int t=i*3;
			
			float* a = edgebuf[mce.tris[t]];
			float* b = edgebuf[mce.tris[t+1]];
			float* c = edgebuf[mce.tris[t+2]];

			d3op3(,bufa,=,b,-,a,);
			d3op3(,bufb,=,c,-,a,);
			d3kreuz(normal,bufa,bufb);

			if(d3prod(normal,normal)<1e-12)
			{
				hasz=true;
			}
			/*
			myinterval<float> nlen=sqrt(d3prod(normal,normal));
			
			x.grad(g,a);
			myinterval<float> dgrad=d3prod(g,normal);
			myinterval<bool> ltz= dgrad<nlen*-0.1;


		
			if(ltz.a)
				wrongs[mce.origCode]++;
			else if(!ltz.b)
				rights[mce.origCode]++;
			else if(nlen.b>0)
			{
				printf("[%.2g %.2g] %g,",dgrad.a,dgrad.b,nlen.b);
			}
			*/
			writeTri(normal,a,b,c);

		}

		if(hasz)
			caseCounts2[mce.origCode]++;

		if(!mce.ntris && code !=0&&code!=255)
			printf(" %d ",code);
	
		ntris+=mce.ntris;

		if(ntris-ontris>1000)
		{
			ontris=ntris;
			printf("\rntris:%d",ntris);
			fflush(stdout);
		}


		return mce.ntris;
	
	}


	void writeTri(float*normal,float*a,float*b,float*c)
	{
		if(binaer)
		{
			assert(sizeof(float)==4);
			fwrite(normal,4,3,out);
			fwrite(a,4,3,out);
			fwrite(b,4,3,out);
			fwrite(c,4,3,out);
			short bla=('|'<<8)|'|';
			fwrite(&bla,2,1,out);	
		}
		else
		{
			fprintf(out,
				"facet normal %f %f %f \n"
				  " outer loop\n"
				   "   vertex %f %f %f\n"
				   "   vertex %f %f %f\n"
				   "   vertex %f %f %f\n"
				   " endloop\n"
				   "endfacet\n",
				   normal[0],normal[1],normal[2],
					a[0],a[1],a[2],
					b[0],b[1],b[2],
					c[0],c[1],c[2]
			);
		}
	}

};


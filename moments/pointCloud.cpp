/*
	operations to create  moment tensors

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
#include<iostream>
#include<cassert>
#include<fstream>
#include<iostream>
#include<vector>
#include<cmath>
#include<sstream>

#include "pointCloud.h"
#include "momentSet.h"
#include "tensor.h"
#include "polynom.h"
#include "optimizedPoly.h"

using namespace std;




typedef faceMomSet::valuetype mynum;







/*
* add 1, x, x x, x x x to the appropriate tensor
* in A with contains the elements in the order
* here example N=2:
* 0A,  1A3, 1A33,    1A2, 2A23,    1A22,
* 1A1, 2A12  ,  2A13,
* 2A11
* so the fastest changing index is the number of 3's
* the second fastest changin is the number of 2's,
* the not-at all fast changing is the number of 1's
* (the number before the a's gives the tensor order)
*
* N is the maximum order here
*/
mytype* add3dcoords(int N,mytype *A,const mytype *x,mytype weight)
{
	mytype a=weight,b,c;

	for(int i=N;i>=0;--i,a*=x[0])
	{
		b=a;
		for(int j=i;j>=0;--j,b*=x[1])
		{
			c=b;
			for(int k=j;k>=0;--k,c*=x[2])
			{
				*A += c;
				++A;			
			}		
		}			
	}

	return A;
}



inline int createadd3dcoords(ostream&out,int N,const char *A,const char*x,bool add,const char*w=0)
{
	out<<"{\n"
		" mytype a="<<(w?w:"1")<<",b,c;";
	
	int ret=0;
	mytype a=1,b,c;

	//index in A
	int reti=0;
	for(int i=N;i>=0;--i,a*=2)
	{
		b=a;
		out<<"\n b=a;";
		for(int j=i;j>=0;--j,b*=2)
		{
			c=b;
			out<<"\n c=b;";
			for(int k=j;k>=0;--k,c*=2)
			{
				out<<A<<"["<<reti<<"]"<<(add?"+=c;":"=c;");++reti;
				if(k>0)
				{
					if( c>1|| w!=0){out<<"c*="<<x<<"[2];";++ret;}
					else out<<"c="<<x<<"[2];";
				}
			}		
			if(j>0)
			{
				if( b>1 || w!=0 ){out<<"b*="<<x<<"[1];";++ret;}
				else out<<"b="<<x<<"[1];";
			}
		}		

		if(i>0)
		{
			if( a>1 || w!=0){out<<"a*="<<x<<"[0];";++ret;}
			else out<<"a="<<x<<"[0];";
		}
	}
	
	out<<"\n}\n";
	return ret;
}






void createTotalFold(ostream &out, int N, int ord, const string& T)
{
	out<<T<<" evalTotalFold<"<<T<<","<<ord<<","<<N
		<<">(const "<<T<<"*A,const "<<T<<"* x){"<<endl;
	out<<T<<" ret,px,py;"<<endl;	

	out<<
		T<<" z["<<ord<<"];\n"
		"z[0]= x[2];\n";
	for(int i=1;i<ord;++i)
		out<<"z["<<i<<"]=z["<<i-1<<"]*x[2];\n";

	out<<"ret = z["<<ord-1<<"]* A["<<calcIndex(N,0,0,ord)<<"];\n\n";
	
	int ordoverny=1;
	
	out<<"py = x[1];\n";
	for(int ny=1;ny<ord;++ny)
	{
			ordoverny = ordoverny * (ord-ny+1) /ny;
			out<<"ret+= "<<ordoverny<<" * py * z["<<ord-ny-1<<"] *A["
				<<calcIndex(N,0,ny,ord-ny)<<"];\n";
			out<<"py*=x[1];\n";
	}	
	out<<"ret+= py * A["<<calcIndex(N,0,ord,0)<<"];\n\n";
	
	out<<"px = x[0];\n";
	int ordovernx=1;
	for(int nx=1;nx<ord;nx++)
	{
		ordovernx = ordovernx * (ord-nx+1) /nx;
		int ordovernxny=ordovernx;
		out<<"py=px;\n";
		for(int ny=0;ny<ord-nx;++ny)
		{
			if(ny>0)ordovernxny = ordovernxny * ( ord-nx-ny+1 )/ ny;
			out<<"ret+=	py * z["<<ord-nx-ny-1<<"]";
			if(ordovernxny>1)out<<" * "<<ordovernxny;
			out	<<"* A["<<calcIndex(N,nx,ny,ord-nx-ny)<<"];\n";			
			out<<"py*=x[1];\n";
		}
		out<<"ret+= py * "<<ordovernx<<" * A["<<calcIndex(N,nx,ord-nx,0)<<"];\n";
		out<<"px*=x[0];\n";
	}
	out<<"ret+= px*A["<<calcIndex(N,ord,0,0)<<"];\n";
	out<<"return ret;\n}"<<endl;
}






void createadd3dcoords2s(ostream& out)
{
	for(int i=1;i<=4;++i)
	{
		out<<"\nadd3dcoords2<"<<i<<">::add3dcoords2(mytype *A,const mytype *x)";
		createadd3dcoords(out,i,"A","x",true);	
	}
}

//partial specializations for add3dcoords2s created by createadd3dcoords2s
#include "computeTsetPartialSpecs.h"

inline void testcalcIndex()
{
	int N=20;

	int num=getAsiz(N);
	tsetiter it(N);
	for(int i=0;i<num;++i,++it)
	{
		int o[3];
		it.getOrders(o);
		int calci=calcIndex(N,o[0],o[1],o[2]);
		cout<<o[0]<<','<<o[1]<<','<<o[2]<<':'
			<<calci<<" =? "<< i<<endl;

		if(calci!=i){cout<<"error"<<endl;break;}


	}

}





/** create an optimized specialization of translate<N>
*in file out
*\return: number of multiplications needed
*/
int createTranslate(ostream & out, int N)
{
	
	out<<"translate<"<<N<<">::translate(mytype*newA, const mytype*A, const mytype *t)\n";
	out<<"{\n";

	out<<"mytype at["<<getAsiz(N)<<"];"<<endl;
	int ret=createadd3dcoords(out,N,"at","t",false);

	int reti=0;
	for(int I=0;I<=N;++I)
	{
		for(int J=0;J<=N-I;++J)
		{
			for(int K=0;K<=N-I-J;++K)
			{
				out<<"newA["<<reti<<"]=";++reti;

				int nu=1; //(n over k) times (n over i) times (n over j)
				for(int i=0;i<=I;++i){
					if(i>0)
						nu=nu*(I-i+1)/i;
					for(int j=0;j<=J;++j)
					{
						if(j>0)
							nu=nu*(J-j+1)/j;

						for(int k=0;k<=K;++k)
						{
							if(k>0)
								nu= nu*(K-k+1)/k;
						
							out<<" + ";
							if(nu>1)
							{ out<<nu<<"*";++ret; }
							out<<"A["<<calcIndex(N,i,j,k)<<"]";

							if(i<I || j<J || k<K ) //if not at[0], which is always 1: multiply it
							{ out<<"*at["<<calcIndex(N,I-i,J-j,K-k)<<"]";++ret;}
						
						}					
					}							
				}

				out<<";"<<endl;
				assert(nu==1);
			}
		}					
	}
	out<<"}\n";
	return ret;
}

//the specializations reated by createtranslate
#include "translatePartialSpecs.h"


inline void testcreatetranslate()
{
	std::ofstream out("translateit.txt");
	for(int N=1;N<=4;++N)
	{
		createTranslate(out,N);
		int mults=createTranslate(cout,N);
		cout<<"number of mults for order "<<N<<" : "<< mults<<endl<<endl;
	}

	out.close();




}

/** number of mults in translate set with max. tensor order N 
*(jus used for documentation)
*/
inline int summandsTrans3D(int N)
{
	int ret=0;
	//add effort to compute the T tensor set 

	//add effort to compute the new vals from the T's and 

	for(int I=0;I<=N;++I){
		int r0=0;
		for(int J=0;J<=N-I;++J)
		{
			int r1=0;
			for(int K=0;K<=N-I-J;++K)
			{
				r1+=K+1;
			}
			r0+=r1*(J+1);
		}
		ret+= r0*(I+1);
	}

	return ret;

}







inline void testcorrectnessadd3dc()
{
	static const int N=10;
	int numc=getAsiz(N);
	mytype *A=new mytype[numc],*B=new mytype[numc];
	mytype c[]={2,3,5};

	mytype* A2;
	memset(A,0,numc*sizeof(*A));
	A2=add3dcoords(N,A,c);
	cout<<(A2-A)<<"==?"<<numc<<endl;

	memcpy(B,A,(numc)*sizeof(*A));
	
	memset(A,0,(numc)*sizeof(*A));
	add3dcoords2<N>(A,c);

	tsetiter it(N);
	for(int i=0;i<numc;++i,++it)
	{
		int o[3];
		it.getOrders(o);
		cout<<A[i]<<" =? "<<B[i]<<" =?";
		if(o[0])cout<<" 2";
		if(o[0]>1)cout<<"^"<<o[0];
		if(o[1])cout<<" 3";
		if(o[1]>1)cout<<"^"<<o[1];
		if(o[2])cout<<" 5";
		if(o[2]>1)cout<<"^"<<o[2];
		cout<<endl;
	}
	
	assert(A2-A==numc);

}



inline void testspeedadd3dc()
{
	int num=2*1000*1000;
	mytype *coords=new mytype[num*3];
	static const int N=4;

	int numc=getAsiz(N);
	mytype *A=new mytype[numc];
	
	clock_t st=clock();
	for(int j=1;j<=10;j++)
	{
		mytype*c=coords,*cend=coords+num*3;
		for(;c!=cend;c+=3){
			add3dcoords2<N>(A,c);
			//add3dcoords(N,A,c);
		}

		double nspt=1e9*(clock()-st)/CLOCKS_PER_SEC/j/num;
		cout<<"time per pt :"<< nspt	
		<<" ns time per component:"<<nspt/numc<<" ns"<<endl;
	}

}


inline void testspeedtranslate()
{
	int num=1000*100;
	mytype *coords=new mytype[num*3];
	static const int N=4;

	int numc=getAsiz(N);
	mytype *A=new mytype[numc*num],*B=new mytype[numc*num];
	
	clock_t st=clock();
	for(int j=1;j<=1000;j++)
	{
		mytype*c=coords,*cend=coords+num*3;
		mytype* Ai=A,*Bi=B;
		for(;c!=cend;c+=3,Ai+=numc,Bi+=numc){
			translate<N>(Ai,Bi,c);
		}

		double nspt=1e9*(clock()-st)/CLOCKS_PER_SEC/j/num;
		cout<<"time per pt :"<< nspt	
		<<" ns time per component:"<<nspt/numc<<" ns"<<endl;
	}

}

inline int effortTrans3D(int N){
return summandsTrans3D(N)*2 + getAsiz(N);
}

inline void  printEffortTranslate3D()
{
	int n=10;
	for(int i=0;i<n;++i)
		cout<<i<<" & ";
	cout<<n<<"\\\\"<<endl;

	for(int i=0;i<n;++i)
		cout<<effortTrans3D(i)<<" & ";

	cout<<effortTrans3D(n)<<"\\\\"<<endl;

}



void createBasisTransform(ostream&o,int N)
{
	
	int siz=getAsiz(N);
	polyn * A=new polyn[siz];
	polyn * B=new polyn[siz];
	polyn M[3][3];

	vector<string> paramnames(siz+9);

	ostringstream s;
	for(int i=0;i<siz;++i)
	{	
		B[i]=polyn(1,i);
		s.str("");s<<"B["<<i<<"]";
		paramnames[i] = s.str();
	}

	for(int i=0;i<9;++i)
	{
		M[0][i]=polyn(1,siz+i);
		s.str("");s<<"M["<<i/3<<"]["<<i%3<<"]";
		paramnames[siz+i] = s.str();
	}

	basisTransform<polyn>(N,A,B,M);
	
	for(int i=0;i<siz;++i)
	{
		o<<"A["<<i<<"]=";
		optimizedPoly op(A[i]);
		o<<op;
		//A[i].printOptCommonCoeff(o,paramnames);
		o<<";\n";
	}

	delete[]A;delete[]B;
}



void testadd3dc()
{
	/*
	ofstream out("specscreateadd3d.txt");
	createadd3dcoords2s(out);
	out.close();
	*/
	//testcreatetranslate();
	//printEffortTranslate3D();
	//testcalcIndex();

	//testcorrectnessadd3dc();

	//testspeedtranslate();
	testspeedadd3dc();


}

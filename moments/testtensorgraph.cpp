/*
	test code for tensor graphs

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

#include<algorithm>
#define _USE_MATH_DEFINES
#include<math.h>
#include <time.h>
#include<fstream>

#include "tensWithSymm.h"
#include "binaryTreesVisitor.h"
#include "tensorGraph.h"
#include "tensorGraphSplitTree.h"
#include "indepPolySet.h"
#include "createSymm3DTens.h"
#include "optimizedFold.h"
#include "optimizedPoly.h"
#include "computeMoments.h"
#include "moments.h"
#include "tensorSet.h"
#include "momStats.h"

/**************************************/
/* the testin functions:*/



#define printme(x)\
	cout<<#x<<":"<<(x)<<endl;




inline double numbinb(int n)
{
	double ret=0;
	if(n<=2)
		return 1;

	//restrict to more balanced trees
	int bord=1;
	for(int i=bord;i<=n/2;++i)
	{
		ret += nuk(n,i)*numbinb(i)*numbinb(n-i);
	}
	return ret;
}

inline void testTGsplitreevisitor()
{

	struct visitor:public binaryTreesVisitor
	{
		int num;
		void accept( const binaryTreesVisitor::node*root) {++num;}		
	} v;

	for(int i=2;i<=10;++i)
	{
		v.num=0;v.visit(i);

		cout<<"num split trees for "<<i<<" nodes: \n"
			<<"iterator says:"<<v.num<<endl;
		cout<<"numbinb says"<<numbinb(i)<<endl;

	}
}

inline void testconncompspeed()
{
	unsigned x[6]=
	{	001,006,016,074,070,070 };
	unsigned x2[4]=
	{ 003,003,014,015};

	vector<unsigned> xx(x,x+6);
	vector<unsigned> xx2(32,0xffffffff);
	vector<unsigned> xx3(x2,x2+4);
	vector<unsigned> grp;


	int numper=1000;int numges;
	clock_t start;
	start=clock();	
	for(numges=0;clock()-start<2*CLOCKS_PER_SEC ; numges+=numper)
		for(int i=0;i<numper;++i)
			tensorGraph::connectedComps(xx,&grp);

	cout<<"time per exec:"<< (clock()-start)*1e6/numges<<"mus"<<endl;

	start=clock();	
	for(numges=0;clock()-start<2*CLOCKS_PER_SEC ; numges+=numper)
		for(int i=0;i<numper;++i)
			tensorGraph::connectedComps(xx3);

	cout<<"time per exec:"<< (clock()-start)*1e6/numges<<"mus"<<endl;


}

inline void testIsConnectedSpeed(const tensorGraph&tg)
{
	int numper=1000;int numges;
	clock_t start;
	bool ret;

	start=clock();	
	for(numges=0;clock()-start<2*CLOCKS_PER_SEC ; numges+=numper)
		for(int i=0;i<numper;++i)
			ret=tg.connectedComps();

	cout<<"time per exec: isConnected"<< (clock()-start)*1e6/numges<<"mus"<<endl;
	cout<<"result was:"<<ret<<endl;

	vector< unsigned > vret;

	start=clock();	
	for(numges=0;clock()-start<2*CLOCKS_PER_SEC ; numges+=numper)
		for(int i=0;i<numper;++i)
			ret=tg.connectedComps(&vret);

	cout<<"time per exec: isConnected"<< (clock()-start)*1e6/numges<<"mus"<<endl;
	cout<<"result was:"<<ret<<endl;


	start=clock();	
	vector<unsigned> mat;
	tg.getConnMatrix(mat);
	
	for(numges=0;clock()-start<2*CLOCKS_PER_SEC ; numges+=numper)
		for(int i=0;i<numper;++i)
		{
			ret=tg.connectedComps(mat);
		}

	cout<<"time per exec: connectedComps"<< (clock()-start)*1e6/numges<<"mus"<<endl;
	cout<<"result was:"<<ret<<endl;
}


void printGraphMLheader(ostream&out)
{
	out<<""
"<?xml version=\"1.0\" encoding=\"UTF-8\"?>                 \n"
"<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\"   \n"
"xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"    \n"
"xsi:schemaLocation=\"http://graphml.graphdrawing.org/xmlns "
"http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd\">   \n";	
}

void createTensorGraph()
{
	ifstream mif("momvecraw.txt");

	indepPolySet::tpolyset vecpolys;
	
	vecpolys.clear();	
	indepPolySet::load(mif,vecpolys);


	indepPolySet::tpolyset::iterator it;

	vector< vector<int> > gsizpertype(6);
	for(int i=1;i<6;++i)
	{
		vector<int> & g=gsizpertype[i];
		if(i>1)
			g.push_back(i-1);
		g.push_back(1);
	}


	vector<const repeatTensor*> tensPerType(6);
	vector<tensor<polyn> > tpt(6);
	for(int i=1;i<=3;++i)
	{
		tensorWithSymm
    p = createSymm3DTens2(i-1,2,1);
		tensPerType[i] = new repeatTensor(p);			
		tpt[i] = *tensPerType[i];
	}


	tensorGraphSplitTree2 tgs2(gsizpertype,3);

	vector<polyn> polys2(vecpolys.size());
	int vind=0;

	for(it=vecpolys.begin();it!=vecpolys.end();++it,++vind)
	{
		cout<<*it->mi<<endl;	

		tensmetainfo&tmi = *dynamic_cast<tensmetainfo*>(it->mi);

		tensorGraph tg;
		
		tg.nodelabels = tmi.sourceTensors;
		//the endpoint label for every index of the tensorproduct
		vector<int> endpointlabels;
		//the tensor id for every index of the tensorproduct
		vector<int> tensorids;
		
		for(int i=0;i<tmi.sourceTensors.size();++i)
		{
			for(int j=0;j<tmi.sourceTensors[i]-1;++j)
				endpointlabels.push_back(0);		

			endpointlabels.push_back( tmi.sourceTensors[i]>1 ? 1: 0 );
			
			for(int j=0;j<tmi.sourceTensors[i];++j)
				tensorids.push_back(i);
		}

		for(int i=0;i<tmi.pairs.size();++i)
		{
			tensorGraph::edge e;
			pr&p=tmi.pairs[i];
			e.d.a= tensorids[p.a];
			e.d.b= tensorids[p.b];
			e.d.u= endpointlabels[p.a];
			e.d.v= endpointlabels[p.b];
			tg.edges.push_back(e);
		}

		tgs2.findOptimalSplit(tg);

		tensor<polyn> buf;
		tg.createPolynTensor(buf,tpt,gsizpertype);
		polys2[vind] = buf[0];
		
		polys2[vind].s.sort();

		polyn diff = *it;

		diff.checkSorted();

		//diff.s.sort();
		diff -= polys2[vind];


		if(diff !=polyn(0) ){
			cout<<" polynomes differ at "<<vind<<":\n"
				<<"new polyn was:"<<polys2[vind]<<endl;
			cout<<"\n old polyn was:\n"
				<<polyn(*it)<<endl;
			cout<<"\n difference was:"
				<<diff<<endl;

			int la=polys2[vind].s.size(),lb=polyn(*it).s.size();
			cout<<"lengths:"<<la<<"  "<<lb
				<<endl;

			if(la + lb < 20)
				break;

		}

		tg.normalize();

		tg.createPolynTensor(buf,tpt,gsizpertype);
		{
		polys2[vind] = buf[0];
		
		polys2[vind].s.sort();

		polyn diff = *it;

		diff.checkSorted();

		//diff.s.sort();
		diff -= polys2[vind];


		if(diff !=polyn(0) ){
			cout<<"normalized tg: polynomialss differ at "<<vind<<":\n"
				<<"new polyn was:"<<polys2[vind]<<endl;
			cout<<"\n old polyn was:\n"
				<<polyn(*it)<<endl;
			cout<<"\n difference was:"
				<<diff<<endl;

			int la=polys2[vind].s.size(),lb=polyn(*it).s.size();
			cout<<"lengths:"<<la<<"  "<<lb
				<<endl;

			if(la + lb < 20)
				break;

		}
		}

	}



	tgs2.updateNodeIds();

	ofstream of("vectorhypergraph_1.ps");
	tgs2.printPS(of);	
	of<<"\nshowpage\n";
	of.close();

	ofstream ofc;
	
	ofc.open("vectorTens.txt");
	tgs2.printForC(ofc,tensPerType);
	ofc.close();
	

	
	tgs2.joinSingleParentNodes();
	
	
	of.open("vectorhypergraph_2.ps");
	tgs2.printPS(of);
	of<<"\nshowpage\n";
	of.close();

	ofc.open("vectorTens2.txt");
	tgs2.printForC(ofc,tensPerType);
	ofc.close();

	ofc.open("vectorhypergraph_2.xml");
	printGraphMLheader(ofc);
	tgs2.printGraphML(ofc);
	ofc<<"</graphml>";
	ofc.close();

}





void createOptimizedCommonSub(const indepPolySet::tpolyset& vecpolys,const char*name)
{

	vector<optimizedPoly> optHorner;
	vector<metainfo*> mi;
	
	indepPolySet::tpolyset::const_iterator it;
	
	for(it=vecpolys.begin();it!=vecpolys.end();++it)
	{
		polyn p(*it);
		optHorner.push_back(p);
		mi.push_back(it->mi);
	}

	optimizedPolyCommonSubMulti optCommonSubMult(optHorner);

	ofstream of(name);
	optCommonSubMult.print(of);	

}

void createScalTensorGraph()
{
	ifstream mif("momScal_raw.txt");

	indepPolySet::tpolyset vecpolys;
	
	vecpolys.clear();	
	indepPolySet::load(mif,vecpolys);


	indepPolySet::tpolyset::iterator it;

	vector< vector<int> > gsizpertype(6);
	for(int i=0;i<6;++i)
	{
		vector<int> & g=gsizpertype[i];
		g.push_back(i);
	}


	vector<const repeatTensor*> tensPerType(6);
	vector<tensor<polyn> > tpt(6);
	for(int i=1;i<=4;++i)
	{
		tensorWithSymm p (createSymm3DTens(i,4) );
		tensPerType[i] = new repeatTensor(p);			
		tpt[i] = *tensPerType[i];
	}


	tensorGraphSplitTree2 tgs2(gsizpertype,3);

	vector<polyn> polys2(vecpolys.size());
	int vind=0;


	createOptimizedCommonSub(vecpolys,"momscal_optcommonsub.txt");


	for(it=vecpolys.begin();it!=vecpolys.end();++it,++vind)
	{
		cout<<*it->mi<<endl;	

		tensmetainfo&tmi = *dynamic_cast<tensmetainfo*>(it->mi);

		tensorGraph tg;
		
		tg.nodelabels = tmi.sourceTensors;
		//the endpoint label for every index of the tensorproduct
		vector<int> endpointlabels;
		//the tensor id for every index of the tensorproduct
		vector<int> tensorids;
		
		for(int i=0;i<tmi.sourceTensors.size();++i)
		{
			for(int j=0;j<tmi.sourceTensors[i];++j)
				endpointlabels.push_back(0);		
			
			for(int j=0;j<tmi.sourceTensors[i];++j)
				tensorids.push_back(i);
		}

		for(int i=0;i<tmi.pairs.size();++i)
		{
			tensorGraph::edge e;
			pr&p=tmi.pairs[i];
			e.d.a= tensorids[p.a];
			e.d.b= tensorids[p.b];
			e.d.u= endpointlabels[p.a];
			e.d.v= endpointlabels[p.b];
			tg.edges.push_back(e);
		}

		tgs2.findOptimalSplit(tg);

		tensor<polyn> buf;
		tg.createPolynTensor(buf,tpt,gsizpertype);
		polys2[vind] = buf[0];
		
		polys2[vind].s.sort();

		polyn diff = *it;

		diff.checkSorted();

		//diff.s.sort();
		diff -= polys2[vind];


		if(diff !=polyn(0) ){
			cout<<" polynomes differ at "<<vind<<":\n"
				<<"new polyn was:"<<polys2[vind]<<endl;
			cout<<"\n old polyn was:\n"
				<<polyn(*it)<<endl;
			cout<<"\n difference was:"
				<<diff<<endl;

			int la=polys2[vind].s.size(),lb=polyn(*it).s.size();
			cout<<"lengths:"<<la<<"  "<<lb
				<<endl;

			if(la + lb < 20)
				break;

		}

		tg.normalize();

		tg.createPolynTensor(buf,tpt,gsizpertype);
		{
		polys2[vind] = buf[0];
		
		polys2[vind].s.sort();

		polyn diff = *it;

		diff.checkSorted();

		//diff.s.sort();
		diff -= polys2[vind];


		if(diff !=polyn(0) ){
			cout<<"normalized tg: polynomialss differ at "<<vind<<":\n"
				<<"new polyn was:"<<polys2[vind]<<endl;
			cout<<"\n old polyn was:\n"
				<<polyn(*it)<<endl;
			cout<<"\n difference was:"
				<<diff<<endl;

			int la=polys2[vind].s.size(),lb=polyn(*it).s.size();
			cout<<"lengths:"<<la<<"  "<<lb
				<<endl;

			if(la + lb < 20)
				break;

		}
		}

	}



	tgs2.updateNodeIds();

	ofstream of("scalarhypergraph_1.ps");
	tgs2.printPS(of);	
	of<<"\nshowpage\n";
	of.close();

	ofstream ofc;
	
	ofc.open("scalarTens.txt");
	tgs2.printForC(ofc,tensPerType);
	ofc.close();
	

	
	tgs2.joinSingleParentNodes();
	
	
	of.open("scalarhypergraph_2.ps");
	tgs2.printPS(of);
	of<<"\nshowpage\n";
	of.close();

	ofc.open("scalarTens2.txt");
	tgs2.printForC(ofc,tensPerType);
	ofc.close();

	ofc.open("scalarTens2.xml");
	printGraphMLheader(ofc);
	tgs2.printGraphML(ofc);
	ofc<<"</graphml>\n";
	ofc.close();


}





//write the moments saved in momvecraw.txt
//as a tex-file and as single postcript files.
void writeVecMoments2Tex()
{
	ifstream mif("momvecraw.txt");

	indepPolySet::tpolyset vecpolys;
	indepPolySet::tpolyset::const_iterator it;
	vecpolys.clear();	
	indepPolySet::load(mif,vecpolys);

	vector<vector<int> > groups(20);
	for(int i=1;i<groups.size();++i)
	{
		vector<int>&g=groups[i];
		if(i>1)
			g.push_back(i-1);
		g.push_back(1);
	}

	ofstream out("vecmomstex.tex");
	ofstream out2("vecmomstex2.tex");

	cout<<"writing texfile vecmomstex.tex"<<endl;;
	int i=0;
	for(it=vecpolys.begin();it!=vecpolys.end();++it,++i)
	{
		cout<<*it->mi<<endl;	
		
		tensmetainfo&tmi = *dynamic_cast<tensmetainfo*>(it->mi);

		for(int j=0;j<tmi.sourceTensors.size();++j)
		{
			vector<int>&g=groups[ tmi.sourceTensors[j] ];
			tmi.tgrp.push_back(g.size());
			for(int k=0;k<g.size();++k)
				tmi.grp.push_back(g[k]);
		}

		out<<"% "<<i<<":\n";
		tmi.printPairsProdForTex(out);
		out<<"&=&";
		tmi.printEsumForTex(out);
		out<<"\\\\%\n%\n";
		out2<<"$M_{"<<i<<"}=";
		tmi.printEsumForTex(out2);
		out2<<"$ , %\n";
		

		ostringstream nam;nam<<"vecmom\\vecmom"<<i<<".eps";
		ofstream pso(nam.str().c_str());
		tmi.printGraphPS(pso);
	}
}



inline void testMomFunSpeed()
{
	int num= 100*1000;
	static const int N=4;

	int numc=getAsiz(N);

	cout<<"number components:"<<numc<<endl;
	int numm=28;
	mytype *A=new mytype[numc*num];
	mytype *M=new mytype[numm*num];

	for(int i=0;i<numc*num;++i)
	{
		A[i]=rand();
	}

	
	clock_t st=clock();
	for(int j=1;j<=100&&(clock()-st) < 2*CLOCKS_PER_SEC ;j++)
	{
		mytype* Ai=A,*Aend= A+numc*num,*Mi=M;
		for(;Ai!=Aend;Ai+=numc,Mi+=numm){
			computeMomenteDim3Ord4_opt(Mi,Ai);
		}

		double nspt=1e9*(clock()-st)/CLOCKS_PER_SEC/j/num;
		cout<<"time per tset :"<< nspt	
		<<" ns time per component:"<<nspt/numc<<" ns"<<endl;
	}


	st=clock();
	for(int j=1;j<=100&&(clock()-st) < 2*CLOCKS_PER_SEC ;j++)
	{
		mytype* Ai=A,*Aend= A+numc*num,*Mi=M;
		for(;Ai!=Aend;Ai+=numc,Mi+=numm){
			computeMomenteDim3Ord4_nopt(Mi,Ai);
		}

		double nspt=1e9*(clock()-st)/CLOCKS_PER_SEC/j/num;
		cout<<"not opt: time per tset :"<< nspt	
		<<" ns time per component:"<<nspt/numc<<" ns"<<endl;
	}

}


inline void testperm2indperm()
{

	int gsiz[3]={1,2,3};
	//index groups:
	//0, 1,2, 3,4,5

	int perm[3]={2,0,1};

	vector<int> indperm;

	tensorGraphSplitTree2::perm2indperm
		(indperm,vector<int>(perm,perm+3),
			vector<int>(gsiz,gsiz+3)
			);


}


inline void testmomsetvectorscorrect()
{
	polyn A[30];
	
	static const int nv=3;
	polyn M[nv][27];

	momentSet<polyn,27,0,0> m0;
	momentSet<polyn,27,0,1> m1;
	momentSet<polyn,27,0,2> m2;
	//momentSet<polyn,27,0,3> m3;
	//momentSet<polyn,27,0,4> m4;


	const char*comments[nv]={m0.comment,
		m1.comment,m2.comment//,m3.comment,m4.comment
	};

	for(int i=0;i<30;++i)
		A[i] = polyn(1,i);

	m0.compute(M[0],A);
	m1.compute(M[1],A);
	m2.compute(M[2],A);
	//m3.compute(M[3],A);
	//m4.compute(M[4],A);


	for(int i=0;i<nv;++i)
		cout<<i<<":\n"<<comments[i]<<endl<<endl;

	for(int i=0;i<nv;++i)
	{
		int k=(i+1)%nv;
		cout<<"\n"<<i<<"  vs.  "<<k<<":\n\n";
		int j;
		for(j=0;j<27;++j)
		{
			cout<<"moment number"<<j<<":";
			if(M[i][j]!=M[k][j])
			{

				cout<<"wrong. details:\n"
					<< m0.metainfo.me[j].des<<"\n"<<endl;


				cout<<M[i][j]<<"\n\n\n"<<M[k][j]
						<<"\n\n\n"<<M[i][j]-M[k][j]<<"\n\n\n";

			}
			else
				cout<<"o.k."<<endl;
		}

	}
}

inline void createspecs()
{
	ofstream of;
	of.open("vecspecs_opt0.txt");
	createVecMomSpecialization(of,3,false,"unoptimized");
	of.close();
	of.open("vecspecs_opt_comonsubex.txt");
	createVecMomSpecialization(of,1,true,
		"optimized for common subexpressions");
	of.close();
}


inline void testprintpolyn()
{
	polynom pp("a+bbc+bbc+d+d+c");
	polyn p(pp);

	const string pname[]={"a","b","c","d"};
	
	cout<<p<<endl;

	p.printOptCommonCoeff(cout,vector<string>(pname,pname+4));



}



inline void testRotInvV()
{
	
	tensorSet3D<2,mynum> A[3],B[3],C[3];

	for(int i=0;i<3;++i)
		for(int j=0;j<A[i].ncomps;++j)
			A[i].A[j] = -1 + 2.0*rand()/RAND_MAX;
	
	mynum R[3][3];
	createRandRot(R);
	mynum S[3][3];
	for(int i=0;i<3;++i)for(int j=0;j<3;++j)
		S[i][j]=R[j][i];

	//do it:
	basisTransformV(2,B[0].A,A[0].A,R);
	
	//undo it:
	basisTransformV(2,C[0].A,B[0].A,S);

	/*
	//test it:
	cout<<"test basis transform:"<<endl;
	for(int i=0;i<3;++i)
		for(int j=0;j<A[i].ncomps;++j)
		{
			cout<<A[i].A[j]<<" = "<<C[i].A[j]<< " != "<<B[i].A[j]<<endl;
		}

	*/

	cout<<"test moments:"<<endl;
	momentSet<mynum,27,0,0> m10,m20;
	momentSet<mynum,27,0,1> m11,m21;
	momentSet<mynum,27,0,2> m12,m22;

	m10.compute(A[0].A);m20.compute(B[0].A);
	m11.compute(A[0].A);m21.compute(B[0].A);
	m12.compute(A[0].A);m22.compute(B[0].A);

	momentSetBase<mynum,27,0> 
		m1[3]={m10,m11,m12},
		m2[3]={m20,m21,m22};

	bool rotinv=true;
	for(int i=0;i<m1[0].n;++i)
	{
		for(int v=0;v<3;++v)
		{
			if(fabs(m1[v].M[i]-m2[v].M[i]) > 
				1e-9* min(fabs(m1[v].M[i]),fabs(m2[v].M[i]) ) )
			{
				cout<<" moment "<<i<<", version" <<v<<": not rot.inv:"
					<<  m1[v].M[i]<<" != "<<m2[v].M[i]<<endl;
				rotinv=false;
			}

			int w=(v+1)%3;
			if( v>0 && fabs(m1[v].M[i]-m1[w].M[i]) > 
				1e-9*min( fabs(m1[v].M[i]) ,fabs(m1[w].M[i]) ) )
			{
				cout.precision(20);
				cout<<" moment "<<i<<": version " <<v<<" != version "<<w
					<<" : "<<  m1[v].M[i]<<" != "<<m1[w].M[i]<<endl
					<<"difference:"<<m1[v].M[i]-m1[w].M[i]<<endl;
			}
		}
	}

	if(rotinv)cout<<"all moment invs were rotationally invariant"<<endl;




}

inline void testRotInvN()
{
	momentSet<mynum,28,0,0> m1,m2;
	
	tensorSet3D<4,mynum> A,B,C;

	for(int j=0;j<A.ncomps;++j)
			A.A[j] = -1 + 2.0*rand()/RAND_MAX;
	
	mynum R[3][3];
	createRandRot(R);
	mynum S[3][3];
	for(int i=0;i<3;++i)for(int j=0;j<3;++j)
		S[i][j]=R[j][i];

	//do it:
	basisTransform(4,B.A,A.A,R);
	
	//undo it:
	basisTransform(4,C.A,B.A,S);

	//test it:
	/*
	cout<<"test basis transform:"<<endl;
		for(int j=0;j<A.ncomps;++j)
		{
			cout<<A.A[j]<<" = "<<C.A[j]<< " != "<<B.A[j]<<endl;
		}
	*/



	cout<<"test moments rot.inv., (scalar moments)" <<endl;

	m1.compute(A.A);m2.compute(B.A);

	for(int i=0;i<m1.n;++i)
		if(fabs(m1.M[i]-m2.M[i]) > 1e-6*fabs(m1.M[i]))
		{
			cout<<i<<":"<<  m1.M[i]<<" != "<<m2.M[i]<<endl;
		}

	cout<<"all moments tested"<<endl;


}



void testtensorgraph()
{
/*
	testRotInvN();
	testRotInvV();
	*/

	//writeVecMoments2Tex();
	
	//testprintpolyn();
	//createspecs();

	//testmomsetvectorscorrect();


  createScalTensorGraph();
  createTensorGraph();


	/*
	testMomFunSpeed2<30,momentSet<long double,27,0,0> >();
	testMomFunSpeed2<30,momentSet<long double,27,0,1> >();
	testMomFunSpeed2<30,momentSet<long double,27,0,2> >();


	testMomFunSpeed2<35,momentSet<long double,28,0,0> >();
	testMomFunSpeed2<35,momentSet<long double,28,0,1> >();
	testMomFunSpeed2<35,momentSet<long double,28,0,2> >();
	testMomFunSpeed2<35,momentSet<long double,28,0,3> >();
	*/

	//testperm2indperm();


	//testTGsplitreevisitor();

	//testconncompspeed();

#if 0
	tensor<polynom> t(4,3);
	tensorWithSymm s(t);
	s.grp.resize(2,2);
	s.tgrp.resize(1,2);

	tensProdWithSymm tpws(s,s),tpws2(tpws,s);

	pr x1[]={{0,4},{1,3},{2,6}};
	pr x2[]={{0,5},{1,3},{2,7}};
	pr x3[]={{0,1},{2,4},{5,8},{9,10},{3,6},{7,11}};

	refp<tensorGraph> tg1(new tensorGraph(tpws2,vpr(x1,x1+3)) );
	tensorGraph tg2(tpws2,vpr(x2,x2+3));
	tensorGraph tg3(tpws2,vpr(x3,x3+6));

	refp<tensorGraph> tt;
	tt=tg1;

	/*
	testIsConnectedSpeed(*tg1);
	testIsConnectedSpeed(tg2);
	testIsConnectedSpeed(tg3);
*/

	tg1->normalize();
	tg2.normalize();
	tg3.normalize();
	printme(*tg1==tg2);
	printme(tg2==tg3);
	printme(tg3==*tg1);

	
	vector<int> orderpernodetype(10);
	
	for(int i=0;i<orderpernodetype.size();++i)
		orderpernodetype[i]=i;


	vector< vector<int> > groupSizPerType(10);
	groupSizPerType[4].resize(2,2);

	tensorGraphSplitTree tgs;
	size_t eff= tgs.findOptimalSplit(tg3,orderpernodetype,3);
cout<<"effort:"<<eff<< endl;

	ofstream tg3ps("tg3.ps");
	/*
	tg1->printPS(tg3ps,groupSizPerType);
	tg2.printPS(tg3ps,groupSizPerType);
	tg3.printPS(tg3ps,groupSizPerType);
	*/

	
	tgs.printPS(tg3ps,groupSizPerType);
	

	tensorGraphSplitTree2 tgs2(groupSizPerType,3);
	eff=tgs2.findOptimalSplit(tg3);
	cout<<"effort:"<<eff<< endl;

	tg3ps<<"showpage\n";
	tgs2.printPS(tg3ps);


	tg3ps.close();



	
#endif
	
}

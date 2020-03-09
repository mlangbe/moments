/*
	analytic representation of the integral over a trilinearly interpolated value in an arbitrarily-shaped hexaeder.

	Copyright (C) 2020 TU Kaiserslautern, Prof.Dr. Hans Hagen (AG Computergraphik und HCI) (hagen@cs.uni-kl.de) et al.

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
#include<fstream>
#include<stdlib.h>
#include<time.h>
#include"polynom.h"
#include"poly.h"
#include "optimizedPoly.h"
#include"bruch.h"
#include<math.h>

typedef bruch<int,0x3fffffff> mybruch;
typedef poly<36,mybruch> mypoly;




//create interpolation function
//with the indices 0..7 are the vertices and 8..10 are the local coords
polyn interpolateHexaeder()
{
	polyn ret=0;

	for(int i=0;i<8;++i)
	{
		polyn x(1,i);
		for(int j=0;j<3;++j)
			if(i&(1<<j))
				x*=polyn(1,8+j);
			else
				x*=((polyn)1)-polyn(1,8+j);
		ret+=x;
	}

	return ret;
}

/*
struct coeffsmallerthan
{	
	double bord;
	coeffsmallerthan(double x=1e-14):bord(x){}

	inline bool operator() ( const mypoly::summand&s )
	{
		return false;//-bord< s.c&&s.c<bord;
	}
};
*/
void cleanmypoly(mypoly &p,double bord)
{
	/*
	p.s.resize( 
		remove_if(p.s.begin(),p.s.end(),coeffsmallerthan(bord))-p.s.begin()
		);
		*/
}



/**integrate scalar field s folded with function f
*of the space
* express 
* \[
*\int_0^1\int_0^1\int_0^1 
*det(d\vec p/d\vec l) s(\vec l) f(\vec p(\vec l))  dl_0dl_1dl_2
* as polynome in the vertex values of a hexaeder cell of the field
\]
*/
inline polyn integrateHexaeder
(  
 /*
	const int vertids[8][3],const int valids[8],
	//polynomial of var ids 0,1,2
	*/
	const polyn& f,
	bool unitcube=false
)
{

	polyn ret;

		polyn matrix[3][3];
			
		polyn l[3];

		polyn verts[8][3];
		polyn vals[8];

		int currid=0;

		vector<polyn> bel0(24+8+3);

		for(int i=0;i<8;++i)
			for(int j=0;j<3;++j,++currid)
				bel0[currid]=verts[i][j]
					= ( unitcube?  polyn( (i>>j)&1 ) : polyn(1,currid) );

		int valstart=currid;
		for(int i=0;i<8;++i,++currid)
			bel0[currid]=vals[i]=polyn(1,currid);

		int lstart=currid;
		for(int i=0;i<3;++i,++currid)
			bel0[currid] = l[i] = polyn(1,currid);

		polyn m[3][3];
		
		polyn trilinear = interpolateHexaeder();

		cout<<"calculating determinant ..."<<endl;	
		//global coords
		vector<polyn> p(3);

		vector<polyn> bel(11);		
		for(int i=0;i<3;++i)
		{
			for(int j=0;j<8;++j)
				bel[j]=verts[j][i];
			copy(l,l+3,bel.begin()+8);
			p[i]=trilinear.eval(bel);
			cout<<"p["<<i<<"]="<<endl<<mypoly(p[i])<<endl;
			for(int j=0;j<3;++j)
			{
				m[i][j] = trilinear.evalderiv(bel,8+j);
				cout<<"\nm"<<i<<j<<"="<<endl<<mypoly(m[i][j])<<endl;
			}
		}

		polyn det;

		//m2 = m0 % m1
		for(int i=0;i<3;++i)
		{
			int j=(i+1)%3,k=(i+2)%3;
			det+= m[2][k] * ( m[0][i]*m[1][j] -m[1][i]*m[0][j] );
		}

		//cout<<det;


		cout<<"calculating polynom s ..."<<endl;	
		//calculate field s
		copy(vals,vals+8,bel.begin());
		copy(l,l+3,bel.begin()+8);
		polyn s = trilinear.eval(bel);


		cout<<"calculating polynom ff ..."<<flush;	
		//calculate function f over the field
		polyn ff= f.eval(p);

		cout<<"ready"<<endl;

		double nff=ff.s.size(),ns=s.s.size(),ndet=det.s.size();
		cout<<"number of summands:"
			<<" ff:"<<nff<<" det:"<<ndet<<" s:"<<ns<< endl;	
		double prod=nff*ns*ndet;

		if(prod>100e6) return 0;
		cout<<"putting together the functions:(expected maximum number of summands: "
			<<prod<<")"<<flush;

		polyn expr = det * s * ff;
		cout<<"ready."<<endl;

		cout<<"resulting polynomial has "<<expr.s.size()<<" summands"<<endl;

		cout<<"calculating integrals..."<<flush;	

		//integrate over l0
		bel0[lstart]=0;
		expr=expr.evalint<polyn>(bel0,lstart,1);
		
		//integrate over l1
		bel0[lstart+1]=0;
		expr=expr.evalint<polyn>(bel0,lstart+1,1);

		//integrate over l2
		bel0[lstart+2]=0;
		expr=expr.evalint<polyn>(bel0,lstart+2,1);
		cout<<"ready."<<endl;

		cout<<"overall result has"<< expr.s.size()<<" summands "<<endl;
					
		return expr;

}







/**integrate scalar field s folded with function f
*of the space
* express 
* \[
*\int_0^1\int_0^1\int_0^1 
*det(d\vec p/d\vec l) s(\vec l) f(\vec p(\vec l))  dl_0dl_1dl_2
* as polynome in the vertex values of a hexaeder cell of the field
\]
*/
inline mypoly integrateHexaeder
(  

	const mypoly verts[8][3],
	const mypoly vals[8],
	//polynomial of var ids 0,1,2
	const mypoly& f
)
{


		mypoly ret;

		mypoly matrix[3][3];


		//ids of the integration parameters
		int lid[3];

		//find free var ids for integration parameters:
		bool vars[mypoly::n],varsused[mypoly::n];
		std::fill(varsused,varsused+mypoly::n,false);

		for(int i=0;i<8;++i){
			for(int j=0;j<3;++j)
			{
				verts[i][j].getVarIdsPresent(vars);
				for(int k=0;k<mypoly::n;++k)
					varsused[k] |= vars[k];
			}
			vals[i].getVarIdsPresent(vars);
			for(int k=0;k<mypoly::n;++k)
					varsused[k] |= vars[k];
		}

		int ll=0;
		for(int i=0;i<3;++i)
		{
			while(ll<mypoly::n && varsused[ll] )++ll;
			assert(ll<mypoly::n);
			lid[i]=ll;
			++ll;
		}

		cout<<"lids"<<lid[0]<<' '<<lid[1]<<' '<<lid[2]<<endl;



		mypoly m[3][3];
		
		mypoly trilinear = (mypoly)interpolateHexaeder();

		cout<<"calculating determinant ..."<<endl;	
		//global coords
		vector<mypoly> p((size_t)mypoly::n,(mypoly)mypoly::tcoeff());
		vector<mypoly> bel((size_t)mypoly::n,(mypoly)mypoly::tcoeff());
		mypoly buf;
		for(int j=0;j<3;++j)
			bel[j+8]=mypoly(1,lid[j]);

		for(int i=0;i<3;++i)
		{
			for(int j=0;j<8;++j)
				bel[j]=verts[j][i];
			
			trilinear.eval<mypoly>(bel,p[i]);
			cout<<"p["<<i<<"]="<<endl<<p[i]<<endl;

			for(int j=0;j<3;++j)
			{
				m[i][j] = p[i];
				m[i][j].differentiate(lid[j]);
				cout<<"\nm"<<i<<j<<"="<<endl<<m[i][j]<<endl;
			}
		}

		mypoly det;

		//m2 = m0 % m1
		for(int i=0;i<3;++i)
		{
			int j=(i+1)%3,k=(i+2)%3;
			det+= m[2][k] * ( m[0][i]*m[1][j] -m[1][i]*m[0][j] );
		}

		//cout<<"determinant:"<<det;


		cout<<"calculating polynom s ..."<<endl;	
		//calculate field s
		copy(vals,vals+8,bel.begin());
		mypoly s = trilinear.eval(bel);

		cout<<"calculating polynom ff ..."<<flush;	
		//calculate function f over the field
		mypoly ff= f.eval(p);

		cout<<"ready"<<endl;

		double nff=ff.s.size(),ns=s.s.size(),ndet=det.s.size();
		cout<<"number of summands:"
			<<" ff:"<<nff<<" det:"<<ndet<<" s:"<<ns<< endl;	
		double prod=nff*ns*ndet;

		cout<<"putting together the functions:(expected maximum number of summands: "
			<<prod<<")"<<flush;
		prod *=sizeof(mypoly::summand)*1e-6;
		cout<<"expected size: "<<prod<<" MB"<<endl;

		if(prod > 1000) return (mypoly)(mybruch)0;


		mypoly expr = det * s * ff;
		cout<<"ready."<<endl;

		cout<<"resulting polynomial has "<<expr.s.size()<<" summands"<<endl;

		cout<<"calculating integrals..."<<flush;	

		//integrate
		expr.integrate(lid[0],0,1);
		expr.integrate(lid[1],0,1);
		expr.integrate(lid[2],0,1);

		cout<<"ready."<<endl;

		cout<<"overall result has"<< expr.s.size()<<" summands "<<endl;
					
		return expr;

}


mypoly integrateHexaeder
(  
	//polynomial of var ids 0,1,2
	const mypoly& f,bool unitcube=false
)
{
	int currid=0;

	mypoly verts[8][3],vals[8];


	for(int i=0;i<8;++i)
		for(int j=0;j<3;++j,++currid)
			verts[i][j] = unitcube?  
			mypoly( mybruch((i>>j)&1) ) 
			: mypoly(1.0l,currid) ;

	for(int i=0;i<8;++i,++currid)
		vals[i]=mypoly(mypoly::tcoeff(1),currid);
	
	return integrateHexaeder(verts,vals,f);		
}










inline polyn randpolyn()
{
	polyn a=1;
	for(int i=0;i<10;++i)
	{
		polyn r=polyn(2-rand()*5/RAND_MAX,rand()*4/RAND_MAX)+1;		
		a*=r;
	}

	return a;
}


#define E(x) #x ## << ## ":\n" ## << ## mypoly(x)

#define outerr(op)\
if( mypoly(c) != _c ){  \
	cout<<" error at operation:"<< #op <<":\n" \
<<"\n---------------------------------\n" \
<<"_a=\n"<<_a \
<<"\n---------------------------------\n" \
<<"_b=\n"<<_b \
<<"\n---------------------------------\n" \
	<<mypoly(c)<<"\n = c != _c =\n"<<_c \
	<<"\n difference:\n"<<_c-mypoly(c)<<endl;return;\
} 


#define testme(op) \
c= a op b;\
_c= _a op _b;\
	outerr(op) 

#define testme2(op) \
	 c= a;c op ## = b;\
	 _c= _a; _c op ## = _b;\
outerr(op ##= )


inline void testpolys()
{

	vector<polyn> aa,bb;
	vector<mypoly> _aa,_bb;
	mypoly _c;polyn c;

	int n=300;
	for(int i=0;i<n;++i)
	{
		polyn a=randpolyn(),b=randpolyn(),c;
		mypoly _a(a),_b(b),_c;
		aa.push_back(a);bb.push_back(b);
		_aa.push_back(_a);_bb.push_back(_b);

		//cout<<E(a)<<"\n--------------------\n"<<E(b)<<endl;

		//cout<<"------------------------------"<<endl;
		cout<<i<<'\r'<<flush;
		testme(+)
		testme(*)
		testme(-)
	}


	//testing speed
	clock_t tstart;
	int effort;int totn;
	double tim;
	

	tstart=clock();
	effort=0;totn=0;
	while( clock()-tstart < 1000 )
		for(int i=0;i<aa.size();++i)
		{
			c=aa[i]+bb[i];
			effort+=aa[i].s.size()+bb[i].s.size();
			++totn;
		}
	tim=1.0*(clock()-tstart)/CLOCKS_PER_SEC;

	cout<<"speed polyn +:" << tim*1e3/totn << "ms per operation and"
			<< tim*1e9/effort <<" ns per summand"<<endl;

	tstart=clock();
	effort=0;totn=0;
	while(clock()-tstart < 1000 )
		for(int i=0;i<_aa.size();++i)
		{
			_c=_aa[i]+_bb[i];
			effort+=_aa[i].s.size()+_bb[i].s.size();
			++totn;
		}
	tim=1.0*(clock()-tstart)/CLOCKS_PER_SEC;

	cout<<"speed mypoly +:" << tim*1e3/totn << "ms per operation and"
			<< tim*1e9/effort <<" ns per summand"<<endl;

	tstart=clock();
	effort=0;totn=0;
	while( clock()-tstart < 1000 )
		for(int i=0;i<aa.size();++i)
		{
			c=aa[i]*bb[i];
			effort+=aa[i].s.size()*bb[i].s.size();
			++totn;
		}
	tim=1.0*(clock()-tstart)/CLOCKS_PER_SEC;

	cout<<"speed polyn *:" << tim*1e3/totn << "ms per operation and"
			<< tim*1e9/effort <<" ns per summand"<<endl;

	tstart=clock();
	effort=0;totn=0;
	while(clock()-tstart < 1000 )
		for(int i=0;i<_aa.size();++i)
		{
			_c=_aa[i]*_bb[i];
			effort+=_aa[i].s.size()*_bb[i].s.size();
			++totn;
		}
	tim=1.0*(clock()-tstart)/CLOCKS_PER_SEC;

	cout<<"speed mypoly *:" << tim*1e3/totn << "ms per operation and"
			<< tim*1e9/effort <<" ns per summand"<<endl;







}

inline int polykgv(mypoly&p)
{
	int ret=1;
	for(int i=0;i<p.s.size();++i)
		ret=kgv(ret,p.s[i].c.n);
	return ret;
}


/** converte the integration over a hex into one over
*a tetrahedron, a pyramid, aprism*
*/
inline mypoly collapsehex(const int*verts,mypoly hex)
{
	vector<mypoly> bel(mypoly::n);
	int id=0;
	for(int i=0;i<8;++i)
		for(int j=0;j<3;++j,++id)
			bel[id]=mypoly(1,verts[i]*3 + j);

	for(int i=0;i<8;++i,++id)
		bel[id]=mypoly(1,24+i);
	mypoly ret;
	hex.eval(bel,ret);
	

	return ret;

}


typedef poly<36,int> mypolyint;



template<class polyType>
inline double testnumstab(const polyType&p,vector<float>& bel)
{
	
double maxerr=0;


	double f0=p.eval(bel);
	double eps=1e-5;

	//cout<<"original value"<<f0<<endl;
	//cout<<"relative errors if ith val is tainted by "
	//		"a relative error of"<<eps<<":\n"<<endl;

	for(int i=0;i<bel.size();++i)
	{
		double obel=bel[i];
		bel[i]*=1+eps;
		
		double f1 = p.eval(bel);

		double e=(f1-f0)/f0;
		if(maxerr<fabs(e))maxerr=fabs(e);
	
		bel[i]=obel;
	}
	//cout<<endl;


	return maxerr;



}





inline void testoptimize()
{

	mypolyint e(polyn(polynom("aaa+aab+ac+acc+abc")));

	
	ifstream in("testpoly.txt");
	e.load(in);
	in.close();

	ofstream out("optimizedpoly.txt");


	polyn epolyn=e.topolyn();

	if(e.s.size()<20)
		cout<<e<<endl;

	optimizedPoly optp(epolyn);
	out<<"\noptimized:\n"<<optp;
	optimizedPolyCommonSub optp2(e.topolyn());

	out<<"\nfurther optimized:\n"<<optp2;
	out.close();

	vector<mypolyint> testbel(mypoly::n);
	for(int i=0;i<mypoly::n;++i)
		testbel[i]=mypolyint(1,i);


	cout<<"testing num. stabiliy:"<<endl;

	double maxerr[3]={0,0,0};
	for(int iter=0;iter<100;++iter)
	{
		vector<float> bel2(mypoly::n);
		for(int i=0;i<mypoly::n;++i)bel2[i]=rand()-RAND_MAX/2;
		double e;
		e=testnumstab(epolyn,bel2);
		if(e>maxerr[0])maxerr[0]=e;
		e=testnumstab(optp,bel2);
		if(e>maxerr[1])maxerr[1]=e;
		e=testnumstab(optp2,bel2);
		if(e>maxerr[2])maxerr[2]=e;

		cout<<'\r'<<iter<<":max.relative error:"
			<<maxerr[0]<<' '<<maxerr[1]<<' '<<maxerr[2]<<flush;
	}



	cout<<" testing equality optim.Poly common sub:"<<endl;
	mypolyint optp2eval=optp2.eval(testbel);
	mypolyint optpeval=optp.eval(testbel);
	mypolyint epolyneval=epolyn.eval(testbel);
	cout<<" ready."<<endl;
	

#define ifeq(a,b) \
	if(a!=b){\
	cout<< #a << "!= " << #b << endl; \
	if(a.s.size()<20 &&b.s.size()<20) \
	 cout<< #a <<"=\n"<<a<<endl \
		   << #b <<"=\n"<<b<<endl;\
	}


	ifeq(optp2eval,optpeval)
	ifeq(optpeval,e)
	ifeq(e,epolyneval)

#undef ifeq




}



inline void testcalculate()
{

//	testpolys();
	mypoly x(interpolateHexaeder());
	cout<< x<<endl;
	

	polyn f=polynom("a");

	//polyn factor = 2*3*4*5; factor=factor*factor*factor +0.1;
	//factor=1;

	cout<<"weight function:"<<f<<endl;

	mypoly e=integrateHexaeder(mypoly(f));

	std::ofstream integ("integral.txt");
	integ.precision(20);
	integ<<"//integral of"<<f<<"times the field in a hexaeder:"<<endl;
	integ<<"//ids 0..23 are the vertices, ids 24..32 are the vertex values "<<endl;
	integ<<"\n not optimized:\n";

	double factor=polykgv(e);
	e*=mypoly((mybruch)factor);
	
	print(integ,e)<<"/"<<factor<<'\n';



	ofstream ou("testpoly.txt");
	e.save(ou);




	vector<mypoly> bel(mypoly::n);

	int aid=0;
	for(int i=0;i<8;++i)
		for(int j=0;j<3;++j,++aid)
			bel[aid]=mybruch( (i>>j)&1);

	for(int i=0;i<8;++i,++aid)
		bel[aid]=( mypoly(1,24+i) );

	mypoly quader = e.eval(bel);

	integ<<"\n if the hexaeder is a unit cube:\n"
		<<quader<<endl;


	int vertstet[8]={0,0,1,1,2,3,2,3} ;			
	mypoly tet = collapsehex(vertstet,e);
	integ<<"\n if the hexaeder is a tetrahedron:\n"<<tet<<'\n';
	cout<<"tet:"<<tet.s.size()<<"summands\n";


	integ<<"optimized:"<<optimizedPoly(tet.topolyn())<<"/"<<factor<<'\n';
	
	int vertspyram[8]={0,1,2,3,4,4,4,4} ;			
	mypoly pyram = collapsehex(vertspyram,e);
	integ<<"\n if the hexaeder is a pyramid:\n"
		<<pyram.s.size()<<"summands:\n"
		<<pyram<<endl;

	integ<<"optimized:"<<optimizedPoly(pyram.topolyn())<<"/"<<factor<<'\n';

	int vertsprism[8]={0,1,2,2,3,4,5,5} ;			
	mypoly prism = collapsehex(vertsprism,e);
	integ<<"\n if the hexaeder is a pyramid:\n"
		<<prism.s.size()<<"summands:\n"
		<<prism<<endl;

	integ<<"optimized:"<<optimizedPoly(prism.topolyn())<<"/"<<factor<< '\n';



	mypoly xx = (mypoly::tcoeff)1;
	for(int i=0;i<4;++i,xx*=mypoly(1,0))
	{		
		mypoly y=xx;
		for(int j=i;j<4;++j,y*=mypoly(1,1))
		{
			mypoly z=y;
			for(int k=j;k<4;++k,z*=mypoly(1,2))
			{
				cout<<"weight function"<<z<<endl;
				cout << integrateHexaeder(z,true);
					
			
			}
		}
	}

	/*
	polyn ee=integrateHexaeder(f);



	if(!(mypoly(ee)==e))
		cout<<"error: no match!"<<endl;


	vector<double> bel(mypoly::n);
	for(int i=0;i<mypoly::n;++i)
		bel[i]= i*i;

	cout<<"eval ee"<<ee.eval(bel)<<endl;
	cout<<"eval e"<<e.eval(bel)<<endl;



	integ<<"alternative:\n"
		<<ee;
	integ<<"\n\ndifference:\n"<<e-ee<<"\n";
	*/
	integ.close();

	
}

void testIntegrateTsets()
{
	//testcalculate();
	testoptimize();
}

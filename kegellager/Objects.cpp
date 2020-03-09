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

#define _USE_MATH_DEFINES
#include<math.h>
#include<iostream>
#include "Implicit.h"
#include "d3op.h"
#include "quat.h"
#include "Objects.h"


//template<class T> T max(const T&x,const T& y) {return x>y?x:y;}
//template<class T> T min(const T&x,const T& y) {return x<y?x:y;}

Cylinder::tval Cylinder::operator()(const tval* x) const
{
	tval dzmin= x[2]-zmin;
	tval dzmax= x[2]-zmax;

	tval a= dzmin*ideltaz;
	tval r=  a*r2 - (a-1)*r1 ;
	tval daxdax=x[0]*x[0]+x[1]*x[1];
	//make gradient at surface be 1:
	tval schale= (daxdax/r - r)*val;


	return maxi(maxi(-dzmin,dzmax),schale);

}

void Cylinder::grad(tval*ret,const tval* x) const
{

	tval dzmin= x[2]-zmin;
	tval dzmax= x[2]-zmax;

	tval a= dzmin*ideltaz;
	tval r=  a*r2 - (a-1)*r1 ;
	tval daxdax=x[0]*x[0]+x[1]*x[1];
	
	//make gradient at surface be 1:
	tval schale= (daxdax/r - r)*val;

	//dr/dx2 =   r2-r1;
	//dschale/dr = -daxdax/r^2 -1


	tval m1=-dzmin,m2=dzmax;

	tval v=maxi(maxi(m1,m2),schale);
	d3op1(,ret,=tval::EMPTY())

		if(containsTrue(v==schale))
		{
			ret[0].join( x[0]/r *(2*val));
			ret[1].join( x[1]/r *(2*val));
			ret[2].join( (daxdax/r.sqr()+1)* ( (r1-r2)*val) );
		}
		if(containsTrue(v==m1))
		{
			ret[0].join(0);
			ret[1].join(0);
			ret[2].join( -1 );
		}
		if(containsTrue(v==m2))
		{
			ret[0].join(0);
			ret[1].join(0);
			ret[2].join( 1 );
		}
}


//--------------------------------------------------------------------------------


inline myinterval<int> getSeg(const CylKreis &k,const CylKreis::tval& x, const CylKreis::tval& y) {
		return myinterval<int>( floor( atan2(y,x) * (k.n/(2*M_PI)) + 0.5)) + k.n  ;
}



CylKreis::CylKreis(int nkugel,float rkreis,float zmin,float zmax, float rr)
	:n(nkugel*2),R(rkreis),zmin(zmin),zmax(zmax)
{
	if(rr<0)
		r=  R * sin(2*M_PI/(2*n)); 
	else
		r=rr;
	cx=new float[n];
	cy=new float[n];
	for(int i=0;i<n;++i)
	{
		double a=2*M_PI*i/n;
		cx[i]=(float) (cos(a)*R);
		cy[i]=(float) (sin(a)*R);
	}
	
}



CylKreis::tval CylKreis::operator()(const tval* x) const{

	
	tval dzmin= x[2]-zmin;
	tval dzmax= x[2]-zmax;


	tval a= dzmin/(zmax-zmin);
	tval daxdax=x[0]*x[0]+x[1]*x[1];
	//make gradient at surface be 1:
	tval ret= (daxdax - R*R) * ( .5 /R);


	myinterval<int> seg(0,n);
	seg=getSeg(*this,x[0],x[1]);

	for(int ii=seg.a;ii<=seg.b;++ii)
	{
		int i=ii%n;

		tval dx=x[0]-cx[i],dy=x[1]-cy[i];
		tval cyli=(dx*dx+dy*dy-r*r)*.5/r;

		if(ret.a>ret.b)
			ret=cyli;
		else
		{
			if((i&1)==0)
				ret=mini(ret,cyli);
			else
				ret=maxi(ret,-cyli);
		}
	}

	//cut away left and right sides
	return maxi(maxi(-dzmin,dzmax),ret);
}

void CylKreis::grad(tval*ret,const tval* x) const
{

}


//--------------------------------------------------------------------------------



inline myinterval<int> getSeg(const Kugelkreis &k,const Kugelkreis::tval& x, const Kugelkreis::tval& y) {
		return myinterval<int>(floor( atan2(y,x)  * k.mula  -k.rotoff)) + k.nkug  ;
}



Kugelkreis::Kugelkreis(int nkugel,float rkugel,float rkreis,float off, float cz)
	:cz(cz),rotoff(off),rkug(rkugel),rkre(rkreis),nkug(nkugel)
{
	cx=new float[nkug];
	cy=new float[nkug];
	for(int i=0;i<nkug;++i)
	{
		double a=2*M_PI*(i+.5+rotoff)/nkug;
		cx[i]=(float) (cos(a)*rkre);
		cy[i]=(float) (sin(a)*rkre);
	}
	mula=(float)(nkug/(2*M_PI));
}



Kugelkreis::tval Kugelkreis::operator()(const tval* x) const{

	myinterval<int> seg=getSeg(*this,x[0],x[1]);

	tval dz=x[2]-cz;
	tval  dzdz_rr= dz*dz-rkug*rkug;

	tval ret =tval::EMPTY();
	for(int ii=seg.a;ii<=seg.b;++ii)
	{
		int i=ii%nkug;
		tval dx=x[0]-cx[i],dy=x[1]-cy[i];
		tval kugi=dx*dx+dy*dy+dzdz_rr;

		if(ret.a>ret.b)
			ret=kugi;
		else 
			ret=mini(ret,kugi);
	}
	//make gradient at surface be 1:
	return ret*.5/rkug;
}

void Kugelkreis::grad(tval*ret,const tval* x) const{

	myinterval<int> seg = getSeg(*this,x[0],x[1]);

	ret[0]=ret[1]=tval::EMPTY();
	for(int ii=seg.a;ii<=seg.b;++ii)
	{
		int i=ii%nkug;
		tval dx=x[0]-cx[i],dy=x[1]-cy[i];
		ret[0].join(dx);
		ret[1].join(dy);
	}
	ret[2]=x[2];
	//make gradient at surface be 1:
	d3op1(,ret,/=rkug);
}

//--------------------------------------------------------------------------------



inline myinterval<int> getSeg(const ArrangedInCircle &k,const ArrangedInCircle::tval& x, const ArrangedInCircle::tval& y) {
		return myinterval<int>(floor( atan2(y,x)  * k.mula  -k.rotoff)) + k.n  ;
}



ArrangedInCircle::ArrangedInCircle(const Implicit&inner,int n,float rotoff)
	:inner(inner),n(n),rotoff(rotoff)
{
	//the quaternion rotation parameters
	cosa2=new float[n];
	sina2=new float[n];
	
	for(int i=0;i<n;++i)
	{
		//only half the angle:
		double a=  M_PI*(i+.5+rotoff)/n;
		cosa2[i]=(float) (cos(a));
		sina2[i]=(float) (sin(a));
	}
	mula=(float)(n/(2*M_PI));
}



ArrangedInCircle::tval ArrangedInCircle::operator()(const tval* x) const{

	myinterval<int> seg=getSeg(*this,x[0],x[1]);

	//seg.a=0;seg.b=n-1;

	float buf[4]={0,0,0,0};

	tval rotx[3];

	tval ret =tval::EMPTY();
	for(int i=seg.a;i<=seg.b;++i)
	{
		int ii= i%n;
		//rotation around z-axis with angle a:
		buf[2]=-sina2[ii];buf[3]=cosa2[ii];

		d3op2(,rotx,=,x,);
		rotateUsingQuat(rotx,buf);



		tval v=inner(rotx);

		if(ret.a>ret.b)
			ret=v;
		else 
			ret=mini(ret,v);


	}

	return ret;

}

void ArrangedInCircle::grad(tval*ret,const tval* x) const{

	myinterval<int> seg=getSeg(*this,x[0],x[1]);

	seg.a=0;seg.b=n-1;

	float buf[4]={0,0,0,0};

	tval rotx[3];

	tval bufret[3];

	ret[0]=ret[1]=ret[2]=tval::EMPTY();

	for(int i=seg.a;i<=seg.b;++i)
	{
		int ii= i%n;
		//rotation around - z-axis with angle a:
		buf[2]=-sina2[ii];buf[3]=cosa2[ii];

		d3op2(,rotx,=,x,);
		rotateUsingQuat(rotx,buf);
		inner.grad(bufret,rotx);
		buf[2]=-buf[2];
		rotateUsingQuat(bufret,buf);

		ret[0].join(buf[0]);
		ret[1].join(buf[1]);
		ret[2].join(buf[2]);

	}
}
//--------------------------------------------------------------------------------

Rotated::Rotated(const Implicit&inner, float qx,float qy,float qz,float qw)
	:inner(inner)
{
	quat[0]=qx;
	quat[1]=qy;
	quat[2]=qz;
	quat[3]=qw;
	normQuat(quat);
	d3op2(,qbar,=-,quat,);
	qbar[3]=quat[3];
}
		
Rotated::tval Rotated::operator()(const tval* x) const
{
	tval buf[3]={x[0],x[1],x[2]};
	rotateUsingQuat(buf,qbar);
	return inner(buf);
}


void Rotated::grad(tval*ret,const tval* x) const 
{
	tval buf[3]={x[0],x[1],x[2]};
	rotateUsingQuat(buf,qbar);
	inner.grad(ret,buf);
	rotateUsingQuat(buf,quat);
}





Kegel::tval Kegel::operator()(const tval* x) const
{
	tval dx=(x[2]-cz);
	
	tval restrict = ratio>0 ? -dx :dx;

	tval d=dx*ratio;

	return   maxi( ((x[0].sqr()+x[1].sqr())/d-d)*val,restrict );

}

void Kegel::grad(tval*ret,const tval* x) const
{
	tval dx=(x[2]-cz);
	

	tval restrict = ratio>0 ? -dx :dx;

	ret[0]=ret[1]=ret[2]=tval::EMPTY();

	// (dax/d-d)*val ->d/dz->    (-dax/d^2  -1) *val * ratio   
	if(restrict.a<0)
	{
		tval d=dx*ratio;	
		tval dax=(x[0].sqr()+x[1].sqr());
		ret[0]=x[0]/d*2*val;ret[1]=x[1]/d*2*val;ret[2]=(dax/d.sqr()+1)*(-ratio*val);
	}
	//could be outside
	if(restrict.b>0)
	{
		ret[0].join(0);
		ret[1].join(0);
		ret[2].join( ratio>0? -1:1);
	}
}

Kugel::tval Kugel::operator()(const tval* x) const
{
	return (d3prod(x,x)/r-r)*.5;
}

void Kugel::grad(tval*ret,const tval* x) const
{
	d3op2(,ret,=,x,/r);
}

Spiral::tval Spiral::operator()(const tval* x) const
{
	//TODO: make gradient at surface 1 everywhere.

	tval r=sqrt(x[0].sqr()+x[1].sqr());
	tval a=(atan2(x[1],x[0])+phase);
	tval rda = ( r/d + a*(n/(M_PI*2)) );
	tval lower=floor(rda+0.5);
	tval rest = (rda-lower)*d;
	
	float dn2=d*n/(M_PI*2);
	float dn2dn2=dn2*dn2;

	//have equal distances
	tval ms2 = sqrt( tval(dn2dn2) / r.sqr() +1 ) *(ms*.5);

	return (abs(rest)-ms2);

}

void Spiral::grad(tval*ret,const tval* x) const
{
	tval rr=x[0].sqr()+x[1].sqr();
	tval r=sqrt(rr);
	tval a=(atan2(x[1],x[0])+phase);
	tval rda = r/d + a/(tval::PI()*2);
	tval lower=floor(rda+.5);
	tval rest = (rda-lower)*d;
	
	// da/dx[0] = -x[1]/rr, da/dx[1]= x[0]/rr
	// drda/da= 1/(2*pi)
	//drda/d r =1/d
	// d rest/ d rda = d  
	// drest/dr = 1 , drest/da = n*d/(2*pi)
	//drr/dx[0] = 2x[0] , dr/drr = 0.5/sqrt(rr)=0.5/r
	// dr/dx[0] = x[0]/r

	float drestda=n*d/(2*M_PI);
	ret[0]  =  (x[0]*r +   -x[1] *drestda) / rr;
	ret[1]  =  (x[1]*r +    x[0] *drestda) / rr;
	ret[2]=0;

	if( rest.b<=0)
	{
		ret[0]=-ret[0];ret[1]=-ret[1];
	}
	else if(rest.a<0)
	{
		ret[0].join(-ret[0]);ret[1].join(-ret[1]);
	}
}



Offset::tval Offset::operator()(const tval* x) const
{

	tval ov=outer(x);
	return maxi(ov,-(ov+ms) );
}

void Offset::grad(tval*ret,const tval* x) const
{
	tval ov=outer(x);
	outer.grad(ret,x);
		
	tval ov2=-(ov+ms); 


	myinterval<bool> lt= ov<ov2;

	//sure that ov2 greater:
	if( lt.a )
	{
		d3op2(,ret,=-,ret,);	
	}
	else if(lt.b) //ov2 could be greater:
	{
		ret[0].join(-ret[0]);
		ret[1].join(-ret[1]);
		ret[2].join(-ret[2]);
	}
}

Shrink::tval Shrink::operator()(const tval* x) const
{

	tval ov=outer(x);
	return ov+ms;
}

void Shrink::grad(tval*ret,const tval* x) const
{
	outer.grad(ret,x);
}


	
	MirrorZ::tval MirrorZ::operator()(const tval* x) const
	{
		tval xx[]={x[0],x[1],-x[2]};
		return a(xx);
	}

	void MirrorZ::grad(tval*ret,const tval* x) const
	{
		tval xx[]={x[0],x[1],-x[2]};
		a.grad(ret,xx);
		ret[2]=-ret[2];
	}



Difference::tval Difference::operator()(const tval* x) const
{
	tval av=a(x),bv=-b(x);
	return maxi(av,bv);
}

void Difference::grad(tval*ret,const tval* x) const
{
	tval av=a(x),bv=-b(x);

	myinterval<bool> lt = av < bv; 

	//surely smaller
	if(lt.a)
	{
		b.grad(ret,x);			
		d3op2(,ret,=-,ret,);
	}
	//surely greater equals
	else if(!lt.b)
	{
		a.grad(ret,x);
	}
	else //both possible
	{
		a.grad(ret,x);
		tval r[3];
		b.grad(r,x);
		d3op2(,r,=-,r,);
		ret[0].join(r[0]);
		ret[1].join(r[1]);
		ret[2].join(r[2]);
	}

}





Intersection::tval Intersection::operator()(const tval* x) const
{
	tval av=a(x),bv=b(x);
	return maxi(av,bv);
}

void Intersection::grad(tval*ret,const tval* x) const
{
	tval av=a(x),bv=b(x);

	myinterval<bool> lt = av < bv; 

	//surely smaller
	if(lt.a)
	{
		b.grad(ret,x);			
	}
	//surely greater equals
	else if(!lt.b)
	{
		a.grad(ret,x);
	}
	else //both possible
	{
		a.grad(ret,x);
		tval r[3];
		b.grad(r,x);
		ret[0].join(r[0]);
		ret[1].join(r[1]);
		ret[2].join(r[2]);
	}

}


Union::tval Union::operator()(const tval* x) const
{
	tval av=a(x),bv=b(x);
	return mini(av,bv);
}

void Union::grad(tval*ret,const tval* x) const
{
	tval av=a(x),bv=b(x);

	myinterval<bool> lt = av < bv; 

	//surely smaller
	if(lt.a)
	{
		a.grad(ret,x);			
	}
	//surely greater equals
	else if(!lt.b)
	{
		b.grad(ret,x);
	}
	else //both possible
	{
		a.grad(ret,x);
		tval r[3];
		b.grad(r,x);
		ret[0].join(r[0]);
		ret[1].join(r[1]);
		ret[2].join(r[2]);
	}

}







Scaled::tval Scaled::operator()(const tval* x) const
{
	tval sc[]={x[0]*iscale,x[1]*iscale,x[2]*iscale};
	return a(sc)*scale;
}

void Scaled::grad(tval*ret,const tval* x) const
{
	tval sc[]={x[0]*iscale,x[1]*iscale,x[2]*iscale};
	a.grad(ret,sc);
}

Translated::tval Translated::operator()(const tval* x) const
{
	tval sc[]={x[0]-t[0],x[1]-t[1],x[2]-t[2]};
	return a(sc);
}

void Translated::grad(tval*ret,const tval* x) const
{
	tval sc[]={x[0]-t[0],x[1]-t[1],x[2]-t[2]};
	a.grad(ret,sc);
}





LTZ::tval LTZ::operator()(const tval* x) const
{
	return x[2]-pz;
}

void LTZ::grad(tval*ret,const tval* x) const
{
	ret[0]=0;ret[1]=0;ret[2]=1;
}





Plane::tval Plane::operator()(const tval* x) const
{
	return d3prod(x,n)-d;
}

void Plane::grad(tval*ret,const tval* x) const
{
	d3op2(,ret,=,n,);
}

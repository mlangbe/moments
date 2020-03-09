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

#ifndef OBJECTS_H
#define OBJECTS_H

#pragma once

#define _USE_MATH_DEFINES
#include<math.h>

#include "Implicit.h"
#include "d3op.h"


//template<class T> T max(const T&x,const T& y) {return x>y?x:y;}
//template<class T> T min(const T&x,const T& y) {return x<y?x:y;}


struct Cylinder :
	 Implicit
{

	float r1,r2;
	float zmin;
	float zmax;
	float ideltaz;
	float val;

/**
 * 
 * @param r radius
 * @param zmin z-koordinate unteres ende
 * @param zmax z-koordinate oberes ende
 */
	 Cylinder( float r1,float r2,float zmin,float zmax)
		 :r1(r1),r2(r2),zmin(zmin),zmax(zmax),ideltaz(1.0/(zmax-zmin))
	{
		float tana = (r2-r1)*ideltaz;
		val=.5/sqrt(1+tana*tana);
	}

	tval operator()(const tval* x) const;
	void grad(tval*ret,const tval* x) const;
};



 struct Kugelkreis : Implicit {
	
	float rkug;
	float rkre;
	float rotoff; //rotation in units of 2PI/nkug;
	
	int nkug;
	
	//centers
	float*cx; 
	float*cy;
	float cz;
	float mula;
	
	Kugelkreis(int nkugel,float rkugel,float rkreis,float off=0,float cz=0);
	
		
	tval operator()(const tval* x) const;

	void grad(tval*ret,const tval* x) const ;

};




 /**
 rotiertes.
 */
struct Rotated : Implicit {

	const Implicit&inner;

	float quat[4],qbar[4];

	Rotated(const Implicit&inner, float qx,float qy,float qz,float qw);
		
	tval operator()(const tval* x) const;

	void grad(tval*ret,const tval* x) const ;

};


/** arranges n rotated instances of inner. it is assumed that inner extends most to the x-axis and is symmetric in y*/
struct ArrangedInCircle:Implicit
{

	
	const Implicit&inner;

	int n;

	//centers
	float*cosa2; 
	float*sina2;

	float mula;

	float rotoff;

	ArrangedInCircle(const Implicit&inner, int n,float rotoff=0);

	tval operator()(const tval* x) const;

	void grad(tval*ret,const tval* x) const ;

	~ArrangedInCircle(){
		delete[]cosa2;delete[]sina2;
	}
};











struct CylKreis : Implicit {
	
	/** radius of cylinder the small cylinders are attached to*/
	float R;
	
	int n;
	
	//centers
	float*cx; 
	float*cy;

	float zmin,zmax;
	float r;
	
	CylKreis(int n,float R,float minz,float maxz,float r=-1);
	
	tval operator()(const tval* x) const;

	void grad(tval*ret,const tval* x) const ;

	~CylKreis(){
		delete[]cx;delete[]cy;
	}
};




 struct Kegel : Implicit {
	
	/**
	 * 
	 * @param cz z-koordinate kegelspitze
	 * @param tana tangens von oeffnungswinkel: >0 =kegel oeffnet nach +z,<0:kegel oeffnet nach -z
	 */
	Kegel(float cz_,float tana)
	{
		ratio=tana;cz=cz_;
		val=.5/sqrt(1+ratio*ratio);
	}

	float ratio,cz;

	//scale ta make gradient 1
	float val;
		
	
	tval operator()(const tval* x) const;
	void grad(tval*ret,const tval* x) const;
};


 struct Kugel : Implicit {
	
	/**
	 * 
	 * @param cz z-koordinate kegelspitze
	 * @param tana tangens von oeffnungswinkel: >0 =kegel oeffnet nach +z,<0:kegel oeffnet nach -z
	 */
	Kugel(float rr)
		:r(rr)
	{
	}

	float r;
		
	
	tval operator()(const tval* x) const;
	void grad(tval*ret,const tval* x) const;
};





 struct Spiral : Implicit
 {
	/** .5*thickness*/
	float ms;

	/** distance in radius direction between two layers*/
	float d;

	/** phase shift*/
	float phase;

	/** number of arms. other rotation dir if negative*/
	int n;

	Spiral(float d_, float ms_,int n_=1,float phase_=0)
		:d(d_),ms(ms_),phase(phase_),n(n_)
	{}

 	
	tval operator()(const tval* x) const;
	void grad(tval*ret,const tval* x) const;
 
 
 
 };






 struct Offset : Implicit {

	const Implicit& outer;
	const float ms;
	
	Offset(const Implicit& o,float mss)
		:outer(o),ms(mss)
	{
	}

	tval operator()(const tval* x) const;

	void grad(tval*ret,const tval* x) const;	
};

 struct Shrink : Implicit {

	const Implicit& outer;
	const float ms;
	
	Shrink(const Implicit& o,float mss)
		:outer(o),ms(mss)
	{
	}

	tval operator()(const tval* x) const;

	void grad(tval*ret,const tval* x) const;	
};

 struct Difference : Implicit {

	
	Difference(const Implicit& aa,const Implicit& bb)
		:a(aa),b(bb)
	{
	}
	
	tval operator()(const tval* x) const;

	void grad(tval*ret,const tval* x) const;


	const Implicit &a;
	const Implicit &b;

};



 struct Intersection : Implicit {

	
	 Intersection(const Implicit& aa,const Implicit& bb)
		:a(aa),b(bb)
	{
	}
	tval operator()(const tval* x) const;

	void grad(tval*ret,const tval* x) const;
	

	const Implicit &a;
	const Implicit &b;
};


struct Union : Implicit {

	
	 Union(const Implicit& aa,const Implicit& bb)
		:a(aa),b(bb)
	{
	}
		
	tval operator()(const tval* x) const;

	void grad(tval*ret,const tval* x) const;

	const Implicit &a;
	const Implicit &b;
};


 struct Scaled : Implicit {

	const Implicit& a;
	float iscale,scale;
	 
	Scaled(const Implicit& aa,float scale)
		 :a(aa),scale(scale),iscale(1/scale)	
	{
	}
	
	
	
	tval operator()(const tval* x) const;
	void grad(tval*ret,const tval* x) const;
	
};

struct MirrorZ : Implicit {

	const Implicit& a;
	 
	MirrorZ(const Implicit& aa)
		:a(aa)
	{
	}
	
	
	
	tval operator()(const tval* x) const;
	void grad(tval*ret,const tval* x) const;
	
};

 
 struct Translated : Implicit {

	const Implicit &a;
	float t[3];
	
	Translated(const Implicit&aa,float ttx,float tty,float ttz)
		:a(aa)
	{
		t[0]=ttx;t[1]=tty;t[2]=ttz;
	}

	tval operator()(const tval* x) const;

	void grad(tval*ret,const tval* x) const;

};


 struct LTZ : Implicit {

	
	const float pz;
	
	LTZ(float z)
		:pz(z)
	{
	}
	
	tval operator()(const tval* x) const;

	void grad(tval*ret,const tval* x) const;

};


 struct Plane : Implicit {

	
	float n[3];float d;
	
	Plane(float ttx,float tty,float ttz, float dist)
	{

		n[0]=ttx;n[1]=tty;n[2]=ttz;
		float l=sqrt(d3prod(n,n));
		d3op1(,n,/=l);
		d=dist/l;
	}
	tval operator()(const tval* x) const;

	void grad(tval*ret,const tval* x) const;
	

};



#endif
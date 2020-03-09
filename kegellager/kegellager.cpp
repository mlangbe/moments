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


#define _USE_MATH_DEFINES
#include<cmath>
#include<iostream>
#include<myinterval.h>
#include "Objects.h"
#include "extractsurface.h"
#include "composites.h"


//using namespace std;


/** laufverhaltnis aussenkegel /mittenkegel*/
double C=3;
/** laufverhaeltnis mittenkegel/achsenkegel*/
double A=3./2;

double sinb=.5*sqrt( (4-pow(C-A,2))/(C*A+1) );
double sina=A*sinb;
double sing=C*sinb;
double beta = asin(sinb);
double alpha= asin(sina);
double _gamma_=alpha+2*beta;
double delta=alpha+beta;
double cosa=cos(alpha);
double cosb=cos(beta);
double cosg=cos(_gamma_);
float tana=(float)( sina/cosa );
float tanb=(float)(sinb/cosb);
float tang=sing/cosg;
float tand=(float)(tan(delta));
double r1=tana;
double l2=1/cosa;
double r2=sinb*l2;
double r3=l2*sing;
float l2cosb=(float)(l2*cosb);
double r5=r3/2;


int ninner=8;
int nax=ninner*A;
int nouter=ninner*C;

//radius fuer kugeln des zahnkranzes
double rkranz=r2*0.75;
//radius loecher
double rmk= rkranz *sin( 2*M_PI /( 2* 2*ninner) ) ;
//radius Kugeln
double rmki=rmk;


//Oeffn. winkel kegel Zahnkranz:
double zka= asin(rmk/(l2*.75));
//kegel +aussparung
double a=zka*0.95,ia=zka*1.05;

double l2s=l2*5/8,l2e=l2*7/8,l2sI=l2s -l2 *.1/8,l2eI=l2e +l2*.1/8;



float scale=20;
float ms=1;

float ms0=ms/scale;



Cylinder stift(ms0,ms0,-2*ms0,2*ms0);
Cylinder hole_0(ms0,ms0,-3,3);
Kegel cap(3*ms0,-1);
Intersection hole(hole_0,cap);




//#define BALLSINVERTED

#define SPIRAL_SUPPORT

#define KEGELZAHNKRANZ




Cylinder cyl(r1/2,r1/2,-1,1);
Cylinder k1(r1/2,r1*3/2,-.5,.5),k2(r1*3/2,r1/2,-.5,.5);
Intersection dkg(k1,k2);
Union ak(dkg,cyl);

#ifndef KEGELZAHNKRANZ
Kugelkreis lak1(nax, rmki, r1*.75,0,-.25),lak1I(nax,rmk,r1*.75,0.5,-.25), lak2(nax, rmki, r1*.75,0,.25),lak2I(nax,rmk,r1*.75,0.5,.25);
#else
KegelKegel lak1(alpha,a,nax,l2s,l2e,-1,0,tana),lak1I(alpha,ia,nax,l2sI,l2eI,-1,.5,tana);
MirrorZ lak2(lak1),lak2I(lak1I);
#endif

#ifndef BALLSINVERTED
Union 
#else
Difference
#endif
	uak1(ak,lak1),uak2(uak1,lak2);

Difference dak1(uak2,lak1I);

Cylinder cylI(4*ms0, 4*ms0, -3, 3);


Difference ak0(dak1,cylI),achsenKegel(ak0,lak2I);


LTZ lt0(0);


Translated holeT(hole,r1*.75-ms0/2,0,0);
ArrangedInCircle holeTA(holeT,12);



Difference achsenKegelUP0(achsenKegel,lt0);
//lower part
Difference achsenKegelUP(achsenKegelUP0,holeTA);




Kegel m1(-r2/tanb,tanb),m2(r2*tand,-1/tand);
LTZ ltr2t(-r2/tanb*.5),lt2((r2 - 2*ms0  )  *tand );
Intersection doppelKegel(m1,m2);
Intersection dk2(doppelKegel,lt2);
Difference mk(dk2,ltr2t);

#ifndef KEGELZAHNKRANZ
Kugelkreis lmk1(ninner, rmki, r2*.75, 0,-r2/tanb*.25 ),lmk1I(ninner,rmk, r2*.75, 0.5,-r2/tanb*.25);
#else
KegelKegel lmk1(beta,a,ninner,l2s,l2e,-l2*cosb,0,tana),lmk1I(beta,ia,ninner,l2sI,l2eI,-l2*cosb,.5,tana);
#endif


#ifndef BALLSINVERTED
Union 
#else	
Difference
#endif
	umk1(mk,lmk1);
Difference mittenKegel(umk1,lmk1I);




MirrorZ mittenKegelLP00(mittenKegel);
Difference mittenKegelLP0(mittenKegelLP00,lt0);
//lower part
Difference mittenKegelLP(mittenKegelLP0,hole_0);

Difference mittenKegelUP0(mittenKegel,lt0);
//upper part
Difference mittenKegelUP(mittenKegelUP0,hole_0);






double h= ms0;

double dck=1-cosg*l2-h/tang,rik=sing*l2+h,rck=rik+h,rok=rck+h;

//min. dist of parts fitted together
double dist=.2/20;



double dkk= 1 -  l2*cosg *0.75;
double ekk = 1 -l2*cosg*.5;
double rkk = l2*sing*0.75;



double dkecut= 1 - l2*cosg*0.5  - .5*l2*sing*tang;
double dkecut2= 1 -l2*cosg -l2*sing*tang;

#ifndef KEGELZAHNKRANZ
Kugelkreis kko1(nouter,rmki,rkk,0,-dkk),kko1I(nouter,rmki,rkk,0.5,-dkk)  ;  
#else
KegelKegel kko1(_gamma_,a,nouter,l2s,l2e,-1,0,tana),kko1I(_gamma_,ia,nouter,l2sI,l2eI,-1,.5,tana);
#endif

Kegel ke1I_0(-1,tang);


#ifndef BALLSINVERTED
Difference
#else
Union
#endif	
	keiv3(ke1I_0,kko1);

Union ke1I(keiv3,kko1I);

Translated ke1(ke1I,0,0,-2*h/tang);



//Kegel ke1I2(-dkecut,-1/tang);
Cylinder ke1I2(r5,r5,-1,0);
Cylinder outer(rok,rok,-1+l2*cosg*.5 -ms0/2,-1+l2*cosg+h/tang);

Spiral sp1( rok/6, ms0/2, 24,M_PI/24);
Spiral sp2( rok/6, ms0/2, -24,M_PI/24);
Union sp12(sp1,sp2);


Union keee0(ke1,sp12);
Cylinder basePlate( rok,rok,-1+l2*cosg*.5-ms0/2,-1+l2*cosg*.5);
Union keee(keee0,basePlate);
Intersection keiv2_00(keee,outer);
Difference keiv2_0(keiv2_00,ke1I);


Difference keiv2(keiv2_0,ke1I2);



Cylinder cok(rok,rok,-dck,dck-dist);

Cylinder cck(rck,rck,-dck,dck-dist);

Cylinder cik(rik,rik,-dck,dck-dist);

Kegel capik(dck-dist+rik+h*.5,-1);

Cylinder knopf000(2*ms0,ms0/2,rck-ms0/2,rck+ms0/2);
Rotated knopf00(knopf000,0,1,0,1);

ArrangedInCircle knopf0(knopf00,2,-0.05),knopf1(knopf00,2,-0.45),knopfI0(knopf00,2,0.05),knopfI1(knopf00,2,0.45);
Union knopf(knopf0,knopf1),knopfI(knopfI0,knopfI1);

float dpla=h;
Plane pl1(-1,0,-1,-dck+dist+dpla),pl2(1,0,-1,-dck+dist+dpla);
Intersection pll(pl1,pl2);

Plane pl3(0,-1,-1,-dck+dist+dpla),pl4(0,1,-1,-dck+dist+dpla);
Intersection plll(pl3,pl4);

Union pla(pll,plll);


Plane pld1(-1,-1,0,h/2),pld2(1,-1,0,h/2),pld3(1,1,0,h/2),pld4(-1,1,0,h/2);
Intersection raute000(pld1,pld2),raute00(raute000,pld3),raute0(raute00,pld4);
Translated raute(raute0,rck,0,0);
ArrangedInCircle rauten(raute,4,0.5);





Plane cutX(1,0,0,0),cutY(0,1,0,0);
Intersection cutQ1(cutX,cutY);

Union cutQ2(cutX,cutY);

Difference cutQ(cutQ2,cutQ1);

Shrink cutQo(cutQ,dck);
Shrink cutQi(cutQ,-dck);


Difference ocyl00(cok,cck);
Difference ocyl01(ocyl00,knopfI);
Intersection ocyl0(ocyl01,cutQ);
Difference ocyl(ocyl0,cutQo);

Difference icyl00(cck,cik);
Difference icyl01(icyl00,cutQ);
Intersection icyl02(icyl01,capik);
Intersection icyl0(icyl02,cutQi);
Union icyl(icyl0,knopf);

Union outerJoin0(icyl,ocyl);

Difference outerJoin1(outerJoin0,pla);
Difference outerJoin(outerJoin1,rauten);







/*
MirrorZ keivX(keiv2);


//Union kkeiv (keivX,keiv4),showIt(kkeiv,achsenKegel);



Cylinder Ok( rok,rok,-dck+dist,dck);

CylKreis ck ( 2, rck,-dck,dck-dist,h*.5);
CylKreis ckm( 2, rck+dist,-dck,dck,h*.5); 
Cylinder ckI( rik,rik, -dck,dck);
Difference dckck(ck,ckI);


Difference dckck2(Ok,ckm);


Union leftOuterPart ( dckck2,keivX),rightOuterPart(dckck,keiv2);
*/


Union leftOuterPart(keiv2,outerJoin);
















Scaled aks(achsenKegelUP,scale),akm(mittenKegelLP,scale),akmU(mittenKegelUP,scale),akl(leftOuterPart,scale),akst(stift,scale);
Offset offmid(akm,ms);
Offset offmidU(akmU,ms);








int main(int argc, char* argv[])
{

	/*
	myinterval<float> x;

	for(float w=M_PI/16;w<M_PI*2;w+=M_PI/8)
	{

		x=w;
		cout<<w/M_PI<<" "<<atan2(sin(x),cos(x))/M_PI<<endl; 
	}
	*/
	myinterval<float> box3[]={myinterval<float>(-4,4),myinterval<float>(-4,4),myinterval<float>(0,dck)};

	
	extractsurface().toStl(outerJoin,"test.stl",1.*ms0/4,box3);
	
	myinterval<float> box[]={myinterval<float>(0,32),myinterval<float>(0,32),myinterval<float>(-32,32),};
	
	
	myinterval<float> box2[]={myinterval<float>(-32,32),myinterval<float>(-32,32),myinterval<float>(-32,32)};	
	//extractsurface().toStl(aks,"achsenkegel_4.stl",0.5,box);
	//extractsurface().toStl(offmid,"mittenkegel_4.stl",0.5,box);
	//extractsurface().toStl(akmU,"mittenkegelU_4.stl",0.5,box);
	extractsurface().toStl(akl,"leftOuterPart_4.stl",0.5,box);
	//extractsurface().toStl(akst,"stift_4.stl",0.5,box);

	//extractsurface().toStl(aks,"achsenkegel.stl",0.25,box2);
	//extractsurface().toStl(offmid,"mittenkegel.stl",0.25);
	//extractsurface().toStl(akmU,"mittenkegelU.stl",0.25,box2);
	extractsurface().toStl(akl,"leftOuterPart.stl",.25,box2);
	//extractsurface().toStl(akst,"stift.stl",0.25);
	
	
	
	
	

	return 0;
}


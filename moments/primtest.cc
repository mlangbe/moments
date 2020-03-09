/*
	test some prime-number operations

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

#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<pthread.h>
#include<unistd.h>

#include "primtest.hh"
#include "infnum.hh"

mythreadvars::mythreadvars(){
  next=0;
  pthread_mutex_init(&mymutex,NULL);
  pthread_cond_init(&nextvalid,NULL);        
  pthread_cond_init(&nextempty,NULL);        
}

mythreadvars::~mythreadvars(){
  if(next)delete next;
  pthread_mutex_destroy(&mymutex);
  pthread_cond_destroy(&nextvalid);        
  pthread_cond_destroy(&nextempty);        
}


inline void mythreadvars::testord()
{
  infnum myn = primz;
  --myn;

  for(int i=0;notready;++i){

    int myfakt;

    //get next prime faktor id to check
    pthread_mutex_lock(&mymutex);

    if( nextfakt < nfakts ){
      myfakt=nextfakt;++nextfakt;
    }
    else{
      pthread_mutex_unlock(&mymutex);
      break;
    }
	
    pthread_mutex_unlock(&mymutex);
      
    cout<<"testing prime factor"<<primfakts[myfakt]<<endl;

    if( powmod
	(a,myn/((infnum::t1elem)primfakts[myfakt]),primz)
	== 1 ){

      pthread_mutex_lock(&mymutex);		
      notready=false;
      pthread_mutex_unlock(&mymutex);		

    }

  }

  cout<< getpid() <<" terminated "<<endl;

}

void* mythreadvars::testord_v(void*x)
{
  ((mythreadvars*)x)->testord();
  return 0;
}



inline void mythreadvars::geninfnum()
{
  infnum * mynum=0;

  static const infnum::t1elem my2=2;
  static const int imax=100000;

  for(int i=0; !mynum && i<imax; ++i){

    mynum = new infnum;      
    genfactorizable(*mynum,log(2)*1024);
    *mynum = *mynum * my2 + infnum(1);

    pthread_mutex_lock(&mymutex);
    {
      //	cout<<"number "<<i<<" ("<<i*100./imax<<"%)"<<primz<< endl;
      if(notready){
	  
	if(next)
	  pthread_cond_wait(&nextempty,&mymutex);

	if(notready){
	  if(next!=0)
	    cout<<"uiuiui"<<endl;
	  next=mynum;mynum=0;
	  pthread_cond_signal(&nextvalid);	
	}
      }       
    }
    pthread_mutex_unlock(&mymutex);
  }

  if(mynum)
    delete mynum;

  else{//i>=imax
    pthread_mutex_lock(&mymutex);      
    {
      primz=0;
      notready=false;
      pthread_cond_broadcast(&nextvalid);	  	
    }
    pthread_mutex_unlock(&mymutex);
  }
}


inline void mythreadvars::testprime()
{
  infnum*mynum=0;

  for(int i=0;!mynum;++i){

    pthread_mutex_lock(&mymutex);
    {
      if(notready){
	if(!next)
	  pthread_cond_wait(&nextvalid,&mymutex);
	//get next number to check

	if(notready){
	  if(!next)
	    cout<<"strange"<<getpid()<<' '<<next<<' '<<primz<<' '<<i<<endl;
	}
	mynum=next;next=0;
	pthread_cond_signal(&nextempty);
      }
    }
    pthread_mutex_unlock(&mymutex);
      
    if(!mynum)
      break;

    if(primzahlen.test(*mynum)!=0 && pt(*mynum) ){
      pthread_mutex_lock(&mymutex);
      {
	if(notready){
	  notready=false;
	  cout<<"yuhuu!"<<*mynum<<endl;
	  primz=*mynum;
	  if(next){delete next;next=0;}
	  //wake up all other threads so they can terminate
	  pthread_cond_broadcast(&nextvalid);	  
	  pthread_cond_broadcast(&nextempty);	  
	}
	else{
	  cout<<"found twice"<<endl;
	  delete mynum;mynum=0;	    
	}
      }
      pthread_mutex_unlock(&mymutex);		
      break;
    }
    else{
      cout<<"no"<<flush;
      delete mynum; mynum=0;
    }

  }

  cout<< getpid() <<" terminated "<<endl;

}

void* mythreadvars::testprime_v(void*x)
{
  ((mythreadvars*)x)->testprime();
  return 0;
}

  

//initiates multiple threads to search for randomly selected primes
infnum& mythreadvars::getprime()
{
  static const int nthreads=4;
  pthread_t testp[nthreads],geni;
  //    pthread_create(&geni,NULL,&geninfnum_v,this);

  do{

    //.............................
    //find 1 possible prime number

    notready=true;
      
    for(int i=0;i<nthreads;++i){
      pthread_create(&(testp[i]),NULL,&testprime_v,this);
      cout<<"created thread"<<testp[i]<<endl;      
    }      
      
    //generate number to be tested
    geninfnum();

    //and test them
    for(int i=0;i<nthreads;++i)
      pthread_join(testp[i],NULL);

    //if too many iterations:
    if(primz==0)
      break;
    //.............................

    //test if it really is a prime number

    infnum n1= primz;--n1;

    nfakts=primfaktors(n1,primfakts,primzahlen.n,10000);

    if( n1!=1 && prim(n1) !=1 )
      continue;        

    //test log(n1)*10 a's if they have order primz-1
    for(int j=0;j<log(primz)*10;j++){

      a = makeEqualDist(primz-infnum(2))+infnum(1);
      notready=true;
      nextfakt=0;

      if(powmod(a,primz-infnum(1),primz)!=1){
	break;
      }

      for(int i=0;i<nthreads;++i){
	pthread_create(&(testp[i]),NULL,&testord_v,this);
	cout<<"created thread"<<testp[i]<<endl;      
      }      
	
      for(int i=0;i<nthreads;++i)
	pthread_join(testp[i],NULL);

      //if ord_primz(a)=primz-1
      if(notready){
	notready=false;
	break;
      }
    }
  }while( notready );

  return primz;
}



int genprime(infnum&n)
{
  infnum::stnkomm=0;


  static const infnum::t1elem  my2=2;
  static const infnum my1=1;
  int nmax=10000;
  char buf[128];
  for(int i=0;i<nmax;i++){
    genfactorizable(n,log(2)*1024);
    n=n*my2+my1;
    
    snprintf(buf,sizeof(buf),"\r%5d (%5.1f%%) bits n: %d",
	     i,i*100.0/nmax,
	     n.Nsh()*infnum::bits1elem);
    cerr<<buf<<flush;
    
    int j=prim(n);

    switch(j)
      {
      case 1:cout<<n<<" is prime"<<endl;return 0;
      case 0:break;
      default:case -1:cout<<n<<"could be prime, I don't know"<<endl;break;
      }  

  }
  return 1;
}


void mydiv(infnum&numer,infnum::t1elem &denom)
{
  denom=numer.t1elemdiv(denom);
}

void mydiv(infnum&numer,infnum&denom)
{
  infnum quotient;
  numer.quot(denom,quotient);
  denom.swap(numer);
  numer.swap(quotient);  
}

void mydiv(long int &numer,long int&denom)
{
  ldiv_t ret=ldiv(numer,denom);
  numer=ret.quot;
  denom=ret.rem;
}

gencandidate::gencandidate(unsigned long level)
{
  if(level==1){
    n=2;phin=1;
    a=new ulint[1];
    a[0]=1;
  }
  else{
    gencandidate gc(level-1);
    ulint prim=gen(gc,1);
    
    n=gc.n*prim;
    phin = gc.phin*(prim-1);
    
    // cout<<"level"<<level<<"prim:"<<prim<<"phin"<<phin<<" n"<<n<<endl;
    a=new ulint[phin];
    ulint *next=a;
    
    for(ulint i=0;i<gc.phin*prim;i++){
      if(gen(gc,i)%prim!=0){
	//	  cout<<gc[i]<<endl;
	*(next++)=gen(gc,i);      
      }
    }
  }
  if(n>infnum::maxt1elem)
    {fprintf(stderr,"level too big for gencandidate and infnum, %u > %u",
	     n,(unsigned)infnum::maxt1elem);}
}

infnum makeEqualDist(const infnum&n)
{
  infnum ret,myn,nn=n;
  do{

    for(infnum::tindex i=0;i<n.N()*infnum::bits1elem/48;++i){

      myn=nn*infnum(drand48());
      nn*=infnum(ldexp(1,-48));
      //round:
      nn.resizebottom(nn.Nsh()-nn.N());
      myn.resizebottom(myn.Nsh()-myn.N());
      
      if(i==0)
	ret=myn;
      else  
	ret+=myn;
      
    }
  }while( ret > n );
    
  return ret;
  
  
    
}

gencandidate gc(5);

primedb<> primzahlen;





/*
	newton algorithm using interval arithmetics

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

#ifndef INTERVALNEWTON_H
#define INTERVALNEWTON_H

#pragma once

#include<vector>
#include<deque>
#include<math.h>
#include<float.h>
#include "myinterval.h"

#if 0
#if defined(GNUCC) && ! defined(isfinite)
#define isfinite(x) _finite(x)
#endif
#endif

/*
*give the solution intervals for the equation f(x)=0.
*\param f: the function
*\param df: the function's derivative, in interval arithmetics
*\param start: the interval where to search for solutions
*\param ret: the solution intervals
*/
template<class T,class tf,class tdf>
void solve(std::vector< myinterval<T> > &ret,
			const myinterval<T> & start,const tf&f,const tdf&df	,
			T eps=1e-12
			)
{
	using namespace std;
	ret.clear();
	cout.precision(16);
	deque<myinterval<T> > todo(1,start);
	myinterval<T> dy,s,r;
	T y,x;
	int iter=0;
	while(!todo.empty())
	{
		s = todo.front(); todo.pop_front();
		x = s.mid();//s.a +s.width()* ( (iter&3) +0.5 )*0.25;

		if(s.width()<=eps*x )
		{	ret.push_back(s); continue;}

		y = f(x);

		dy  = df(s);
		//cout<<"\ns="<<s<<" f(x)="<<y<<" f'(s)="<<dy<<endl;
		//cout<<"s.width()="<<s.width()<<endl;
		++iter;
		//cout<<"NO."<<iter<<endl;


		if( dy.b<dy.a)
		{	ret.push_back(s); continue;}

		if(dy.contains(0))
		{

			if(y > 0)
			{
				r.b = x - y / dy.b;
				r.a = x - y / dy.a;
			}
			else
			{
				r.b = x - y / dy.a;
				r.a = x - y / dy.b;
			}
			
			//cout<<"Cr="<<r<<endl;
			
			if( 
				isfinite(r.a)
				&& r.a <= s.b )
			{
				//cout<<"news:"<<myinterval<T>(r.a,s.b)<<endl;
				if(r.a>s.a)
					todo.push_back(myinterval<T>(r.a,s.b));
				else
					ret.push_back(s);
			}
			if( 
				isfinite(r.b)
				&& r.b >= s.a )
			{
				//cout<<"news:"<<myinterval<T>(s.a,r.b)<<endl;
				if(r.b<s.b)
					todo.push_back(myinterval<T>(s.a,r.b));			
				else
					ret.push_back(s);
			}
		}
		else
		{
			if(y==0)
			{ret.push_back(s);continue;}

			r = -y; r/=dy; r+=x;
			
			//cout<<"r="<<r<<endl;

			r.intersect(s);
			//cout<<"r^s="<<r<<endl;


			if(r.a<=r.b)
			{
				if(r.width()<s.width())		
					todo.push_back(r);
				else
					ret.push_back(r);
			}
		}
	}
}



/*
*give the solution intervals for the equation f(x)=0.
*\param f: the function
*\param df: the function's derivative, in interval arithmetics
*\param start: the interval where to search for solutions
*\param ret: the solution intervals
*/


template<class T,class tf,class tdf>
void solve2(std::vector< myinterval<T> > &ret,
			const myinterval<T> & start,const tf&f,const tdf&df	,
			T eps=1e-12
			)
{
	using namespace std;
	ret.clear();
	//cout.precision(16);
	vector< myinterval<T> > todo;
	
	todo.push_back(start);
	
	myinterval<T> dy,invdy, s,r;
	T y,x;
	int iter=0;
	while(!todo.empty())
	{
		s = todo.back(); todo.pop_back();

		if(s.width()<=eps*s.mid() )
		{	ret.push_back(s); continue;}

		dy  = df(s);
	
		if( dy.b<dy.a)
		{	ret.push_back(s); continue;}
	
	
		if(dy.contains(0))
		{
			invdy.a=1./dy.a;
			invdy.b=1./dy.b;
	
			bool usea = isfinite(invdy.a);
			bool useb = isfinite(invdy.b);
		  
			if(!(usea|useb)){	
				if(f(s.mid()).mid()==0)
					ret.push_back(s); 
				continue;
				}
				

			for(int i=1;;++i)
			{
				x = s.mid();
				y = f(x).mid();

			  //cout<<"\ns="<<s<<" f(x)="<<y<<" f'(s)="<<dy<<endl;
			  //cout<<"s.width()="<<s.width()<<endl;
				++iter;
			  //cout<<"NO."<<iter<<endl;

				r.a=-HUGE_VAL; r.b=HUGE_VAL;
				
				if(usea) (y>0? r.b : r.a) = x - y * invdy.a ;
				if(useb) (y>0? r.a : r.b) = x - y * invdy.b;
				
				
				//it doesn't change any more
				if(s.b==r.a || s.a==r.b)
				{
					ret.push_back(s);break;
				}
				
				
				bool salera = (s.a<=r.a),
						 rblesb = (r.b<=s.b); 
				
				//a split
				if(salera & rblesb )
				{
					todo.push_back(myinterval<T>(s.a,r.a));
					todo.push_back(myinterval<T>(r.b,s.b));
					break;
				}
				
			  //no zero inside s
				if(!(salera|rblesb))
					break;
					
				//empty interval on one side of s
				if(salera)
					s.b=r.a;
				else
					s.a=r.b;
					
				if(i>=4)
				{	
					todo.push_back(s);
					break;
				}
				
			}			
			
		}
		else //if dy does not intersect 0 : max. 1 solution.
		{
			invdy.a=1./dy.a;
			invdy.b=1./dy.b;
							

			for(int i=1;;++i)
			{
				if( ! f(s).contains(0) )
					break;

				x=s.mid();
				y=f(x).mid();
			
				if(y==0)
				{ret.push_back(s);break;}
				
				r.a = x - y * invdy.a ;
				r.b = x - y * invdy.b;
			
				if(r.a>r.b)
				{x=r.a;r.a=r.b;r.b=x;}
				//cout<<"r="<<r<<endl;


				//as we do not do exact rounding:
				//still find a solution on presence of rounding  errors.
				if(r.b <= s.a  && r.b + s.width() * eps >=s.a)
				{
					ret.push_back(s.a); break;
				}
				if(r.a>=s.b  &&r.a - s.width() * eps <=s.a)
				{
					ret.push_back(s.b); break;
				}

				r.intersect(s);
				//cout<<"r^s="<<r<<endl;


				if(r.a>r.b)
					break;
					
			
				if(
				  (r.width()>=s.width())
				  |
				  (r.width()<=eps*r.mid())
				  )
				{ ret.push_back(r); break; }		
				
				
				if(i>=4)
				{
					todo.push_back(r);
					break;
				}
				
				s=r;
			}
		}
	}
}

template<class T,class tf,class tdf>
struct solve3
{
  tf f;
  tdf df;
  std::vector<myinterval<T> > ret;   
  T eps;

	myinterval<T> dy,invdy;
	T y,x;

  solve3(){}
  
	solve3(std::vector< myinterval<T> > &re,
			const myinterval<T> & start	,
			T eps=1e-12
			)
			:eps(eps) 
  {
  	re.swap(ret);
  	
  	ret.clear();
  	rec(start);
  	
  	re.swap(ret);  	
  }
  
  void solve(std::vector<myinterval<T> > &re, const myinterval<T>&start,T leps=1e-12 )
  {
    eps=leps;
  	ret.clear();
  	rec(start);
  	re.swap(ret);	  
  }
  
  
  void rec(myinterval<T> s)
	{
		myinterval<T> r;
	
		int iter=0;
		for(;;)
		{	

			if(s.width()<=eps*s.mid() )
			{	ret.push_back(s); return;}

			dy  = df(s);
	
			if( dy.b<dy.a)
			{	ret.push_back(s); return;}
	
	
			if(dy.contains(0))
			{
				invdy.a=1./dy.a;
				invdy.b=1./dy.b;
				bool usea = isfinite(invdy.a);
				bool useb = isfinite(invdy.b);
		  
				if(!(usea|useb)){	
					if(f(s.mid())==0)ret.push_back(s); 
					return;
				}
				

				for(int i=1;i<=4;++i)
				{
					x = s.mid();
					y = f(x);

					//cout<<"\ns="<<s<<" f(x)="<<y<<" f'(s)="<<dy<<endl;
					//cout<<"s.width()="<<s.width()<<endl;
					++iter;
					//cout<<"NO."<<iter<<endl;

					r.a=-HUGE_VAL; r.b=HUGE_VAL;
				
					if(usea) (y>0? r.b : r.a) = x - y * invdy.a ;
					if(useb) (y>0? r.a : r.b) = x - y * invdy.b;
				
				  //cout<<"->r="<<r<<endl;
				
					//it doesn't change any more
					if(s.b==r.a || s.a==r.b)
					{
						ret.push_back(s);return;
					}
				
				
					bool salera = (s.a<=r.a),
							 rblesb = (r.b<=s.b); 
				
					//a split
					if(salera & rblesb )
					{
						rec(myinterval<T>(s.a,r.a));
						s.a=r.b;
						//cout<<"split: s="<<s<<endl;
						break;
					}
				
					//no zero inside s
					if(!(salera|rblesb))
						return;
					
					//empty interval on one side of s
					if(salera)
						s.b=r.a;
					else
						s.a=r.b;
				
				}			
			
			}
			else //if dy doesnt intersect 0
			{

				invdy.a=1./dy.a;
				invdy.b=1./dy.b;
							
				for(int i=1;i<=4;++i)
				{
					x=s.mid();
					y=f(x);

					//cout<<"\ns="<<s<<" f(x)="<<y<<" f'(s)="<<dy<<endl;
					//cout<<"s.width()="<<s.width()<<endl;
					++iter;
					//cout<<"NO."<<iter<<endl;
			
					if(y==0)
					{ret.push_back(s);return;}
				
					r.a = x - y * invdy.a;
					r.b = x - y * invdy.b;
					if(r.a>r.b){x=r.a;r.a=r.b;r.b=x;}
			
				  //cout<<"r="<<r<<endl;

					r.intersect(s);
				  //cout<<"r^s="<<r<<endl;


					if(r.a>r.b)
						return;
					
			
					if(
				  r.width()>=s.width()
				  |
				  r.width()<=eps*r.mid()
				  )
				  { ret.push_back(r); return; }		
				
				
					s=r;
				}
			}
		}
	}
};


//a version that accounts for rounding errors
template<class T,class tf,class tdf>
struct solve4
{
	tf f;tdf df;
  std::vector<myinterval<T> > ret;   
  T eps;

	myinterval<T> y,dy,invdy;
	T x;

  solve4(){}
  
	solve4(std::vector< myinterval<T> > &re,
			const myinterval<T> & start	,
			T eps=1e-12,
			const tf&f_=tf(),
			const tdf&df_=tdf()
			)
			:eps(eps) ,f(f_),df(df_)
  {
  	re.swap(ret);
  	
  	ret.clear();
  	rec(start);
  	
  	re.swap(ret);  	
  }
  
  void solve(std::vector<myinterval<T> > &re, const myinterval<T>&start,T leps=1e-12 )
  {
    eps=leps;
  	ret.clear();
  	rec(start);
  	re.swap(ret);	  
  }
  
	inline void addorjoin(const myinterval<T>&s)
	{
		if(!ret.empty()&&ret.back().b == s.a)
			ret.back().join(s);
		else
			ret.push_back(s); 	
	}
  
  void rec(myinterval<T> s)
	{
		y=f(s);
		if(!y.contains(0)) return;
		myinterval<T> r;
	
		int iter=0;
		for(;;)
		{	

			if(s.width() <= (eps>0 ? eps*s.mid() : -eps)  )
			{	
				addorjoin(s);
				return;
			}

			dy  = df(s);
	
			if( !(dy.a<dy.b) )
			{	ret.push_back(s); return;}
	
	
			if(dy.contains(0))
			{
				//round it so inwards, so the resulting  interval gets as small as possible
				//round it up:
				invdy.a=-(1./-dy.a);				
				//this is rounded down.
				invdy.b=1./dy.b;

				bool usea = isfinite(invdy.a);
				bool useb = isfinite(invdy.b);
		  
				if(!(usea|useb)){	
					if(f(myinterval<T>(s.mid())).contains(0))
						addorjoin(s);
					return;
				}
				

				for(int i=1;i<=1;++i)
				{
					x = s.mid();
					y = f(myinterval<T>(x));

					//cout<<"\ns="<<s<<" f(x)="<<y<<" f'(s)="<<dy<<endl;
					//cout<<"s.width()="<<s.width()<<endl;
					++iter;
					//cout<<"NO."<<iter<<endl;

					r.a=-HUGE_VAL; r.b=HUGE_VAL;
				
					//always use the absolute-value-wise smaller value of y,
					//so the exclusion interval gets as small as possible
					if(y.a>=0)
					{ 
						if(useb) r.a = x - y.a * invdy.b ;
						if(usea) r.b = x + y.a * (-invdy.a) ;
					}
					else if(y.b<=0)
					{ 
						if(usea) r.a = x - y.b * invdy.a ;
						if(useb) r.b = x + y.b * (-invdy.b) ;
					}
					else r.a=r.b=x;

					if(r.a>r.b)r.a=r.b=r.mid();

				
				  //cout<<"->r="<<r<<endl;
				
					//it doesn't change any more
					if(s.b==r.a || s.a==r.b)
					{
						if(f(s).contains(0))
							addorjoin(s);
						return;
					}
				
				
					bool salera = (s.a<=r.a),
							 rblesb = (r.b<=s.b); 
				
					//no zero inside s
					if(!(salera|rblesb))
						return;

					//a split
					if(salera & rblesb )
					{
						rec(myinterval<T>(s.a,r.a));
						s.a=r.b;
						//cout<<"split: s="<<s<<endl;
						break;
					}
				
					
					//empty interval on one side of s
					if(salera)
						s.b=r.a;
					else
						s.a=r.b;
				
				}			
			
			}
			else //if dy doesn't intersect 0 function is strictly monotonous
			{

				//round outwards
				invdy.b = -(1./(-dy.a));
				invdy.a = 1./(dy.b);				

				bool usea=isfinite(invdy.a);
				bool useb=isfinite(invdy.b);

				//if derivatives were so small that they lead to infinities
				if(!(usea|useb)){	
					if(f(myinterval<T>(s.mid())).contains(0))
							addorjoin(s);
					return;
				}
							
				for(int i=1;i<=2;++i)
				{
					x=s.mid();
					y=f(myinterval<T>(x));

					//cout<<"\ns="<<s<<" f(x)="<<y<<" f'(s)="<<dy<<endl;
					//cout<<"s.width()="<<s.width()<<endl;
					++iter;
					//cout<<"NO."<<iter<<endl;
							
					r = -y*invdy + x;
				  //cout<<"r="<<r<<endl;

					r.intersect(s);
				  //cout<<"r^s="<<r<<endl;

					if(!(r.a<=r.b))
						return;

					if(y.contains(0) && (usea&useb) )
					{
						addorjoin(r);
						return;											
					}

					
			
					if(
				  r.width()>=s.width()
				  |
					r.width()<=(eps>0? eps*r.mid() :-eps)
				  )
					{ 
						if(f(r).contains(0))
							addorjoin(r);
						return; 
					}		
				
				
					s=r;
				}
			}
		}
	}
};


#endif

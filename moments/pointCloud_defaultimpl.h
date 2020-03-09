/*
	partial specializations for transformations of moment tensors

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

/** here go default implementations of 
functions from pointCloud.h
where better spcialized functions are available
*/
/* translate an A of 3d total symmetric tensors */

template<int N>
translate<N>::translate(mytype*newA, const mytype*A, const mytype *t)
{
	mytype at[(((N+6)*N+11)*N)/6 +1];
	memset(at,0,sizeof(at));
	//calculate the tensors 1,t,t x t t x t x t , ... into at
	add3dcoords2<N>(at,t);


	for(int I=0;I<=N;++I)
	{
		for(int J=0;J<=N-I;++J)
		{
			for(int K=0;K<=N-I-J;++K, ++newA)
			{
				mytype ret = 0;
				int nu=1; //(n over k) * (n over i) * (n over j)
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
						
							ret += nu * A[ calcIndex(N,i,j,k) ] *at[ calcIndex(N,I-i,J-j,K-k) ];
							
						}					
					}							
				}

				*newA=ret;

				assert(nu==1);
			}
		}					
	}
}



template<int N>
add3dcoords2<N>::
add3dcoords2(mytype *A,const mytype *x)
{
	mytype a,b,c;

	*(A++)+=1; 

	c=x[2];
	*(A++) += c; 
	for(int k=N-1;k>0;--k)
	{
		c*=x[2];
		*(A++) += c;
	}		

	b=x[1];
	*(A++) += b;
	
	c=b;
	for(int k=N-1;k>0;--k)
	{
		c*=x[2];
		*(A++) += c;
	}		

	for(int j=N-1;j>0;--j)
	{
			b*=x[1];
			*(A++) += b;
			c=b;
			for(int k=j-1;k>0;--k)
			{
				c*=x[2];
				*(A++) += c;
			}		
	}			

	a=x[0];
	*(A++)+=a; 

	c=a;
	for(int k=N-1;k>0;--k)
	{
		c*=x[2];
		*(A++) += c;
	}		

	b=a;
	for(int j=N-1;j>0;--j)
	{
			b*=x[1];
			*(A++) += b;
			c=b;
			for(int k=j-1;k>0;--k)
			{
				c*=x[2];
				*(A++) += c;
			}		
	}			


	for(int i=N-1;i>0;--i)
	{
		a*=x[0];
		*(A++)+=a;
		c=a;
		for(int k=i-1;k>0;--k)
		{
			c*=x[2];
			*(A++) += c;
		}		

		b=a;
		for(int j=i-1;j>0;--j)
		{
			b*=x[1];
			*(A++)+=b;
			c=b;
			for(int k=j-1;k>0;--k)
			{
				c*=x[2];
				*(A++) += c;
			}		
		}			
	}
}
};




/**
* evaluate the folding of  the Tensor of order ord
* in the set of the tensor with max order N
* with the vector x
* e.g order 3: sum_ijk x_i x_j x_k A(3)_ijk
*/
template<class T,int ord,int N>
T evalTotalFold(const T*A,const T*x)
{
	T px,py,ret;
	
	T z[ord]; //z^index
	z[0] = x[2]; 
	for(int i=1;i<ord;++i)
		z[i]=z[i-1]*x[2];

	ret = z[ord-1]* A[calcindex(N,0,0,ord)];
	
	int ordoverny=1;
	py = x[1];
	for(int ny=1;ny<ord;++ny,py *= x[1])
	{
			ordoverny = ordoverny * (ord-ny+1) /ny;
			ret+= ordoverny * py * z[ord-ny-1] *A[calcIndex(N,0,ny,ord-ny)];
	}	
	ret+= py * A[calcIndex(N,0,ord,0)];

	px = x[0];
	int ordovernx=1;
	for(int nx=1;nx<ord;nx++,px*=x[0])
	{
		ordovernx = ordovernx * (ord-nx+1) /nx;
		int ordovernxny=ordovernx;
		py=px;
		for(int ny=0;ny<ord-nx;++ny,py *= x[1];)
		{
			ret+=	pxy * z[ord-nx-ny-1] * ordovernxny 
				* A[calcIndex(N,nx,ny,ord-nx-ny)];
			
			ordovernxny = ordovernxny * ( ord-nx-ny+1 )/ ny;
		}
		ret+= py * ordovernx * A[calcIndex(N,nx,ord-nx,0)];
	}

	ret+= px*A[calcIndex(N,ord,0,0)];

	return ret;
}


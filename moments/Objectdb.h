/*
    Object database using sets of moment invariants as representation.

	This file is part of the source code used for the calculation of the moment invariants
	as described in the dissertation of Max Langbein
	https://nbn-resolving.org/urn:nbn:de:hbz:386-kluedo-38558

	Copyright (C) 2020 TU Kaiserslautern,
	              Prof.Dr. Hans Hagen (AG Computergraphik und HCI) (hagen@cs.uni-kl.de)

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

#ifndef OBJECTDB_H
#define OBJECTDB_H
#ifndef GNUCC
#pragma once
#endif

#include<map>

#ifndef OLD_VC

	#include<unordered_map>
	#include<array>

	namespace std{
	template<>
	struct hash< std::array<int,28> >
	{
		std::size_t operator() (const std::array<int,28>&a ) const noexcept
		{
				size_t r=0;
				for(int i=0;i<28;++i)
					r=r*13+a[i];
				return r;
		};
	};
	}

	#define _HASH_MAP_ std::unordered_map
	//__gnu_cxx::hash_map
#else
	#include<hash_map>
	#define _HASH_MAP_ std::hash_map
#endif

#include<string>
#include<stdio.h>
#include<iostream>


#include "Momkdtree.h"
#include "momentSet.h"
#include "vecstats.h"
#include "outputflow.h"
#include "tensorSet.h"

class TsTree;


class Objectdb
{

public:
	typedef unsigned int uint32;

	typedef Momkdtree::mynum mynum;
	typedef faceMomSet tmom;
	//for debug output if needed.
	static std::ostream * & dbgout;

	struct stats{

	private:
		mynum avg,var;
	
	public:
		mynum min,max;

		//number of values
		int n;

		stats(){n=0;};

		//add values
		void add(mynum v)
		{
			if(n==0){
				min=max=v; avg=var=0; 			
			}
			else
			{
				if(min>v)min=v; 
				if(max<v)max=v;
				avg+=v;var+=v*v;
			}
			++n;
		}

		void operator +=(const stats&o){
			if(n==0)*this=o;
			else{
				avg+=o.avg;var+=o.var;
				if(o.min<min)min=o.min;
				if(o.max>max)max=o.max;
				n+=o.n;
			}
		}
		
		mynum average() const {return avg/n;}
		mynum variance()const { return (var -avg *avg/n)/(n-1);}


	};

	struct obj
	{
		//the ids of the moment sets belonging to this object
		typedef std::map<mynum,std::vector<int> > mapmoms;
		mapmoms momentSets;

		typedef vecstatscollector<tmom::n,mynum> stlist;
		//statistics only for this Object
		typedef std::map<mynum,stlist > mapstats;
		mapstats momstats;

		std::string name;
	};

	struct momsetMeta
	{
		mynum radius;
		mynum center[3];
		int objectid;
	};


	//for visualization & statistics of our method
	struct match{
		momsetMeta testmom;
		momsetMeta dbmom;
		mynum dist;
		mynum seconddist;	
	};



	void getradii(mynum&rmin_,mynum&rmax_) const
	{rmin_=rmin;rmax_=rmax;}

	//a kdtree of the moment sets of the current Objects stored.
	//is empty if momentSets are Not empty.
	struct momPerRadius{
		Momkdtree tree;
		//if tree is not initialized: the moment sets.
		std::vector<mynum> momentSets;

		//metainformation for every moment set
		std::vector<momsetMeta> meta;

		//if more than one object per moment set(i.e. reduce has been called): the additional object ids per moment set id,
		//and their distance to the original one in the distance measure that comes from pca of the currently inserted set.
		typedef std::vector<std::pair<int,mynum> > tvpim;
		typedef _HASH_MAP_<int,tvpim> tobjectsnearby;
		tobjectsnearby  objectsnearby;

		//for every moment set component, the statistics
		//stats momstats[tmom::n];


		//the current scale of the moments stored in the tree
		//mynum scale[tmom::n];

		//the covariance matrix + other statistics
		vecstatscollector<tmom::n,mynum> vstats;

		//matrix to convert moments to tree format.
		mynum convert[tmom::n][tmom::n];

		//number of independent dimensions of moments. now always set to tmom::n
		int treedim;

		void tosearchable(mynum*out,const mynum*in) const;
		//euclidean distance after coord. transform
		mynum distance(const mynum*a,const mynum*b) const;

		//init convert matrix +treedim from vs (get its with vstats.getstatscov())
		void initConvert(const vecstats<tmom::n,mynum>& vs);

		//initialize tree from moment sets if not yet done
		const Momkdtree & getTree();
	
		//restore moment sets and remove tree scale.
		std::vector<mynum> & getMomSets();	
		
		//reduce number of moments. calls getMomentsets.
		void reduce(Objectdb::mynum clusterSize);	

		void save(FILE*f);
		void load(FILE*f);
	};

	typedef std::map<mynum,momPerRadius> mapmom;

	const std::vector<obj>&getObjects()const {return objects;}
	const mapmom & getMomsets() const {return moments;}

	//build tree data structures in db
	void buildTrees();

private:

	//for every radius a structure of type momPerRadius
	mapmom moments;

	std::vector<obj> objects;

	//update the information in objects from moments
	void updateObjects();

	//range of the radii 
	//of the balls the moments were calculated from.
	mynum rmin,rmax;


public:


		//should tsets from the border be sorted out
		//in getMoms ?
		static bool dosortouttsets;

		/*
		get the invariant moments for 
		balls with radius ra around all
		leaves of the octree of the points of the pointcloud
		describing the object
		\retval m
			a vector containing for every ball the invariant moments
		\retval pt 
			a vector containing the center of gravity for every pointset
			from which a moment set was computed
		\retval stats:
		 a set of statistics (variance, mean,...)
		 for every component of the moment sets
		 \param t: 
		 an octree of the pointcloud having the tensorset
		 0A,1A,2A,3A,4A of the points inside the node 
		 stored in each node
		 \param ra:
			the radius of the balls
		 \param minh:
			minimum height at which to stop the 
			summing up of the tensorsets for each ball
	*/
	static int getMoms(
		//the invariant moments
		std::vector<mynum>& m,
		//the corresponding points
		std::vector<mynum>& pt,
		//statistik
		vecstatscollector<tmom::n,mynum> &s,
		//Baum, der die tensorsets enthaelt				
		const TsTree&t,
		//radius of the balls
		mynum ra,
		//minimum height to go down in the tree
		int minh=0
	);


	/*add an object to the database
	*( it is assumed that its in the same scale/coordinates
	*  as the other objects )
	*\param points: the pointcloud describing that object
	* should we go through ?
	*\param clusterit: should moments too near in one object be
	*thrown out
	* be thrown out ?
	*/
	const obj& addObject(const std::string&name,
								const std::vector<mynum> & points,								
								bool clusterit=false
								);

	/** reduce number of moments in db.
	*(only one moment per clusterSize part of std deviation)
	*\param clusterSize: size of one bucket 
	* as fraction of the std deviation
	*/
	void reduce(Objectdb::mynum clusterSize);

	/**
	* get the likelihood for the different objects of being contained in 
	*the pointcloud given.
	*
	*\param ret:
	*vector with a probability for each object
	*being the searched one 
	*(if sum can is < 1  it is likely that no known object
	* is inside, if it is > 1 it is likely 
	*that there are multiple known objects in the pointcloud
	*
	*\param points:
	* the coordinates of the points of the pointcloud
	*
	*\param matches:
	*this vector will be filled if not zero.
	*if u want to do stats: the matches that were found
	*and that contribute to the likelihood	
	*\param classes: if not zero, gives an object class for every object.
	*ret will give probabilities for the classes instead of the object.
	*/
	void findObjects(
		std::vector<mynum> &ret,
		const std::vector<mynum> & points		
		,std::vector<match> * matches=0
		,const std::vector<int> * classes=0
		,bool ignoreRadii=false
	);

	int nmoments();


	void save(const char*name);

	void load(const char*name);


	Objectdb(void);
	~Objectdb(void);

	static void kdtreeEntryFromTset(
	mynum* ret,
	faceMomSet & bufm,
	const tensorSet3D<> &ts,
	mynum radius
	);

	static int getInv( double ra,	const std::vector<mynum> & points, const TsTree&t, mynum &pnpr,
		mynum*m);



};

#endif

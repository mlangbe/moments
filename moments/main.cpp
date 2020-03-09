/*
test the things.

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

#include<fstream>
#include "pointCloud.h"
#include"moments.h"
#include "momStats.h"
#include "PointCloudIO.h"
#include "myhash.h"
#include "zindex.h"
using namespace std;

void testtensorgraph();
int vecMomAppMain(int argc, const char*argv[]);

void smoothfractaldim();

namespace di{
void testImplicits();
};
namespace gr{
	void testgroebner();
};




int main(int argc, const char* argv[])
{ 

	//myhash<int> mh;

	//smoothfractaldim();
	//vecMomAppMain(argc,argv);
	/*
	ofstream mout("basistransform.txt");
	createBasisTransform(mout,4);
	mout.close();
*/
	momStatsmain(argc,argv);

	//getMomentsScalarField(cerr,cout);

	//testtensorgraph();
	//testOptimizedFold();
	/*
	testMomFunSpeed();
	testMomFunSpeed2<35,28,0>();

	cout<<"vector functions:"<<endl;
	cout<<"optimized with quasi-horner-scheme:"<<endl;
	testMomFunSpeed2<30,27,2>();	
	cout<<"optimized for common subexpr.:"<<endl;
	testMomFunSpeed2<30,27,1>();
	*/

	/*
	ofstream of("vecmomspec_dnf.txt");
	createVecMomSpecialization(of,0,false,false,"polynomials in dnf");
	of.close();
*/

//	testSavePolynom();

/*
for(int i=0;i<argc;++i)
		printf("argv[%d]:%s\n",i,argv[i]);
		*/
//	testIntegrateTsets();
//	testgetpairs3();
// testpoly();
//return 0;
/*	
  infnum a=13892560913865,b=2385627345;
  cout<<a<<' '<<b<<' ';
  a*=b;
  cout<<a<<endl;
  a*=b;
  cout<<a<<endl;
  a*=b;
  cout<<a<<endl;
  a*=b;
  cout<<a<<endl;
  a*=b;
  cout<<a<<endl;
  a*=b;
  cout<<a<<endl;
  a*=b;
  cout<<a<<endl;
  a*=b;
  cout<<a<<endl;


	return 0;
	*/

//	testdependencytest();
//	return 0;

/*
	ofstream myout("pairingsScalarDim3Ord4.txt");
	ofstream mypolys("polysScalarDim3Ord4.txt");
	getMomentsScalarField(myout,mypolys);
	myout.close();mypolys.close();
*/
	
//	ofstream myout("pairingsVectorDim3Ord4.txt");
	//ofstream mypolys("polysVectorDim3Ord4.txt");
	

	





		
	


//	myout.close();
	//mypolys.close();

	//testgetpairs3();


	//doOptimize();



	
	//testMomFuns();
	//testadd3dc();	
	//testMomFunSpeed();
	//createTotalFolds();
//throw;
/*	
	string prefix="c:\\u\\max\\models\\da\\";

	static const int 	argc2=8;

	const char* names[8]={ 
		0,
		"dino",
		"huhn",
		"handy",
		"kaestchen2",
		"kaestchen",
		"steckdose",
		"pig"};

  string nam[argc2];
	
	const char* argv2[argc2];
	
	argv2[0]=argv[0];
	
	for(int i=1;i<argc2;++i){
		nam[i]=prefix+names[i]+".raw";
		argv2[i]=nam[i].c_str();
	}


*/
  
	//momStatsmain(argc,argv);

	//testspeedzindex();
	//testclustertree();
	//testMomkdtree();

	//di::testImplicits();

	//gr::testgroebner();
	return 0;
	
}


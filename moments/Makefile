SOURCES = $(wildcard *.cpp) $(wildcard *.cc)

#SOURCES = clustertree.cpp computeMoments.cpp  \
#infnum.cc moments.cpp momentSet.cpp Momkdtree.cpp momStats.cpp Objectdb.cpp pointCloud.cpp \
#polynom.cc symmtensor.cpp TsTree.cpp TsTree2.cpp infbruch.cpp \
#optimizedPoly.cpp PointCloudIO.cpp \
#optimizedFold.cpp \
#momentSetVectors.cpp momentSetVectors2.cpp indepPolySet.cpp \
#tensWithSymm.cpp tensorGraph.cpp tensorGraphSplitTree.cpp outputflow.cpp \
#computeMomentsPolyn.cpp 


OBJ = $(SOURCES:.cpp=.o)
OBJECTS = $(OBJ:.cc=.o)
INCLUDES = -I../ #-I../../libsvm-3.17
OPTIONS =  -w -DGNUCC -std=c++17 -pthread #-fabi-version=1 #-DGCC33
LIBS =-lpthread 
CXX = g++ #-4.8 #g++-3.3



all: a.out deps.d


deps.d: $(SOURCES)
	$(CXX) -MM $(OPTIONS) $(INCLUDES) $(SOURCES) > deps.d

include deps.d

%.o:%.cc 
	$(CXX) $(OPTIONS) $(INCLUDES) -c -o $@ $<

%.o:%.cpp 
	$(CXX) $(OPTIONS) $(INCLUDES) -c -o $@ $<

momStats.pre: momStats.cpp
	$(CXX) $(OPTIONS) $(INCLUDES)  -E momStats.cpp > momStats.pre



a.out: $(OBJECTS)
	$(CXX) $(OPTIONS) $(LIBS) $(OBJECTS)

clean: 
	rm -f $(OBJECTS) deps.d a.out gmon.out core.*

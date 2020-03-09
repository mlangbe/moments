SOURCES = $(wildcard *.cpp) $(wildcard *.cc)



OBJ = $(SOURCES:.cpp=.o)
OBJECTS = $(OBJ:.cc=.o)
INCLUDES = -I../moments 
OPTIONS =  -w -DGNUCC -std=c++11  
CXX = g++ #-4.8 #g++-3.3



all: kegellager.exe kugel.exe deps.d


deps.d: $(SOURCES)
	$(CXX) -MM $(OPTIONS) $(INCLUDES) $(SOURCES) > deps.d

include deps.d

%.o:%.cpp 
	$(CXX) $(OPTIONS) $(INCLUDES) -c -o $@ $<

%.exe:%.o Objects.o
	$(CXX) $(OPTIONS) -o $@ $< Objects.o

clean: 
	rm -f $(OBJECTS) deps.d a.out gmon.out core.*

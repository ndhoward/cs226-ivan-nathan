all: point.o particle-filter.o geomModel.o
	g++ -lglut point.o particle-filter.o geomModel.o -o point
point.o: point.cpp point.h particle-filter.h geomModel.h
	g++ -lglut -c point.cpp
particle-filter.o: particle-filter.cpp particle-filter.h
	g++ -c particle-filter.cpp
geomModel.o: geomModel.cpp geomModel.h
	g++ -c geomModel.cpp
clean:
	rm -f point *.o *~

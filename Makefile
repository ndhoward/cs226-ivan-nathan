all: point.o particle-filter.o geomModel.o
	g++ -lglut point.o particle-filter.o geomModel.o -o point
point.o: point.cpp
	g++ -lglut -c point.cpp
particle-filter.o: particle-filter.cpp
	g++ -c particle-filter.cpp
geomModel.o: geomModel.cpp
	g++ -c geomModel.cpp
clean:
	rm -f point *.o *~

# This is a makefile
all: point.o particle-filter.o geomModel.o
#	g++ -lglut -Iann/include -Lann/lib -lANN point.o particle-filter.o geomModel.o -o point
	g++ point.o particle-filter.o geomModel.o -o point -Lann/lib -lANN -lglut
point.o: point.cpp point.h particle-filter.h geomModel.h
#	g++ -lglut -Iann/include -Lann/lib -lANN -c point.cpp
	g++ -Iann/include -c point.cpp
particle-filter.o: particle-filter.cpp particle-filter.h
	g++ -c particle-filter.cpp
geomModel.o: geomModel.cpp geomModel.h
	g++ -c geomModel.cpp
clean:
	rm -f point *.o *~

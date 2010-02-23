point: point.cpp
	g++ -lglut point.cpp -o point
	g++ particle-filter.cpp -o particle
clean:
	rm -f point particle

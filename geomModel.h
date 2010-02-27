#ifndef GEOM_MODEL_H
#define GEOM_MODEL_H

#include <vector>

#ifndef PI
#define PI 3.14159265
#endif

using namespace std;

// sample point belonging to a geometric model (e.g., pedestrian, bicyclist)
typedef struct GMPoint {
	float x, y, z;
} GMPoint;

// a point of measurement type (raw Velodyne data)
typedef struct ZPoint {
	float x, y, z;
	int t;
} ZPoint;

float likelihoodPerson(vector < ZPoint > &frame,
					   float xPos, float yPos);

#endif

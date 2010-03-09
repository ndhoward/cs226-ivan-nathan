#include <cmath>
#include <algorithm>	// min
#include <iostream>
#include "geomModel.h"

// for pca
#include <Eigen/Core>
#include <Eigen/Eigen>

using namespace std;

/*
void generatePersonModelSamples(vector < GMPoint > &samples,
							    int num_samples, // apprixmate, check size of returned vector
							    float radius,
							    float height,
							    float x_ctr,
							    float y_ctr) {
	
	// sample the cylinder uniformly
	float samplesPerArea = num_samples / (2.0*PI*radius*height + 2.0*PI*radius*radius);
	float samplesPerSide = samplesPerArea * 2.0*PI*radius*height;
	float samplesPerTopBtm = samplesPerArea * PI*radius*radius;
	float samplesPerCircle = samplesPerSide/height;
	
	
	// sample the sides
	float deltaTheta = 2*PI/samplesPerCircle;
	float deltaH = height/(samplesPerSide/(2*PI*radius));
	for(float z = 0.0; z <= height; z+= deltaH)
	{
		for(float theta = 0.0; theta < 2*PI; theta += deltaTheta)
		{
			GMPoint pt;
			pt.x = radius*cos(theta) + x_ctr;
			pt.y = radius*sin(theta) + y_ctr;
			pt.z = z;
			samples.push_back(pt);
		}
	}
	
	// sample the top
	float samplesInTopSquare = samplesPerTopBtm*4.0/PI;
	float deltaXY = 2.0*radius/sqrt(samplesInTopSquare);
	cout << samplesInTopSquare / (2*radius*2*radius) << endl;
	cout << samplesPerArea << endl;
	for(float x = -radius; x <= radius; x+= deltaXY)
	{
		for(float y = -radius; y <= radius; y+= deltaXY)
		{
			if(x*x+y*y <= radius*radius)
			{
				GMPoint ptTop;
				ptTop.x = x + x_ctr;
				ptTop.y = y + y_ctr;
				ptTop.z = height;
				samples.push_back(ptTop);
				
				GMPoint ptBtm;
				ptBtm.x = x + x_ctr;
				ptBtm.y = y + y_ctr;
				ptBtm.z = 0.0;
				samples.push_back(ptBtm);
			}
		}
	}
}

float distanceToSet(vector < GMPoint > &samples,
					float x_pt, float y_pt, float z_pt)
{
	float x = samples[0].x;
	float y = samples[0].y;
	float z = samples[0].z;
	float minDist = sqrt((x-x_pt)*(x-x_pt) + (y-y_pt)*(y-y_pt) + (z-z_pt)*(z-z_pt));
	for(int i = 1; i < samples.size(); i++)
	{
		x = samples[i].x;
		y = samples[i].y;
		z = samples[i].z;
		float dist = sqrt((x-x_pt)*(x-x_pt) + (y-y_pt)*(y-y_pt) + (z-z_pt)*(z-z_pt));
		if(dist < minDist)
			minDist = dist;
	}
	return minDist;
}
*/

// float distanceToCylinder(float cylR, float cylH, // cylinder size
// 						 float cylX, float cylY, // center at bottom cylinder face
// 						 float testX, float testY, float testZ)	// point from which to measure distance to cylinder
// {
// 	// distance to sides = abs(||(x,y) - (cylX, cylY)|| - cylR)
// 	float distToSides = abs(sqrt((testX-cylX)*(testX-cylX) + (testY-cylY)*(testY-cylY)) - cylR);
// 	float distToTop = abs(testZ - cylH);
// 	float distToBtm = abs(testZ - 0.0);
// 	float distToTopEdg = sqrt(distToTop*distToTop + distToSides*distToSides);
// 	float distToBtmEdg = sqrt(distToBtm*distToBtm + distToSides*distToSides);
	
// 	// return the correct min distance depending on where the test point is
// 	if(0.0 <= testZ && testZ <= cylH)
// 		return min(min(distToSides, distToTop), distToBtm);
// 	else if(testX*testX + testY*testY <= cylR*cylR)
// 		return min(distToTop, distToBtm);
// 	else
// 		return min(distToTopEdg, distToBtmEdg);
// }

// --------------------------------------------------------------------
// returns the ratio of the largest eigenvalue to the second largest
// eigenvalue of the covariance matrix when doing PCA on points belonging
// to a blob.
//
// inspired by PCA code from 
// http://codingplayground.blogspot.com/2010/01/pca-dimensional-reduction-in-eigen.html
// --------------------------------------------------------------------
float pVal(vector < ZPoint > &blob)
{
	using namespace Eigen;
	unsigned int const m = 3;			// dimension of each point
	unsigned int n = blob.size();		// number of points
	
	MatrixXf DataPoints = MatrixXf::Zero(m,n);	// matrix (m x n)
	for(int r = 0; r < blob.size(); r++)
	{
		DataPoints(r,0) = blob[r].x;
		DataPoints(r,1) = blob[r].y;
		DataPoints(r,2) = blob[r].z;
	}
	
	double mean;
	VectorXf meanVector;
	
	// for each point, 
	// center the point with the mean along all the coordinates
	for(int i = 0; i < DataPoints.cols(); i++)
	{
		mean = (DataPoints.col(i).sum())/m;			// compute mean
		meanVector = VectorXf::Constant(m, mean);	// create vector with constant value = mean
		DataPoints.col(i) -= meanVector;			// subtract away the mean
	}
	
	// get the covariance matrix
	MatrixXf Covariance = MatrixXf::Zero(m,m);
	Covariance = (1.0/(float)n) * DataPoints * DataPoints.transpose();
	
	// compute eigenvalues of the covariance matrix
	EigenSolver<MatrixXf> m_solve(Covariance);
	VectorXf eigenvalues = VectorXf::Zero(m);
	eigenvalues = m_solve.eigenvalues().real();
	
	// return the ratio of the two largest eigenvalues
	return (eigenvalues(m-1)/eigenvalues(m-2));
}


//hack distanceToCylinder till more complex one works
float distanceToCylinder(float cylR, float cylH, // cylinder size
			 float cylX, float cylY, // center at bottom cylinder face
			 float testX, float testY, float testZ)	// point from which to measure distance to cylinder
{
  return sqrt(pow((cylX-testX),2)+pow((cylY-testY),2));
}

// not normalized
float likelihoodPerson(vector < ZPoint > &frame,
					   float xPos, float yPos) {
	
	// tunable model parameters
	const float EPSILON = 1.0/3.0;	
	const float PERSON_RADIUS = 1.0/3.0;
	const float PERSON_HEIGHT = 1.5;
	
	int countWithinEps = 0;
	for(int i = 0; i < frame.size(); i++)
	{
		float dist = distanceToCylinder(PERSON_RADIUS, PERSON_HEIGHT,
						xPos, yPos,	// person center
						frame[i].x, frame[i].y, frame[i].z);
		
		if(dist <= EPSILON)
			countWithinEps++;
	}
	return (float)countWithinEps;
}

/*
int main(int argc, char *argv[])
{
	// test for person samples
	vector<GMPoint> personSet;
	generatePersonModelSamples(personSet, 1000,
							   2, 5, 0, 0);
	cout << "#Num person samples: " << personSet.size() << endl;
	for(int i = 0; i < personSet.size(); i++)
	{
		//cout << personSet[i].x << ',' << personSet[i].y << ',' << personSet[i].z << endl;
	}
	return 0;
}
*/






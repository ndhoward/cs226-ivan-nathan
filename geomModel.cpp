#include <cmath>
#include <algorithm>	// min, max
#include <iostream>
#include "geomModel.h"

// for pca
#include <Eigen/Core>
#include <Eigen/Eigen>

using namespace std;

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
	// import most common Eigen types 
	USING_PART_OF_NAMESPACE_EIGEN
	
	unsigned int n = blob.size();	// number of points
	unsigned int m = 3;				// dimension of each point
	
	MatrixXf DataPoints = MatrixXf::Zero(m,n);	// matrix (m x n)
	
	for(int i = 0; i < n; i++)
	{
		DataPoints(0,i) = blob[i].x;
		DataPoints(1,i) = blob[i].y;
		DataPoints(2,i) = blob[i].z;
	}
	
	//double mean;
	//VectorXf meanVector;
	// for each point, 
	// center the point with the mean along all the coordinates
	//for(int i = 0; i < DataPoints.cols(); i++)
	//{
	//	mean = (DataPoints.col(i).sum())/m;			// compute mean
	//	meanVector = VectorXf::Constant(m, mean);	// create vector with constant value = mean
	//	DataPoints.col(i) -= meanVector;			// subtract away the mean
	//}
	
	double mean;
	RowVectorXf meanVector;
	for(int i = 0; i < DataPoints.rows(); i++)
	{
		mean = (DataPoints.row(i).sum())/n;
		meanVector = RowVectorXf::Constant(n,mean);
		DataPoints.row(i) -= meanVector;
	}
	
	// get the covariance matrix
	MatrixXf Covariance = MatrixXf::Zero(m,m);
	Covariance = (1.0/(float)n) * DataPoints * DataPoints.transpose();
	
	// compute eigenvalues of the covariance matrix
	EigenSolver<MatrixXf> m_solve(Covariance);
	VectorXf eigenvalues = VectorXf::Zero(m);
	eigenvalues = m_solve.eigenvalues().real();
	
	// largest and second largest eigenvalues
	// http://stackoverflow.com/questions/1582356/fastest-way-of-finding-the-middle-value-of-a-triple/1582406#1582406
	float a = eigenvalues(0);
	float b = eigenvalues(1);
	float c = eigenvalues(2);
	float midev = (a <= b) 
    				? ((b <= c) ? b : ((a < c) ? c : a)) 
    				: ((a <= c) ? a : ((b < c) ? c : b));
	float maxev = max(max(a,b),c);
	
	// return the ratio of the two largest eigenvalues
	return (maxev/midev);
}


//hack distanceToCylinder till more complex one works
float distanceToCylinder(float cylR, float cylH, // cylinder size
			 float cylX, float cylY, // center at bottom cylinder face
			 float testX, float testY, float testZ)	// point from which to measure distance to cylinder
{
  return sqrt(pow((cylX-testX),2)+pow((cylY-testY),2));
}

float likelihoodPerson(vector <ZPoint> &blob,
					   float xPos, float yPos) {
	// p-value mean and std deviation
	const float MEANP = 7.04;
	const float SIGP  = 2.80;	
	
	// P-value likelihood
	float pV = pVal(blob);
	//float pVlike = 1.0/abs((pV - MEANP)/SIGP + SIGP/10000.0); // regularize to make sure not dividing by zero
	float pVlike = 20.0/sqrt(2*PI*SIGP*SIGP)*exp(-pow(pV-MEANP,2)/(2*SIGP*SIGP));
	
	// density likelihood
	// tunable model parameters
	const float EPSILON = 1.0/3.0;	
	const float PERSON_RADIUS = 1.0/3.0;
	const float PERSON_HEIGHT = 1.5;
	
	int countWithinEps = 0;
	for(int i = 0; i < blob.size(); i++)
	{
		float dist = distanceToCylinder(PERSON_RADIUS, PERSON_HEIGHT,
						xPos, yPos,	// person center
						blob[i].x, blob[i].y, blob[i].z);
		
		if(dist <= EPSILON)
			countWithinEps++;
	}
	float distLike = (float)countWithinEps/100.0;
	
	return pVlike*distLike;
}

float likelihoodBike(vector <ZPoint> &blob,
					 float xPos, float yPos) {
	// p-value mean and std deviation
	const float MEANP = 2.01;
	const float SIGP  = 0.33;
	
	// P-value likelihood
	float pV = pVal(blob);
	//float pVlike = 1.0/abs((pV - MEANP)/SIGP + SIGP/10000.0); // regularize to make sure not dividing by zero
	float pVlike = 1.0/sqrt(2*PI*SIGP*SIGP)*exp(-pow(pV-MEANP,2)/(2*SIGP*SIGP));
	
	// density likelihood
	// tunable model parameters
	const float EPSILON = 1.0/3.0;	
	const float PERSON_RADIUS = 1.0/1.0;
	const float PERSON_HEIGHT = 2.0;
	
	int countWithinEps = 0;
	for(int i = 0; i < blob.size(); i++)
	{
		float dist = distanceToCylinder(PERSON_RADIUS, PERSON_HEIGHT,
						xPos, yPos,	// person center
						blob[i].x, blob[i].y, blob[i].z);
		
		if(dist <= EPSILON)
			countWithinEps++;
	}
	float distLike = (float)countWithinEps/100.0;
	
	return pVlike*distLike;
}

/*
// not normalized
float likelihoodPerson(vector < ZPoint > &blob,
					   float xPos, float yPos) {
	// tunable model parameters
	const float EPSILON = 1.0/3.0;	
	const float PERSON_RADIUS = 1.0/3.0;
	const float PERSON_HEIGHT = 1.5;
	
	int countWithinEps = 0;
	for(int i = 0; i < blob.size(); i++)
	{
		float dist = distanceToCylinder(PERSON_RADIUS, PERSON_HEIGHT,
						xPos, yPos,	// person center
						blob[i].x, blob[i].y, blob[i].z);
		
		if(dist <= EPSILON)
			countWithinEps++;
	}
	return (float)countWithinEps;
}

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






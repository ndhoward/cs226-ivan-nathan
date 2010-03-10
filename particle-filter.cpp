#include <cmath>
#include <cstdlib>
#include <cstdio>

#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <boost/math/distributions/normal.hpp>

#include "particle-filter.h"

using namespace std;

Particle_filter::Particle_filter(float max_x, float min_x, float max_y, float min_y) {
  particles = new Particle[NUM_PARTICLES];
  oldParticles = new Particle[NUM_PARTICLES];
  weights = new float[NUM_PARTICLES];
  cumulativeWeightsIndex = new float[NUM_PARTICLES];
  /* initialize random seed: */
  srand ( time(NULL) );
  maxX = max_x;
  minX = min_x;
  maxY = max_y;
  minY = min_y;
  float xRange = maxX-minX;
  float yRange = maxY-minY;
  for(int i=0; i<NUM_PARTICLES; i++) {
    particles[i].x = rand()*xRange/((float)RAND_MAX) + minX;
    particles[i].y = rand()*yRange/((float)RAND_MAX) + minY;
    particles[i].theta = 2*PI*((float)rand())/RAND_MAX;
    //particles[i].theta = PI/2;
    //cout << "initial points: " << particles[i].x << ", " << particles[i].y << ", theta: " << particles[i].theta << endl;
  }
}

Particle_filter::~Particle_filter() {
  delete particles;
  delete oldParticles;
  delete weights;
  delete cumulativeWeightsIndex;
}

void Particle_filter::bound_particle(Particle &p) {
  if (p.x > maxX) {
    p.x = maxX;
  } else if (p.x < minX) {
    p.x = minX;
  }
  if (p.y > maxY) {
    p.y = maxY;
  } else if (p.y < minY) {
    p.y = minY;
  }
}

// function modified to make sure inheritance works correctly
// this should never be called
void Particle_filter::update_Xt(Particle &p) {
  //p.x = p.x + cos(p.theta)*motionModelStepSize;
  //p.y = p.y + sin(p.theta)*motionModelStepSize;
  //bound_particle(p);
  p.x = 0;
  p.y = 0;
}


float Particle_filter::setup_importance_sample() {
  float sum = 0;
  for (int i=0; i<NUM_PARTICLES; i++) {
    sum += weights[i];
    cumulativeWeightsIndex[i] = sum;
  }
  return sum;
}

void Particle_filter::updateMeasurments(vector< vector< ZPoint > > *blobs) {
  connectedComponents = blobs;
}

void Particle_filter::setConnectedComponentMeans(vector< ZPoint > *blobMeans) {
  connectedComponentsMean = blobMeans;
}


float distance(ZPoint &point, Particle &particle) {
  float xval = point.x - particle.x;
  float yval = point.y - particle.y;
  return xval*xval+yval*yval;
}

vector< ZPoint > *Particle_filter::findClosestConnectedCompoenent(Particle &p) {
  int idx = -1;
  float minDist = 1000;
  for (int i=0; i<connectedComponentsMean->size(); i++) {
    ZPoint mean = connectedComponentsMean->at(i);
    float dist = distance(mean, p);
    if (dist < minDist) {
      minDist = dist;
      idx = i;
    }
  }
  return &connectedComponents->at(idx);
}

//this function should never be called
float Particle_filter::likelihood(Particle &p) {
  p.likelihood = -666;
  return -666;
}


// update takes in X_{t-1},u_t,z_t and returns X_t
void  Particle_filter::update() {
  // swap particles and oldParticles lists
  Particle *temp = particles;
  particles = oldParticles;
  oldParticles = temp;
 
  for (int i=0; i<NUM_PARTICLES; i++) {
    // sample x_t^m ~ p(x_t | u_t, x_{t-1}^m)
    update_Xt(oldParticles[i]);
    // w_t^m = p(z_t | x_t^m)
    weights[i] = likelihood(oldParticles[i]);
    //cout << "weights["<< i << "] ("<< oldParticles[i].x << "," << oldParticles[i].y << "," << oldParticles[i].theta << "): computed by likelihood as: " << weights[i] << endl;
  }

  float cumulative = setup_importance_sample();
  
  for (int i=0; i<NUM_PARTICLES; i++) {
    // draw i with probability ~ w_t^i
    float weightIdx = rand()*(cumulative/RAND_MAX);
    //cout << "weightIdx: " << weightIdx << endl;
    int j = binarySearch(cumulativeWeightsIndex, 0, NUM_PARTICLES-1, weightIdx);
    // add x_t^i to X_t
    particles[i] = oldParticles[j];
    jiggle_particle(particles[i]);
  }

  //randomly resample 1/3 of the particles
  //random_resample(0.1);
}

// For now this is a uniform probabilty centered at the particle
void Particle_filter::jiggle_particle(Particle &p) {
  const float POSJIGGLE=0.04; // in meers
  const float THETAJIGGLE=90.0; // in degrees
  p.x += POSJIGGLE*(float)rand()/RAND_MAX;
  p.y += POSJIGGLE*(float)rand()/RAND_MAX;
  p.theta = fmodf(p.theta - THETAJIGGLE/2.0*PI/180.0 + THETAJIGGLE*PI/180.0*((float)rand())/RAND_MAX, 2*PI);
}

// Randomly resample some fraction of the particles with a uniform distribution
void Particle_filter::random_resample(float prob) {
  float xRange = maxX-minX;
  float yRange = maxY-minY;
  for (int i=0; i<NUM_PARTICLES; i++) {
    float t = ((float)rand())/RAND_MAX;
    if (t <= prob) {
      particles[i].x = rand()*xRange/((float)RAND_MAX) + minX;
      particles[i].y = rand()*yRange/((float)RAND_MAX) + minY;
      particles[i].theta = 2*PI*((float)rand())/RAND_MAX;
    }
  }
}

Particle *Particle_filter::getParticles() {
  return particles;
}


int Particle_filter::binarySearch(float sortedArray[], int first, int last, float key) {
  int origFirst = first;
  int origLast = last;
   while (first <= last) {
     if (first == last)
       return first;
     int mid = (first + last) / 2;  // compute mid point.
     if (key > sortedArray[mid])
       first = mid + 1;  // repeat search in top half.
     else if (mid == 0) // special case where 0<= key <= sortedArray[0]
       return mid;
     else if (key <= sortedArray[mid-1]) 
       last = mid-1; // repeat search in bottom half.
     else
       return mid;     // found it. return position
   }
   cout << "binarySearch failed..." << endl;
   cout << "first: "<<first<<" last: "<<last<<" key: "<<key<<endl;
   for (int i=origFirst; i<=origLast;i++){
     cout<< "sortedArray["<<i<<"]: "<<sortedArray[i]<<endl;
   }
   return -(first + 1);    // failed to find key
}


// Person functions
Person_filter::Person_filter(float max_x, float min_x, float max_y, float min_y) : Particle_filter(max_x, min_x, max_y, min_y) { }

float Person_filter::likelihood(Particle &p) {
  vector<ZPoint> *blob = findClosestConnectedCompoenent(p);
  float likelihood = likelihoodPerson(*blob, p.x, p.y);
  p.likelihood = likelihood;
  return likelihood;
}


void Person_filter::update_Xt(Particle &p) {
  p.x = p.x + cos(p.theta)*0.17;
  p.y = p.y + sin(p.theta)*0.17;
  //bound_particle(p);
}

void Person_filter::jiggle_particle(Particle &p) {
  const float POSJIGGLE=0.04; // in meters
  const float THETAJIGGLE=30.0; // in degrees
  p.x += POSJIGGLE*(float)rand()/RAND_MAX;
  p.y += POSJIGGLE*(float)rand()/RAND_MAX;
  p.theta = fmodf(p.theta - THETAJIGGLE/2.0*PI/180.0 + THETAJIGGLE*PI/180.0*((float)rand())/RAND_MAX, 2*PI);
}

// Bike functions
Bike_filter::Bike_filter(float max_x, float min_x, float max_y, float min_y) : Particle_filter(max_x, min_x, max_y, min_y) { }

float Bike_filter::likelihood(Particle &p) {
  //find the nearest blob
  vector<ZPoint> *blob = findClosestConnectedCompoenent(p);
  float likelihood = likelihoodBike(*blob, p.x, p.y);
  p.likelihood = likelihood;
  return likelihood;
}

void Bike_filter::update_Xt(Particle &p) {
  p.x = p.x + cos(p.theta)*0.32;
  p.y = p.y + sin(p.theta)*0.32;
  //bound_particle(p);
}


void Bike_filter::jiggle_particle(Particle &p) {
  const float POSJIGGLE=0.15; // in meters
  const float THETAJIGGLE=30.0; // in degrees
  p.x += POSJIGGLE*(float)rand()/RAND_MAX;
  p.y += POSJIGGLE*(float)rand()/RAND_MAX;
  p.theta = fmodf(p.theta - THETAJIGGLE/2.0*PI/180.0 + THETAJIGGLE*PI/180.0*((float)rand())/RAND_MAX, 2*PI);
}

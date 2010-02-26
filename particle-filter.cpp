#include <cmath>
#include <cstdlib>
#include <cstdio>

#include <vector>
#include <fstream>
#include <iostream>
#include <string>

#include "particle-filter.h"

using namespace std;

Particle_filter::Particle_filter() {
  particles = new struct particle[NUM_PARTICLES];
  oldParticles = new struct particle[NUM_PARTICLES];
  weights = new float[NUM_PARTICLES];
  cumulativeWeightsIndex = new float[NUM_PARTICLES];
  /* initialize random seed: */
  srand ( time(NULL) );
  for(int i=0; i<NUM_PARTICLES; i++) {
    particles[i].x = rand() * ((MAX_X-MIN_X)/RAND_MAX)*rand()+MIN_X;
    particles[i].y = rand() * ((MAX_Y-MIN_Y)/RAND_MAX)*rand()+MIN_Y;
  }
}

Particle_filter::~Particle_filter() {
  delete particles;
  delete oldParticles;
  delete weights;
  delete cumulativeWeightsIndex;
}

void Particle_filter::bound_particle(particle &p) {
    if (p.x > MAX_X) {
      p.x = MAX_X;
    } else if (p.x < MIN_X) {
      p.x = MIN_X;
    }
    if (p.y > MAX_Y) {
      p.y = MAX_Y;
    } else if (p.y < MIN_Y) {
      p.y = MIN_Y;
    }
}

//TODO: this is crap, replace!
//Dummy motion model - make realistic later
void Particle_filter::update_Xt(particle &p) {
  p.x = p.x + 0.5;
  p.y = p.y + 0.5;
}

float Particle_filter::setup_importance_sample() {
  float sum = 0;
  for (int i=0; i<NUM_PARTICLES; i++) {
    sum += weights[i];
    cumulativeWeightsIndex[i] = sum;
  }
  return sum;
}

void Particle_filter::updateMeasurments(vector<point> *curFrame) {
  measurments = curFrame;
}

// fill this in with Ivan's code
float Particle_filter::likelihood(float blah, struct particle p) {
  return p.x * p.y;
}


// update takes in X_{t-1},u_t,z_t and returns X_t
void  Particle_filter::update() {
  // swap particles and oldParticles lists
  struct particle *temp = particles;
  particles = oldParticles;
  oldParticles = temp;

  for (int i=0; i<NUM_PARTICLES; i++) {
    // sample x_t^m ~ p(x_t | u_t, x_{t-1}^m)
    update_Xt(oldParticles[i]);
    // w_t^m = p(z_t | x_t^m)
    weights[i] = likelihood(13, oldParticles[i]);
  }

  float cumulative = setup_importance_sample();
  
  for (int i=0; i<NUM_PARTICLES; i++) {
    // draw i with probability ~ w_t^i
    float weightIdx = rand()*(cumulative/RAND_MAX);
    int j = binarySearch(weights, 0, NUM_PARTICLES-1, weightIdx);
    // add x_t^i to X_t
    particles[i] = oldParticles[j];
  }
}

//TODO: implement this later - probably using BOOST library
void Particle_filter::jiggle_particle(particle &p) {
  
}

// Taken from: http://www.fredosaurus.com/notes-cpp/algorithms/searching/binarysearch.html
int Particle_filter::binarySearch(float sortedArray[], int first, int last, float key) {
   // function:
   //   Searches sortedArray[first]..sortedArray[last] for key.  
   // returns: index of the matching element if it finds key, 
   //         otherwise  -(index where it could be inserted)-1.
   // parameters:
   //   sortedArray in  array of sorted (ascending) values.
   //   first, last in  lower and upper subscript bounds
   //   key         in  value to search for.
   // returns:
   //   index of key, or -insertion_position -1 if key is not 
   //                 in the array. This value can easily be
   //                 transformed into the position to insert it.
   
   while (first <= last) {
       int mid = (first + last) / 2;  // compute mid point.
       if (key >= sortedArray[mid+1]) 
           first = mid + 1;  // repeat search in top half.
       else if (key < sortedArray[mid]) 
           last = mid - 1; // repeat search in bottom half.
       else
           return mid;     // found it. return position /////
   }
   return -(first + 1);    // failed to find key
}

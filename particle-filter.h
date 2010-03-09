#ifndef PARTICLE_FILTER_H
#define PARTICLE_FILTER_H

#include <vector>
#include "geomModel.h"

#define NUM_PARTICLES 300

typedef struct Particle {
  float x;
  float y;
  float theta;
  float likelihood;
} Particle;

class Particle_filter {
  
 public:
  struct Particle *particles;
  struct Particle *oldParticles;
  float *weights;
  float *cumulativeWeightsIndex;
  vector <ZPoint> *measurments;
  Particle_filter(float,float,float,float);
  ~Particle_filter();
  void update();
  void updateMeasurments(vector< ZPoint > *);
  Particle *getParticles();
  virtual float likelihood(Particle &);
  virtual void update_Xt(Particle &);
  float maxX;
  float maxY;
  float minX;
  float minY;

 protected:
  void bound_particle(Particle &);
  virtual void jiggle_particle(Particle &);
  float setup_importance_sample();
  void random_resample(float);
  int binarySearch(float[], int, int, float);
};

class Person_filter : public Particle_filter {
 
 public:
  Person_filter(float,float,float,float);
  float likelihood(Particle &);
  void update_Xt(Particle &);
  void jiggle_particle(Particle &);
};

class Bike_filter : public Particle_filter {
 
 public:
  Bike_filter(float,float,float,float);
  float likelihood(Particle &);
  void update_Xt(Particle &);
  void jiggle_particle(Particle &);
};

#endif

#ifndef PARTICLE_FILTER_H
#define PARTICLE_FILTER_H

#include <vector>
#include "geomModel.h"

#define NUM_PARTICLES 200

typedef struct Particle {
  float x;
  float y;
  float theta;
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
  float maxX;
  float maxY;
  float minX;
  float minY;

private:
  void bound_particle(Particle &);
  void jiggle_particle(Particle &);
  void update_Xt(Particle &);
  float setup_importance_sample();
  float likelihood(Particle &);
  int binarySearch(float[], int, int, float);
};

#endif

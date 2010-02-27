#ifndef PARTICLE_FILTER_H
#define PARTICLE_FILTER_H

#include <vector>
#include "geomModel.h"

#define MAX_X 200
#define MAX_Y 200
#define MIN_X -200
#define MIN_Y -200
#define NUM_PARTICLES 100

typedef struct Particle {
  float x;
  float y;
} Particle;

class Particle_filter {
  
public:

  struct Particle *particles;
  struct Particle *oldParticles;
  float *weights;
  float *cumulativeWeightsIndex;
  vector <ZPoint> *measurments;
  Particle_filter();
  ~Particle_filter();
  void update();
  void updateMeasurments(vector< ZPoint > *);
  Particle *getParticles();

private:
  void bound_particle(Particle &);
  void jiggle_particle(Particle &);
  void update_Xt(Particle &);
  float setup_importance_sample();
  float likelihood(Particle &);
  int binarySearch(float[], int, int, float);
};

#endif

#ifndef PARTICLE_FILTER_H
#define PARTICLE_FILTER_H

#define MAX_X 200
#define MAX_Y 200
#define MIN_X -200
#define MIN_Y -200
#define NUM_PARTICLES 500

using namespace std;

class Particle_filter {
  
public:
  struct particle {
    float x;
    float y;
  };

  struct particle *particles;
  struct particle *oldParticles;
  float *weights;
  float *cumulativeWeightsIndex;
  Particle_filter();
  ~Particle_filter();
  void update();

private:
  void bound_particle(struct particle &);
  void jiggle_particle(struct particle &);
  void update_Xt(struct particle&);
  float setup_importance_sample();
  float likelihood(float, struct particle);
  int binarySearch(float[], int, int, float);
};



#endif

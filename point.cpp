#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <GL/glut.h>

#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <set>
#include <map>
#include <utility> 

#include "ANN/ANN.h"

#include "point.h"
#include "particle-filter.h"
#include "geomModel.h"

#define CONNECTED_COMPONENT_RADIUS 0.3
#define MIN_CONNECTED_COMPONENT_SIZE 10

using namespace std;


// data
vector < vector < ZPoint > > frames;
int curFrame;

// keeps track of view info
float ctrX = 0.0, ctrY = 0.0, ctrZ = 0.0;
float upX = 0.0, upY = 0.0, upZ = 1.0;
float viewRho   = 50.0;	// eye distance from center
float viewTheta = 80*PI/180.0;	// from Z axis
float viewPhi   = 45*PI/180.0;	// from X axis
float markerX = 0.0, markerY = 0.0, markerZ = 0.0;
float meanX, meanY, meanZ;
float maxX, minX, maxY, minY, minZ;

//display options
bool dispRawPoints = false;
bool dispParticlesTracked = true;
bool dispParticlesPotential = false;
bool dispBlobs = true;

// list of connected components
vector< vector<ZPoint> > connectedComponents;
vector< ZPoint > connectedComponentsMean;

// person data structures
vector< Person_filter * > potentialPerson;
vector< Particle > potentialPersonMean;
vector< int > potentialPersonExistenceTime;
vector< Person_filter * > trackedPerson;
vector< Particle > trackedPersonMean;

//bike data structures
vector< Bike_filter * > potentialBike;
vector< Particle > potentialBikeMean;
vector< int > potentialBikeExistenceTime;
vector< Bike_filter * > trackedBike;
vector< Particle > trackedBikeMean;

// different visualization modes
enum ViewMode {VIEW_NORMAL, VIEW_MOVE_MARKER};
ViewMode mode = VIEW_NORMAL;


//find the mean (x,y) of a particle a filter
void findParticleFilterMean(Particle &mean, Particle_filter *filter) {
  mean.x = 0;
  mean.y = 0;
  mean.theta = 0;
  Particle* particles = filter->getParticles();
  for (int i=0; i<NUM_PARTICLES; i++) {
    mean.x += particles[i].x;
    mean.y += particles[i].y;
    mean.theta += particles[i].theta;
  }
  mean.x = mean.x / NUM_PARTICLES;
  mean.y = mean.y / NUM_PARTICLES;
  mean.theta = mean.theta / NUM_PARTICLES;
}

void findConnectedComponentMean(ZPoint &mean, vector<ZPoint> &component) {
  mean.x = 0;
  mean.y = 0;
  mean.z = 0;
  int numPoints = component.size();
  for (int i=0; i<numPoints; i++) {
    mean.x += component[i].x;
    mean.y += component[i].y;
    mean.z += component[i].z;
  }
  mean.x = mean.x / numPoints;
  mean.y = mean.y / numPoints;
  mean.z = mean.z / numPoints;  
}

float computeParticleFilterMeanLikelihood(Particle_filter *filter) {
  double sum = 0;
  Particle *particles = filter->getParticles();
  for (int i=0; i<NUM_PARTICLES; i++) {
    sum += particles[i].likelihood;
  }
  return sum / NUM_PARTICLES;
}

float distSq(ZPoint &point, Particle &particle) {
  float xval = point.x - particle.x;
  float yval = point.y - particle.y;
  return xval*xval+yval*yval;
}


// Create Kd-tree, remove isolated points and return lists of
// connected components.
void findConnectedComponents(vector< ZPoint > &frame)
{
  //clear the previous list of connected componenets from the previous frame
  connectedComponents.clear();
  connectedComponentsMean.clear();

  //load the points in the frame into the Kd-tree
  int numberOfNearestNeighbors = 40;
  int                  nPts = frame.size();// actual number of data points
  ANNpointArray        dataPts;        // data points
  ANNpoint             queryPt;        // query point
  ANNidxArray          nnIdx;          // near neighbor indices
  ANNdistArray         dists;          // near neighbor distances
  ANNkd_tree*          kdTree;         // search structure
  double eps = 0.0001;
  queryPt  = annAllocPt(3);        // allocate query point
  dataPts  = annAllocPts(nPts, 3); // allocate data points
  // allocate near neigh indices  
  nnIdx =  new ANNidx[numberOfNearestNeighbors];
  // allocate near neighbor dists
  dists =  new ANNdist[numberOfNearestNeighbors];

  for (int i=0; i<frame.size(); i++) {
    double curPoint[3] = {frame[i].x, frame[i].y, frame[i].z};
    dataPts[i][0] = curPoint[0];
    dataPts[i][1] = curPoint[1];
    dataPts[i][2] = curPoint[2];
  }
  
  // construct the kd-tree
  kdTree = new ANNkd_tree( // build search structure
			  dataPts, // the data points
			  nPts,    // number of points
			  3);    // dimension of space
  
  //find all connected components
  
  // set of points (their index) already assigned to a conneted components
  set<int> accountedForPoints;
  
  for (int i=0; i<frame.size(); i++) {
    //test if the point has already been added to a connected component
    if (accountedForPoints.count(i) > 0) {
      continue;
    }
    accountedForPoints.insert(i);
    
    // create a new connected component using curPoint as the seed
    vector< ZPoint > connectedComponent;
    set<int> pointsInConnectedComponent;
    vector<int> pointsToSearchFrom;
    pointsInConnectedComponent.insert(i);
    pointsToSearchFrom.push_back(i);
    while (pointsToSearchFrom.size() > 0) {
      int curPtIdx = pointsToSearchFrom.back();
      pointsToSearchFrom.pop_back();
      connectedComponent.push_back(frame[curPtIdx]);
      queryPt[0] = frame[curPtIdx].x;
      queryPt[1] = frame[curPtIdx].y;
      queryPt[2] = frame[curPtIdx].z;
      kdTree->annkSearch( // search
			 queryPt,    // query point
			 numberOfNearestNeighbors,// number of near neighbors
			 nnIdx,      // nearest neighbors (returned)
			 dists,      // distance (returned)
			 eps);       // error bound
      
      // use all points that are with in  CONNECTED_COMPONENT_RADIUS of query point
      for (int j=0; j < numberOfNearestNeighbors; j++) {
	if (dists[j] < CONNECTED_COMPONENT_RADIUS * CONNECTED_COMPONENT_RADIUS) {
	  int index = nnIdx[j];
	  // add this point the connected component
	  if (pointsInConnectedComponent.count(index) == 0) {
	    pointsInConnectedComponent.insert(index);
	    pointsToSearchFrom.push_back(index);
	    accountedForPoints.insert(index);
	    connectedComponent.push_back(frame[index]);
	  }
	} 
      }
    }
    // Only add the connectedComponent if it contains more than MIN points
    if (connectedComponent.size() > MIN_CONNECTED_COMPONENT_SIZE) {
      connectedComponents.push_back(connectedComponent);
    }
  }

  // clean up the kd-tree
  annDeallocPt(queryPt);
  annDeallocPts(dataPts);
  delete [] nnIdx;
  delete [] dists;
  delete kdTree;
  annClose();

  //compute the means for each connectedComponents
  for (int i=0; i<connectedComponents.size(); i++) {
    ZPoint mean;
    findConnectedComponentMean(mean, connectedComponents[i]);
    connectedComponentsMean.push_back(mean);
  }
}

void identifyNewConnectedComponentsToClassifyPerson(vector< Particle > &potentialMean,
						    vector< Particle > &trackedMean,
						    vector<Person_filter*> &potential,
						    vector<Person_filter*> &tracked,
						    vector<int> &existenceTime,
						    string debugName) {
  // for each connected componet find the nearest particle filter
  // if the nearest is too far away, create a new particle filter
  // centered on the connected component
  potentialMean.clear();
  trackedMean.clear();
  for (int j=0; j<potential.size(); j++) {
    Particle mean;
    findParticleFilterMean(mean, potential[j]);
    potentialMean.push_back(mean);
  }
  for (int j=0; j<tracked.size(); j++) {
    Particle mean;
    findParticleFilterMean(mean, tracked[j]);
    trackedMean.push_back(mean);
  }
  
  for (int i=0; i < connectedComponentsMean.size(); i++) {
    ZPoint blobMean = connectedComponentsMean[i];
    
    cout << "blob " << i << " has mean: (" << blobMean.x
	 << ", " << blobMean.y << ")" << endl;

    //find the nearest particle filter mean
    bool alreadyTracked = false;
    for (int j=0; j<potentialMean.size(); j++) {
      Particle particle = potentialMean[j];
      float d = distSq(blobMean, particle);
      if (d < 1) {
	alreadyTracked = true;
	break;
      }
    }

    if (!alreadyTracked) {
      for (int j=0; j<trackedMean.size(); j++) {
	Particle particle = trackedMean[j];
	float d = distSq(blobMean, particle);
	if (d < 1) {
	  alreadyTracked = true;
	  break;
	}
      }
    }

    // The connected component does not have a particle filter assigned to
    // it. Create a new particle filter and center it around the blob.
    if (!alreadyTracked) {
      Person_filter *filter = new Person_filter(blobMean.x+0.6,
						blobMean.x-0.6,
						blobMean.y+0.6,
						blobMean.y-0.6);
      potential.push_back(filter);
      existenceTime.push_back(1);

      cout << "created new " << debugName << " particle filter centered at: (" << blobMean.x
	   << ", " << blobMean.y << ")" << endl;
    }
  }
}

void identifyNewConnectedComponentsToClassifyBike(vector< Particle > &potentialMean,
						  vector< Particle > &trackedMean,
						  vector<Bike_filter*> &potential,
						  vector<Bike_filter*> &tracked,
						  vector<int> &existenceTime,
						  string debugName) {
  // for each connected componet find the nearest particle filter
  // if the nearest is too far away, create a new particle filter
  // centered on the connected component
  potentialMean.clear();
  trackedMean.clear();
  for (int j=0; j<potential.size(); j++) {
    Particle mean;
    findParticleFilterMean(mean, potential[j]);
    potentialMean.push_back(mean);
  }
  for (int j=0; j<tracked.size(); j++) {
    Particle mean;
    findParticleFilterMean(mean, tracked[j]);
    trackedMean.push_back(mean);
  }
  
  for (int i=0; i < connectedComponentsMean.size(); i++) {
    ZPoint blobMean = connectedComponentsMean[i];
    
    cout << "blob " << i << " has mean: (" << blobMean.x
	 << ", " << blobMean.y << ")" << endl;

    //find the nearest particle filter mean
    bool alreadyTracked = false;
    for (int j=0; j<potentialMean.size(); j++) {
      Particle particle = potentialMean[j];
      float d = distSq(blobMean, particle);
      if (d < 1) {
	alreadyTracked = true;
	break;
      }
    }

    if (!alreadyTracked) {
      for (int j=0; j<trackedMean.size(); j++) {
	Particle particle = trackedMean[j];
	float d = distSq(blobMean, particle);
	if (d < 1) {
	  alreadyTracked = true;
	  break;
	}
      }
    }

    // The connected component does not have a particle filter assigned to
    // it. Create a new particle filter and center it around the blob.
    if (!alreadyTracked) {
      Bike_filter *filter = new Bike_filter(blobMean.x+0.6,
					    blobMean.x-0.6,
					    blobMean.y+0.6,
					    blobMean.y-0.6);	
      potential.push_back(filter);
      existenceTime.push_back(1);

      cout << "created new " << debugName << " particle filter centered at: (" << blobMean.x
	   << ", " << blobMean.y << ")" << endl;
    }
  }
}



//remove particle filters whose average likelihood is too low
void removeStalePersonParticleFilters(vector< Person_filter * > &potential,
				vector< Person_filter * > &tracked,
				vector<int> &existenceTime,
				string debugName)
{
  for (int i=0; i<potential.size();) {
    float likelihood = computeParticleFilterMeanLikelihood(potential[i]);    
    cout << "mean likelihood for filter " << debugName << ": " << i
	 << " is " << likelihood << endl;
    if (existenceTime[i] > 3) {
      //This particle is no longer tracking a person, remove it
      if (likelihood < 50) {
	delete potential[i];
      } else {
	// This particle is being graduated to the tracked person list
	tracked.push_back(potential[i]);
      }

      potential.erase(potential.begin()+i);
      existenceTime.erase(existenceTime.begin() + i);
    } else {
      i++;
    }
  }

  for (int i=0; i<tracked.size();) {
    float likelihood = computeParticleFilterMeanLikelihood(tracked[i]);
    cout << "mean likelihood for " << debugName << " filter: " << i
	 << " is " << likelihood << endl;    
    if (likelihood < 20) {
      //delete trackedPerson[i];
      potential.push_back(tracked[i]);
      existenceTime.push_back(3);
      tracked.erase(tracked.begin()+i);
    } else {
      i++;
    }
  }
}

//remove particle filters whose average likelihood is too low
void removeStaleBikeParticleFilters(vector< Bike_filter * > &potential,
				vector< Bike_filter * > &tracked,
				vector<int> &existenceTime,
				string debugName)
{
  for (int i=0; i<potential.size();) {
    float likelihood = computeParticleFilterMeanLikelihood(potential[i]);    
    cout << "mean likelihood for filter " << debugName << ": " << i
	 << " is " << likelihood << endl;
    if (existenceTime[i] > 3) {
      //This particle is no longer tracking a person, remove it
      if (likelihood < 50) {
	delete potential[i];
      } else {
	// This particle is being graduated to the tracked person list
	tracked.push_back(potential[i]);
      }

      potential.erase(potential.begin()+i);
      existenceTime.erase(existenceTime.begin() + i);
    } else {
      i++;
    }
  }

  for (int i=0; i<tracked.size();) {
    float likelihood = computeParticleFilterMeanLikelihood(tracked[i]);
    cout << "mean likelihood for " << debugName << " filter: " << i
	 << " is " << likelihood << endl;    
    if (likelihood < 20) {
      //delete trackedPerson[i];
      potential.push_back(tracked[i]);
      existenceTime.push_back(3);
      tracked.erase(tracked.begin()+i);
    } else {
      i++;
    }
  }
}


void resize(int w, int h)
{
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glViewport(0, 0, w, h);
	// may want to decrease zNear for better visual performance
	gluPerspective(45.0f, 1.0f * w / h, 0.1f, 100000.0f);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

// draws a string in 3d space
// from http://www.lighthouse3d.com/opengl/glut/index.php?bmpfont
void drawString(float x, float y, float z,  char *str) 
{
	//glPushMatrix();
	glRasterPos3f(x,y,z);
	for(char *c = str; *c != '\0'; c++)
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10, *c);
	//glPopMatrix();
}

void drawCoordinateAxes(void)
{
	// define the range over which to draw the ground plane
	int const PLANE_RANGE = 200.0;
	int const LABEL_DELTA = 10.0;
	
	// draw ground plane
	glBegin(GL_POLYGON);
	glColor4f(1,1,1,0.05);	// transparent
	glVertex3f(-PLANE_RANGE,-PLANE_RANGE, 0);
	glVertex3f(-PLANE_RANGE, PLANE_RANGE, 0);
	glVertex3f(PLANE_RANGE, PLANE_RANGE, 0);
	glVertex3f(PLANE_RANGE, -PLANE_RANGE, 0);
	glEnd();
	
	// draw coordinates
	glBegin(GL_LINES);
	glColor4f(1, 1, 1, 0.1);
	for(int x = -PLANE_RANGE; x <= PLANE_RANGE; x+= LABEL_DELTA) {	
		glVertex3f((float)x, -PLANE_RANGE, 0);
		glVertex3f((float)x, PLANE_RANGE,  0);
	}
	for(int y = -PLANE_RANGE; y <= PLANE_RANGE; y += LABEL_DELTA) {	
		glVertex3f(-PLANE_RANGE, (float)y, 0);
		glVertex3f( PLANE_RANGE, (float)y, 0);
	}
	glEnd();
	
	// label coordinate axes
	char buf[10];
	for(int x = -PLANE_RANGE; x <= PLANE_RANGE; x += LABEL_DELTA) {
		sprintf(buf, "%d", x);
		drawString(x, 0.0, 0.0, buf);
	}
	for(int y = -PLANE_RANGE; y <= PLANE_RANGE; y += LABEL_DELTA) {
		sprintf(buf, "%d", y);
		drawString(0.0, y, 0.0, buf);
	}
	
	// draw origin
	//glColor3f(0.75, 0.1, 0.1);
	//glPointSize(10.0);
	//glBegin(GL_POINTS);
	//glVertex3f(0.0, 0.0, 0.0);
	//glEnd();
	
	// draw origin markers
	glBegin(GL_LINES);
	glColor3f(1.0, 0.0, 0.0);
	glVertex3f(0.0, 0.0, 0.0);	// x axis (red)
	glVertex3f(1.0, 0.0, 0.0);
	glColor3f(0.0, 1.0, 0.0);
	glVertex3f(0.0, 0.0, 0.0);	// y axis (green)
	glVertex3f(0.0, 1.0, 0.0);
	glColor3f(0.0, 0.0, 1.0);
	glVertex3f(0.0, 0.0, 0.0);	// z axis (blue)
	glVertex3f(0.0, 0.0, 1.0);
	glEnd();
}

void drawMarker(float x, float y, float z)
{
  glPointSize(5.0);
  glColor4f(0.0, 1.0, 0.0, 0.5);
  glBegin(GL_POINTS);
  glVertex3f(x, y, z);
  glEnd();
}

void draw(void)
{
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  // enable transparency
  // from http://www.opengl.org/resources/faq/technical/transparency.htm
  glEnable (GL_BLEND); 
  glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  
  // calculate where we're looking
  glLoadIdentity();
  gluLookAt(ctrX + viewRho*sin(viewTheta)*cos(viewPhi), ctrY + viewRho*sin(viewTheta)*sin(viewPhi), ctrZ + viewRho*cos(viewTheta),
	    ctrX, ctrY, ctrZ,
	    upX, upY, upZ);
  
  drawCoordinateAxes();
  drawMarker(markerX, markerY, markerZ);
  
  // The raw points from the frame
  if (dispRawPoints) {
    glColor3f(1,1,1);
    glPointSize(1.0);
    glBegin(GL_POINTS);
    for(int i = 0; i < frames[curFrame].size(); i++) {
      glVertex3f(frames[curFrame][i].x,
		 frames[curFrame][i].y,
		 frames[curFrame][i].z);
    }
    glEnd();
  }
  
  //The blobs found as connected components
  if (dispBlobs) {
    int numComponents = connectedComponents.size();
    for (int i=0; i < numComponents; i++) {
      vector< ZPoint > component = connectedComponents[i];
      double colorVal = ((double) i)/numComponents;
      if (i % 3 == 0)
	glColor3f(colorVal, 0, 0);
      else if (i % 3 == 1)
	glColor3f(0, colorVal, 0);
      else if (i % 3 == 2)
	glColor3f(0, 0, colorVal);
      glPointSize(2.0);
      glBegin(GL_POINTS);
      for (int j=0; j < component.size(); j++) {
	glVertex3f(component[j].x,
		   component[j].y,
		   component[j].z);
	    }
      glEnd();
    }
  }

  //particle filter points
  if (dispParticlesPotential) {
    glColor3f(1,0,0);
    glPointSize(3);
    glBegin(GL_POINTS);
    for (int n=0; n<potentialPerson.size(); n++) {
      Particle *particles = potentialPerson[n]->getParticles();
      for (int i=0; i<NUM_PARTICLES; i++) {
	Particle p = particles[i];
	glVertex3f(p.x, p.y, minZ);
      }
    }
    glEnd();
  }
  
  if (dispParticlesTracked) {
    glColor3f(0,1,0);
    glPointSize(3);
    glBegin(GL_POINTS);
    for (int n=0; n<trackedPerson.size(); n++) {
      Particle *particles = trackedPerson[n]->getParticles();
      for (int i=0; i<NUM_PARTICLES; i++) {
	Particle p = particles[i];
	glVertex3f(p.x, p.y, minZ);
      }
    }
    glEnd();
  }
  glFlush();
  glutPostRedisplay();
  glutSwapBuffers();	
}

void nextFrame()
{
	
	
}

void specialKey(int k, int x, int y)
{
  switch (k)
    {
    case GLUT_KEY_UP: // up arrow
      viewRho -= 0.05*viewRho;
      break;
    case GLUT_KEY_DOWN: // down arrow
      viewRho += 0.05*viewRho;
      break;
    case GLUT_KEY_RIGHT: // right arrow
      viewPhi += 3.0*PI/180.0;
      break;
    case GLUT_KEY_LEFT: // left arrow
      viewPhi -= 3.0*PI/180.0;
      break;
    }
}

void updatePersonParticleFilters(vector<Person_filter*> &potential,
				 vector<Person_filter*> &tracked,
				 vector<int> &existenceTime)
{
  for (int i=0; i<tracked.size(); i++) {
	tracked[i]->updateMeasurments(&connectedComponents);
	tracked[i]->setConnectedComponentMeans(&connectedComponentsMean);
	tracked[i]->update();
      }
  for (int i=0; i<potential.size(); i++) {
    potential[i]->updateMeasurments(&connectedComponents);
    potential[i]->setConnectedComponentMeans(&connectedComponentsMean);
    potential[i]->update();
    existenceTime[i] += 1;
  }
}

void updateBikeParticleFilters(vector<Bike_filter*> &potential,
			       vector<Bike_filter*> &tracked,
			       vector<int> &existenceTime)
{
  for (int i=0; i<tracked.size(); i++) {
    tracked[i]->updateMeasurments(&connectedComponents);
    tracked[i]->setConnectedComponentMeans(&connectedComponentsMean);
    tracked[i]->update();
  }
  for (int i=0; i<potential.size(); i++) {
    potential[i]->updateMeasurments(&connectedComponents);
    potential[i]->setConnectedComponentMeans(&connectedComponentsMean);
    potential[i]->update();
    existenceTime[i] += 1;
  }
}

void key(unsigned char k, int x, int y)
{
  switch (k)
    {
    case 'q':
      exit(0); 
      break;
    case 'p':
      dispRawPoints =  !dispRawPoints;
      break;
    case 'o':
      dispBlobs = !dispBlobs;
      break;
    case 'i':
      dispParticlesPotential = !dispParticlesPotential;
      break;
    case 'u':
      dispParticlesTracked = !dispParticlesTracked;
      break;
    case 'd':
      ctrX -= 0.05*viewRho*cos(viewTheta)*sin(viewPhi);
      ctrY += 0.05*viewRho*cos(viewTheta)*cos(viewPhi);
      break;
    case 'a':
      ctrX += 0.05*viewRho*cos(viewTheta)*sin(viewPhi);
      ctrY -= 0.05*viewRho*cos(viewTheta)*cos(viewPhi);
      break;
    case 'w':
      ctrX -= 0.05*viewRho*sin(viewTheta)*cos(viewPhi);
      ctrY -= 0.05*viewRho*sin(viewTheta)*sin(viewPhi);
      //ctrZ -= 0.05*viewRho*cos(viewTheta);
      break;
    case 's':
      ctrX += 0.05*viewRho*sin(viewTheta)*cos(viewPhi);
      ctrY += 0.05*viewRho*sin(viewTheta)*sin(viewPhi);
      //ctrZ += 0.05*viewRho*cos(viewTheta);
      break;
    case 'e':
      viewTheta -= 3.0*PI/180.0;
      if(viewTheta < 0.1)	// prevent from discontinuous viewing
	viewTheta = 0.01;
      break;
    case 'c':
      viewTheta += 3.0*PI/180.0;
      if(viewTheta > 179.9*PI/180.0)
	viewTheta = 179.99*PI/180.0;
      break;
    case '1':
      if (curFrame < frames.size() - 1) curFrame++;			
      printf("Current frame: %d\n", frames[curFrame][0].t);
      findConnectedComponents(frames[curFrame]);
      cout << "there are " << connectedComponents.size() << " blobs" << endl; 
      cout << "there are " << trackedPerson.size() << " tracked persons" << endl;
      cout << "there are " << potentialPerson.size() << " potential persons" << endl;
      updatePersonParticleFilters(potentialPerson,
			    trackedPerson,
			    potentialPersonExistenceTime);
      updateBikeParticleFilters(potentialBike,
			    trackedBike,
			    potentialBikeExistenceTime);      
      removeStalePersonParticleFilters(potentialPerson,
				       trackedPerson,
				       potentialPersonExistenceTime,
				       "person");
      removeStaleBikeParticleFilters(potentialBike,
				     trackedBike,
				     potentialBikeExistenceTime,
				     "bike");
      identifyNewConnectedComponentsToClassifyPerson(potentialPersonMean,
						     trackedPersonMean,
						     potentialPerson,
						     trackedPerson,
						     potentialPersonExistenceTime,
						     "person");
      identifyNewConnectedComponentsToClassifyBike(potentialBikeMean,
						   trackedBikeMean,
						   potentialBike,
						   trackedBike,
						   potentialBikeExistenceTime,
						   "bike");
      break;
    case '!':
      if(curFrame > 0) curFrame--;
      printf("Current frame: %d\n", frames[curFrame][0].t);
      break;
    }
}

void readDataFile(vector < vector<ZPoint> > &frames, char *filename)
{
	ifstream fin(filename);
	string line;
	float x,y,z;
	int t;
	
	printf("Reading file %s...", filename);
	int lastT = -1;
	while(getline(fin, line))
	{
		if(line[0] == '#') continue;
		sscanf(line.c_str(), "%f,%f,%f,%d",&x,&y,&z,&t);
		ZPoint pt;
		pt.x = x;
		pt.y = y;
		pt.z = z;
		pt.t = t;
		
		if(t > lastT)
		{
			vector <ZPoint> newFrame;
			frames.push_back(newFrame);
			lastT = t;
		}
		
		frames.back().push_back(pt);
	}
	printf("done!\n");
	fin.close();
}

void preProcessData(vector< vector<ZPoint> > &frames)
{
  meanX = 0.0;
  meanY = 0.0;
  meanZ = 0.0;
  int totalPoints = 0;
  
  printf("-- Preprocessing data --\n");
  printf("Subtracting off the mean..."); fflush(stdout);
  // center the data set so that the mean is at 0,0,0
  // calculate mean
  for(int fr = 0; fr < frames.size(); fr++) {
    for(int i = 0; i < frames[fr].size(); i++) {
      meanX += frames[fr][i].x;
      meanY += frames[fr][i].y;
      meanZ += frames[fr][i].z;
      totalPoints++;
    }
  }
  meanX /= (float)totalPoints;
  meanY /= (float)totalPoints;
  meanZ /= (float)totalPoints;
  
  // subtract mean (and smallest value of z)
  for(int fr = 0; fr < frames.size(); fr++) {
    for(int i = 0; i < frames[fr].size(); i++) {
      frames[fr][i].x -= meanX;
      frames[fr][i].y -= meanY;
      // also subtract away the z-mean (ground plane will cut through dataset)
      frames[fr][i].z -= meanZ;
    }
  }
  printf("done!\n");
  
  printf("Finding an interesting square...");
  fflush(stdout);
  vector< vector<ZPoint> > newFrames;
  for(int fr = 0; fr < frames.size(); fr++)
    {
      vector< ZPoint > thisFrame;
      newFrames.push_back(thisFrame);
      for(int i = 0; i < frames[fr].size(); i++)
  	{
	  if(-25.0 <= frames[fr][i].x && frames[fr][i].x <= 0.0 &&
	     -25.0 <= frames[fr][i].y && frames[fr][i].y <= 0.0)
	    {
	      newFrames.back().push_back(frames[fr][i]);
	    }
  	}
    }
  printf("done!\n");
  
  maxX = 0;
  minX = -25;
  maxY = 0;
  minY = -25;
  minZ = 0;
  
  frames = newFrames;
  
  printf("Moving points up..."); 
  float minZFound = 1000.0;
  fflush(stdout);
  for(int fr = 0; fr < frames.size(); fr++) {
  	for(int i = 0; i < frames[fr].size(); i++) {
  		if(frames[fr][i].z < minZFound)
  			minZFound = frames[fr][i].z;
  	}
  }
  // hack
  minZFound = -1.0;
  for(int fr = 0; fr < frames.size(); fr++) {
    for(int i = 0; i < frames[fr].size(); i++) {
      frames[fr][i].z -= minZFound;
    }
  }
  printf("done!");
}


int main(int argc, char **argv)
{
  // read the points
  if(argc < 2) 
    {
      cout << argv[0] << " <point_cloud_file.txt>" << endl;
      exit(0);
    }
  
  // read and display data
  readDataFile(frames, argv[1]);
  preProcessData(frames); 	// subtracts away the mean
  
  // initialize particle filter
  //filter = new Particle_filter(maxX, minX, maxY, minY);
  
  curFrame = 0;
  findConnectedComponents(frames[curFrame]);
  printf("Current frame: %d\n", frames[curFrame][0].t);
  
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
  glutInitWindowSize(600, 500);
  
  glutCreateWindow("Points");
  
  glutSpecialFunc(specialKey);
  glutKeyboardFunc(key);
  //glutVisibilityFunc(visible);
  glutDisplayFunc(draw);
  glutReshapeFunc(resize);
  
  glutMainLoop();
  return 0;
}

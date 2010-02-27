#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <GL/glut.h>

#include <vector>
#include <fstream>
#include <iostream>
#include <string>

#include "point.h"
#include "particle-filter.h"
#include "geomModel.h"

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
float minZ;

// different visualization modes
enum ViewMode {VIEW_NORMAL, VIEW_MOVE_MARKER};
ViewMode mode = VIEW_NORMAL;

Particle_filter *filter;

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
	int const PLANE_RANGE = 1000.0;
	int const LABEL_DELTA = 100.0;
	
	// draw ground plane
	glBegin(GL_POLYGON);
	glColor4f(1,1,1,0.05);
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
	
	// test points
	glColor3f(1,1,1);
	glPointSize(1.0);
	glBegin(GL_POINTS);
	// test data
	//for(int i = 0; i < 10; i++)
	//	for(int j = 0; j < 10; j++)
	//		for(int k = 0; k < 10; k++)
	//			glVertex3f(i,j,k);	

	// point data
	for(int i = 0; i < frames[curFrame].size(); i++)
	{
		glVertex3f(frames[curFrame][i].x, frames[curFrame][i].y, frames[curFrame][i].z);
		//cout << "point: " << frames[curFrame][i].x << ", " << frames[curFrame][i].y << ", " <<  frames[curFrame][i].z << endl;
	}
	glEnd();

	//particle filter points
	glColor3f(1,0,0);
	glPointSize(5.0);
	glBegin(GL_POINTS);
	Particle *particles = filter->getParticles();
	for (int i=0; i<NUM_PARTICLES; i++) {
	  Particle p = particles[i];
	  glVertex3f(p.x, p.y, minZ);
	  cout << "particle filter: " << p.x << ", " << p.y << endl;
	}
	glEnd();
	
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

void key(unsigned char k, int x, int y)
{
	switch (k)
	{
		case 'q':
			exit(0); 
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
			if(curFrame < frames.size() - 1) curFrame++;
			printf("Current frame: %d\n", frames[curFrame][0].t);
			filter->updateMeasurments(&frames[curFrame]);
			filter->update();
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
  minZ = 1000;
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
			if (frames[fr][i].z < minZ) {
			  minZ = frames[fr][i].z;
			}
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
			//frames[fr][i].z -= meanZ;
		}
	}
	printf("done!\n");
	
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
	curFrame = 0;
	printf("Current frame: %d\n", frames[curFrame][0].t);
	
	// initialize particle filter
	filter = new Particle_filter();

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

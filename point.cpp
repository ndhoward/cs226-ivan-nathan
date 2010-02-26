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

#define PI 3.14159265

using namespace std;

GLfloat viewDist = 2.3;
GLfloat viewAngle = -5.93;
GLfloat eyeXdelta = -75.0;
GLfloat eyeYdelta = 3.16;
GLfloat eyeZ = 1.4;
GLfloat centerX = -75.0;
GLfloat centerY = 3.16;
GLfloat centerZ = 1.4;


vector < vector < struct point > > frames;
int curFrame;

// keeps track of view info
float ctrX = 0.0, ctrY = 0.0, ctrZ = 0.0;
float upX = 0.0, upY = 0.0, upZ = 1.0;
float viewRho   = 1.0;	// eye distance from center
float viewTheta = 45*PI/180.0;	// from Z axis
float viewPhi   = 45*PI/180.0;	// from X axis

Particle_filter *filter;

void resize(int w, int h)
{
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glViewport(0, 0, w, h);
	gluPerspective(45.0f, 1.0f * w / h, 1.0f, 100.0f);
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
	// draw coordinates
	glBegin(GL_LINES);
	glColor4f(1, 1, 1, 0.1);
	for(int x = -10; x <= 10; x++) {	
		glVertex3f((float)x, -10.0, 0);
		glVertex3f((float)x, 10.0,  0);
	}
	for(int y = -10; y <= 10; y++) {	
		glVertex3f(-10.0, (float)y, 0);
		glVertex3f( 10.0, (float)y, 0);
	}
	glEnd();
	
	// label coordinate axes
	char buf[10];
	for(int x = -10; x <= 10; x++) {
		sprintf(buf, "%d", x);
		drawString(x, 0.0, 0.0, buf);
	}
	for(int y = -10; y <= 10; y++) {
		sprintf(buf, "%d", y);
		drawString(0.0, y, 0.0, buf);
	}
	
	// draw origin
	glColor3f(0.75, 0.1, 0.1);
	glPointSize(10.0);
	glBegin(GL_POINTS);
	glVertex3f(0.0, 0.0, 0.0);
	glEnd();
	
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
	
	// test points
	glColor3f(1,1,1);
	glPointSize(2.0);
	glBegin(GL_POINTS);
	//	for(int i = 0; i < 10; i++)
	//	for(int j = 0; j < 10; j++)
	//		for(int k = 0; k < 10; k++)
	//			glVertex3f(i,j,k);	
	//glEnd();
	for(int i = 0; i < frames[curFrame].size(); i++)
	{
		glVertex3f(0.1*frames[curFrame][i].x, 0.1*frames[curFrame][i].y, 0.1*frames[curFrame][i].z);
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
			viewAngle -= 5.0*3.14159/180;
			break;
		case GLUT_KEY_LEFT: // left arrow
			viewAngle += 5.0*3.14159/180;
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
			centerX -= 0.1*sin(viewAngle);
			centerY += 0.1*cos(viewAngle);
			eyeXdelta -= 0.1*sin(viewAngle);
			eyeYdelta += 0.1*cos(viewAngle);
			break;
		case 'a':
			centerX += 0.1*sin(viewAngle);
			centerY -= 0.1*cos(viewAngle);
			eyeXdelta += 0.1*sin(viewAngle);
			eyeYdelta -= 0.1*cos(viewAngle);
			break;
		case 'w':
			centerX -= 0.1*cos(viewAngle);
			centerY -= 0.1*sin(viewAngle);
			eyeXdelta -= 0.1*cos(viewAngle);
			eyeYdelta -= 0.1*sin(viewAngle);
			printf("viewDist=%f\n viewAngle%f\n eyeXdelta=%f\n eyeYdelta=%f\n eyeZ=%f\n centerX=%f\n centerY=%f\n centerZ=%f\n", 
		        viewDist,     viewAngle,    eyeXdelta,     eyeYdelta,     eyeZ,     centerX,     centerY,     centerZ);
			break;
		case 's':
			centerX += 0.1*cos(viewAngle);
			centerY += 0.1*sin(viewAngle);
			eyeXdelta += 0.1*cos(viewAngle);
			eyeYdelta += 0.1*sin(viewAngle);
			break;
		case 'e':
			eyeZ += 0.1;
			centerZ -= 0.1;
			break;
		case 'c':
			eyeZ -= 0.1;
			centerZ += 0.1;
			break;
			
		case '1':
			if(curFrame < frames.size() - 1) curFrame++;
			printf("Current frame: %d\n", frames[curFrame][0].t);
			filter->updateMeasurments(&frames[curFrame]);
			break;
		case '!':
			if(curFrame > 0) curFrame--;
			printf("Current frame: %d\n", frames[curFrame][0].t);
			break;
	}
}

int main(int argc, char **argv)
{
	
	// read the points
	if(argc < 2) 
	{
		cout << argv[0] << " <point_cloud_file.txt>" << endl;
		exit(0);
	}
	
	ifstream fin(argv[1]);
	string line;
	float x,y,z;
	int t;
	
	printf("Reading file %s...", argv[1]);
	int lastT = -1;
	while(getline(fin, line))
	{
		if(line[0] == '#') continue;
		sscanf(line.c_str(), "%f,%f,%f,%d",&x,&y,&z,&t);
		struct point pt;
		pt.x = x;
		pt.y = y;
		pt.z = z;
		pt.t = t;
		
		if(t > lastT)
		{
			vector <struct point> newFrame;
			frames.push_back(newFrame);
			lastT = t;
		}
		
		frames.back().push_back(pt);
	}
	printf("done!\n");
	fin.close();
	
	// initialize particle filter
	filter = new Particle_filter();

	// initialize current frame
	curFrame = 0;
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

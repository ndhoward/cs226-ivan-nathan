#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <GL/glut.h>

#include <vector>
#include <fstream>
#include <iostream>
#include <string>

using namespace std;

GLfloat viewDist = 2.3;
GLfloat viewAngle = -5.93;
GLfloat eyeXdelta = -75.0;
GLfloat eyeYdelta = 3.16;
GLfloat eyeZ = 1.4;
GLfloat centerX = -75.0;
GLfloat centerY = 3.16;
GLfloat centerZ = 1.4;

struct point {
	float x,y,z;
	int t;
};
vector < vector < struct point > > frames;
int curFrame;

void resize(int w, int h)
{
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glViewport(0, 0, w, h);
	gluPerspective(45.0f, 1.0f * w / h, 1.0f, 100.0f);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}


void draw(void)
{
	int i,j,k;
	
	glClear(GL_COLOR_BUFFER_BIT);
	glColor3f(1.0, 1.0, 1.0);
	glLoadIdentity();
	gluLookAt(viewDist*cos(viewAngle) + eyeXdelta,viewDist*sin(viewAngle) + eyeYdelta,eyeZ, 
			  centerX, centerY, centerZ,
			  0.0,0.0,1.0);
	
	glBegin(GL_POINTS);
	/*
	for(i = 0; i < 10; i++)
		for(j = 0; j < 10; j++)
			for(k = 0; k < 10; k++)
				glVertex3f(0.1*i,0.1*j,0.1*k);	
	*/
	
	for(i = 0; i < frames[curFrame].size(); i++)
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
			viewDist -= 0.1;
			break;
		case GLUT_KEY_DOWN: // down arrow
			viewDist += 0.1;
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
	
	// initialize current frame
	curFrame = 0;
	printf("Current frame: %d\n", frames[curFrame][0].t);

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
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

#include <cstdlib>
#include <sys/time.h>

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#endif

#include <glvu.h>

#include "WALL.h"
#include "PARTICLE.h"
#include "PARTICLE_SYSTEM.h"
#include <vector>
#include "CH1.h"


GLVU glvu;

PARTICLE_SYSTEM *particleSystem;

double dt = 1.0 / 100.0;
bool animate = false;

int iterationCount = 0;

double arUtilTimer(void);
void arUtilTimerReset(void);


// Draw coordinate axes
void drawAxes()
{

  glPushMatrix();
  glTranslatef(-0.1f, -0.1f, -0.1f);
  glLineWidth(3.0f);
  glBegin(GL_LINES);
  // x axis is red
  glColor4f(10.0f, 0.0f, 0.0f, 1.0f);
  glVertex3f(0.0f, 0.0f, 0.0f);
  glColor4f(10.0f, 0.0f, 0.0f, 0.0f);
  glVertex3f(1.0f, 0.0f, 0.0f);
  
  // y axis is green 
  glColor4f(0.0f, 10.0f, 0.0f, 1.0f);
  glVertex3f(0.0f, 0.0f, 0.0f);
  glColor4f(0.0f, 10.0f, 0.0f, 0.0f);
  glVertex3f(0.0f, 1.0f, 0.0f);
  
  // z axis is blue
  glColor4f(0.0f, 0.0f, 10.0f, 1.0f);
  glVertex3f(0.0f, 0.0f, 0.0f);
  glColor4f(0.0f, 0.0f, 10.0f, 0.0f);
  glVertex3f(0.0f, 0.0f, 1.0f);
  glEnd();
  glLineWidth(1.0f);
  glPopMatrix();
}



// The drawing function
void displayCallback()
{
  glvu.BeginFrame();
  
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glEnable(GL_DEPTH_TEST);

    
  particleSystem->draw();

  
  glvu.EndFrame();
}


// The projection function
void reshapeCallback(int width, int height)
{
  cout << "reshape" << endl;
  glViewport(0, 0, (GLsizei) width, (GLsizei) height);

  glMatrixMode(GL_PROJECTION);

  glLoadIdentity();
  gluPerspective(65.0, (float)width / height, 0.01, 1000.0);

  // set the matric mode back to modelview
  glMatrixMode(GL_MODELVIEW);

  glLoadIdentity();
  gluLookAt(0, 0, 0.5, 0, 0, 0, 0, 1, 0);
}


// Keyboard commands
void keyboardCallback(unsigned char key, int x, int y)
{
  switch (key)
  {
    // quit entirely
    case 'q':
    case 'Q':
      exit(0);
      break;

    case 'a':
      animate = !animate;
      break;
      
    case 'g':
#ifndef BRUTE
      particleSystem->toggleGridVisble();
#endif
      break;

    case '1':
      particleSystem->loadScenario(SCENARIO_DAM);
      break;
      
    case 'f':
      printf("*** %f (frame/sec)\n", (double)iterationCount/arUtilTimer());
      break;
  }
  glvu.Keyboard(key, x, y);
  
  glutPostRedisplay();
  
}


void specialKeys(int key, int x, int y) {
    float moveStep = 0.01f;
    float minX = -1.0f, maxX = 1.0f; // Define min and max x
    float minY = -1.0f, maxY = 1.0f; // Define min and max y
    float minZ = -1.0f, maxZ = 1.0f; // Define min and max z

    switch (key) {
        case GLUT_KEY_LEFT:
            if (PARTICLE_SYSTEM::rectX > minX)
                PARTICLE_SYSTEM::rectX -= moveStep;
            break;
        case GLUT_KEY_RIGHT:
            if (PARTICLE_SYSTEM::rectX < maxX)
                PARTICLE_SYSTEM::rectX += moveStep;
            break;
        case GLUT_KEY_UP:
            if (PARTICLE_SYSTEM::rectY < maxY)
                PARTICLE_SYSTEM::rectY += moveStep;
            break;
        case GLUT_KEY_DOWN:
            if (PARTICLE_SYSTEM::rectY > minY)
                PARTICLE_SYSTEM::rectY -= moveStep;
            break;
        case GLUT_KEY_PAGE_UP:
            if (PARTICLE_SYSTEM::rectZ < maxZ)
                PARTICLE_SYSTEM::rectZ += moveStep;
            break;
        case GLUT_KEY_PAGE_DOWN:
            if (PARTICLE_SYSTEM::rectZ > minZ)
                PARTICLE_SYSTEM::rectZ -= moveStep;
            break;
    }

    glutPostRedisplay(); // Update the display after a key press
}





void glutMouseClick(int button, int state, int x, int y)
{
  glvu.Mouse(button,state,x,y);
}


void glutMouseMotion(int x, int y)
{
  glvu.Motion(x,y);
  glvuVec3f viewVector = glvu.GetCurrentCam()->Y;
  particleSystem->setGravityVectorWithViewVector(VEC3D(viewVector[0],viewVector[1],viewVector[2]));  
}


// Idle command processing
void idleCallback()
{
  if (!animate) return;

#ifdef BRUTE
  particleSystem->stepVerletBrute(dt);
#else
  if (iterationCount == 0) arUtilTimerReset();
  particleSystem->stepVerlet(dt);
  iterationCount++;
#endif

  glutPostRedisplay();
}

int main(int argc, char** argv)
{
  char title[] = "sph";

  glutInit(&argc, argv);
  glvu.Init(title, GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH, 0, 0, 800, 800);
  glShadeModel(GL_SMOOTH);
  
  glvu.SetInertiaEnabled(0);

  glutDisplayFunc(displayCallback); 
  glutIdleFunc(idleCallback); 
  glutKeyboardFunc(keyboardCallback);
  glutMouseFunc(glutMouseClick);
  glutMotionFunc(glutMouseMotion);
  glutSpecialFunc(specialKeys);

  

  // set background to black
  glClearColor(1.0, 1.0, 1.0, 1.0);

  // enable lights
  GLfloat ambient[] = {0.7,0.7,0.7};
  GLfloat diffuse[] = {1.0,1.0,1.0};
  GLfloat specular[] = {0.0, 0.0, 0.0};

  glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse);
  glLightfv(GL_LIGHT0, GL_SPECULAR, specular);
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  
  glvuVec3f ModelMin(-10,-10,-10), ModelMax(10,10,10), 
  Eye(0, 0, 1.5),  LookAtCntr(0, 0, 0),  Up(0, 1, 0);
  
  float Yfov = 45;
  float Aspect = 1;
  float Near = 0.001f;
  float Far = 10.0f;
  glvu.SetAllCams(ModelMin, ModelMax, Eye, LookAtCntr, Up, Yfov, Aspect, Near, Far);
  
  glvuVec3f center(0.0, 0.0, 0.0);
  glvu.SetWorldCenter(center);
  
  particleSystem = new PARTICLE_SYSTEM();

  glutMainLoop();
  
  return 0;
}

static int      ss, sms;

double arUtilTimer(void)
{
#ifdef _WIN32
  struct _timeb sys_time;
  double             tt;
  int                s1, s2;
  
  _ftime(&sys_time);
  s1 = sys_time.time  - ss;
  s2 = sys_time.millitm - sms;
#else
  struct timeval     time;
  double             tt;
  int                s1, s2;
  
#if defined(__linux) || defined(__APPLE__)
  gettimeofday( &time, NULL );
#else
  gettimeofday( &time );
#endif
  s1 = time.tv_sec  - ss;
  s2 = time.tv_usec/1000 - sms;
#endif
  
  tt = (double)s1 + (double)s2 / 1000.0;
  
  return( tt );
}

void arUtilTimerReset(void)
{
#ifdef _WIN32
  struct _timeb sys_time;
  
  _ftime(&sys_time);
  ss  = sys_time.time;
  sms = sys_time.millitm;
#else
  struct timeval     time;
  
#if defined(__linux) || defined(__APPLE__)
  gettimeofday( &time, NULL );
#else
  gettimeofday( &time );
#endif
  ss  = time.tv_sec;
  sms = time.tv_usec / 1000;
#endif
}


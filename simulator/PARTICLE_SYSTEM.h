#ifndef PARTICLE_SYSTEM_H
#define PARTICLE_SYSTEM_H


#include "PARTICLE.h"
#include "WALL.h"
#include <vector>
//#include <tr1/tuple>
//#include <map>
#include "FIELD_3D.h"

#include "CH1.h"

#define h 0.0457 //0.02 //0.045

#define GAS_STIFFNESS 1.0 
#define REST_DENSITY 1050 // kg/m^3 is rest density of blood particle
#define PARTICLE_MASS 0.02 // kg
#define VISCOSITY 10.0 // Ns/m^2 or Pa*s viscosity of blood
#define SURFACE_TENSION 0.0728 // N/m 
#define SURFACE_THRESHOLD 7.065
#define KERNEL_PARTICLES 20.0

#define GRAVITY_ACCELERATION -9.80665


#define WALL_K 10000.0 // wall spring constant
#define WALL_DAMPING -0.9 // wall damping constant

#define BOX_SIZE 0.4
#define MAX_PARTICLES 3000

#define INITIAL_SCENARIO SCENARIO_DAM

using namespace std;

class PARTICLE_SYSTEM {
  
public:
  static float rectX, rectY, rectZ;
  static float width;
  static float height;
  static float depth;
  static std::vector<VEC3D> attractionPoints1;
  static std::vector<VEC3D> attractionPoints2;  // Position of the attraction point
  static float attractionDistance;  // Effective distance \( d \)
  PARTICLE_SYSTEM();
  ~PARTICLE_SYSTEM();

  void updateGrid();
  
  // draw to OGL
  void draw();
  
  void addParticle(const VEC3D& position);
  
  void addParticle(const VEC3D& position, const VEC3D& velocity);
  
  void stepVerlet(double dt);
  
  void stepVerletBrute(double dt);
    
  void calculateAcceleration();
  
  void calculateAccelerationBrute();
  
  void collisionForce(PARTICLE& particle, VEC3D& f_collision);
  
  void getNeighborParticles(vector<PARTICLE>& totalNeighborParticles, int x, int y, int z);
  
  double Wpoly6(double radiusSquared);
  
  void Wpoly6Gradient(VEC3D& diffPosition, double radiusSquared, VEC3D& gradient);
  
  double Wpoly6Laplacian(double radiusSquared); 
  
  void WspikyGradient(VEC3D& diffPosition, double radiusSquared, VEC3D& gradient);
  
  double WviscosityLaplacian(double radiusSquared);
  
  void toggleGridVisble();
  
  void toggleSurfaceVisible();
  
  void toggleGravity();
  
  void toggleArrows();
  
  void toggleTumble();
  
  void generateFaucetParticleSet();
  
  void setGravityVectorWithViewVector(VEC3D viewVector);

  
  
  void loadScenario(int scenario);
  
  ////Stick
void drawSolidRectangle(GLfloat width, GLfloat height, GLfloat depth, GLubyte R, GLubyte G, GLubyte B, GLubyte alpha) {
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    // Set the color of the rectangle
    glColor4ub(R, G, B, alpha);

    GLfloat halfWidth = width / 2.0f;
    GLfloat halfHeight = height / 2.0f;
    GLfloat halfDepth = depth / 2.0f;

    // Start drawing the rectangle
    glBegin(GL_QUADS);

    glDisable(GL_BLEND);

    // Front face
    glVertex3f(-halfWidth, -halfHeight, halfDepth);
    glVertex3f(halfWidth, -halfHeight, halfDepth);
    glVertex3f(halfWidth, halfHeight, halfDepth);
    glVertex3f(-halfWidth, halfHeight, halfDepth);

    // Back face
    glVertex3f(-halfWidth, -halfHeight, -halfDepth);
    glVertex3f(-halfWidth, halfHeight, -halfDepth);
    glVertex3f(halfWidth, halfHeight, -halfDepth);
    glVertex3f(halfWidth, -halfHeight, -halfDepth);

    // Top face
    glVertex3f(-halfWidth, halfHeight, -halfDepth);
    glVertex3f(-halfWidth, halfHeight, halfDepth);
    glVertex3f(halfWidth, halfHeight, halfDepth);
    glVertex3f(halfWidth, halfHeight, -halfDepth);

    // Bottom face
    glVertex3f(-halfWidth, -halfHeight, -halfDepth);
    glVertex3f(halfWidth, -halfHeight, -halfDepth);
    glVertex3f(halfWidth, -halfHeight, halfDepth);
    glVertex3f(-halfWidth, -halfHeight, halfDepth);

    // Right face
    glVertex3f(halfWidth, -halfHeight, -halfDepth);
    glVertex3f(halfWidth, halfHeight, -halfDepth);
    glVertex3f(halfWidth, halfHeight, halfDepth);
    glVertex3f(halfWidth, -halfHeight, halfDepth);

    // Left face
    glVertex3f(-halfWidth, -halfHeight, -halfDepth);
    glVertex3f(-halfWidth, -halfHeight, halfDepth);
    glVertex3f(-halfWidth, halfHeight, halfDepth);
    glVertex3f(-halfWidth, halfHeight, -halfDepth);

    glEnd();
}



///// End stick
  
  
  
  FIELD_3D* grid;
  FIELD_3D* grid_new;
  double surfaceThreshold;
  VEC3D gravityVector;

private:
  // list of particles, walls, and springs being simulated
  vector<PARTICLE> _particles;
  vector<WALL>     _walls;

  //unsigned int _particleCount;
  bool _isGridVisible;
  bool _tumble;
  
  VEC3D boxSize;

};

#endif

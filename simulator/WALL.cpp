#include "WALL.h"

#define thickness 0.02

// Constructor
WALL::WALL(const VEC3D& normal, const VEC3D& point, bool isLarge) :
  _normal(normal), _point(point), isLargeWall(isLarge)
{
  _normal.normalize();
}

// OGL drawing
void WALL::draw() 
{
  glPushMatrix();
    // translate to the point
    glTranslated(_point[0], _point[1], _point[2]);
    
    // apply a rotation
    double angle1 = asin(_normal[0]) / (2 * M_PI) * 360.0;
  double angle2 = asin(_normal[1]) / (2 * M_PI) * 360.0;

  float width, height, depth;
    if (isLargeWall) {
        width = 20.0;  // Large wall width
        height = 20.0; // Large wall height
        depth = 1;   // Same depth for simplicity
    } else {
        width = 1.0;  // Small wall width
        height = 1.0; // Small wall height
        depth = 1.0;   // Same depth for simplicity
    }

    glRotatef(-angle1, 0, 1, 0);
  glRotatef(-angle2, 1, 0, 0);

    glTranslated(0, 0, thickness/2.0);

  
    
  glPopMatrix();
}



#include "bbox.h"

#include "GL/glew.h"

#include <algorithm>
#include <iostream>

namespace CMU462 {

bool BBox::intersect(const Ray& r, double& t0, double& t1) const {

  // TODO:
  // Implement ray - bounding box intersection test
  // If the ray intersected the bouding box within the range given by
  // t0, t1, update t0 and t1 with the new intersection times.


  Vector3D orig = r.o;
  Vector3D d = r.d;
  double mins[3], maxs[3];
  for (int dim = 0; dim < 3; dim++) {
    double a = (this->min[dim] - orig[dim]) * r.inv_d[dim];
    double b = (this->max[dim] - orig[dim]) * r.inv_d[dim];
    mins[dim] = std::min(a,b);
    maxs[dim] = std::max(a,b);
  }
 
  double min_t = std::max(mins[0], std::max(mins[1], mins[2]));
  double max_t = std::min(maxs[0], std::min(maxs[1], maxs[2]));

  if (max_t <= min_t) return false;

  Vector3D inters = r.o + (min_t+max_t)/2*d;
  if (inters[0] > max[0] || inters[0] < min[0] ||
      inters[1] > max[1] || inters[1] < min[1] ||
      inters[2] > max[2] || inters[2] < min[2]) return false;
  if (min_t > r.max_t) return false;

  t0 = min_t;
  t1 = max_t;

  return true;
  
}

void BBox::draw(Color c) const {

  glColor4f(c.r, c.g, c.b, c.a);

	// top
	glBegin(GL_LINE_STRIP);
	glVertex3d(max.x, max.y, max.z);
  glVertex3d(max.x, max.y, min.z);
  glVertex3d(min.x, max.y, min.z);
  glVertex3d(min.x, max.y, max.z);
  glVertex3d(max.x, max.y, max.z);
	glEnd();

	// bottom
	glBegin(GL_LINE_STRIP);
  glVertex3d(min.x, min.y, min.z);
  glVertex3d(min.x, min.y, max.z);
  glVertex3d(max.x, min.y, max.z);
  glVertex3d(max.x, min.y, min.z);
  glVertex3d(min.x, min.y, min.z);
	glEnd();

	// side
	glBegin(GL_LINES);
	glVertex3d(max.x, max.y, max.z);
  glVertex3d(max.x, min.y, max.z);
	glVertex3d(max.x, max.y, min.z);
  glVertex3d(max.x, min.y, min.z);
	glVertex3d(min.x, max.y, min.z);
  glVertex3d(min.x, min.y, min.z);
	glVertex3d(min.x, max.y, max.z);
  glVertex3d(min.x, min.y, max.z);
	glEnd();

}

std::ostream& operator<<(std::ostream& os, const BBox& b) {
  return os << "BBOX(" << b.min << ", " << b.max << ")";
}

} // namespace CMU462

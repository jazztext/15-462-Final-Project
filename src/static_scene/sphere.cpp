#include "sphere.h"

#include <cmath>

#include  "../bsdf.h"
#include "../misc/sphere_drawing.h"

namespace CMU462 { namespace StaticScene {

bool Sphere::test(const Ray& r, double& t1, double& t2) const {

  Vector3D co = r.o - this->o;
  double tmp = -dot(r.d, co);
  double det = tmp*tmp - dot(co,co) + this->r2;
  if (det < 0.) return false;
  
  double d = sqrt(det);
  t1 = tmp - d;
  t2 = tmp + d;
  return (t2 > 0); 
}

bool Sphere::intersect(const Ray& r) const {
  double t1, t2;
  if (test(r, t1, t2)) {
    if ((t1 <= r.max_t && t1 >= r.min_t) || (t2 <= r.max_t && t2 >= r.min_t)) 
      return true;
  }
  return false;

}

bool Sphere::intersect(const Ray& r, Intersection *i) const {
  double t1 = 0, t2 = std::numeric_limits<double>::infinity();
  Vector3D N;
  if (test(r, t1, t2)) {
    if (t1 < r.max_t && t1 > r.min_t && t1 < i->t) {
      i->t = t1;
      i->primitive = this;
      N = this->normal(r.o + t1*r.d);
      i->n = (dot(N,r.d) > 0) ? -N : N;
      i->bsdf = this->get_bsdf();
      return true;
    } else if (t2 < r.max_t && t2 > r.min_t && t2 < i->t) {
      i->t = t2;
      i->primitive = this;
      N = this->normal(r.o + t2*r.d);
      i->n = (dot(N,r.d) > 0) ? -N : N;
      i->bsdf = this->get_bsdf();
      return true;
    }
  }
  return false;

}

void Sphere::draw(const Color& c) const {
  Misc::draw_sphere_opengl(o, r, c);
}

void Sphere::drawOutline(const Color& c) const {
    //Misc::draw_sphere_opengl(o, r, c);
}


} // namespace StaticScene
} // namespace CMU462

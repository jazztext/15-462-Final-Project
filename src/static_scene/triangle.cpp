#include "triangle.h"

#include "CMU462/CMU462.h"
#include "GL/glew.h"

namespace CMU462 { namespace StaticScene {

Triangle::Triangle(const Mesh* mesh, size_t v1, size_t v2, size_t v3) :
    mesh(mesh), v1(v1), v2(v2), v3(v3) { }

BBox Triangle::get_bbox() const {
  
  // TODO: 
  // compute the bounding box of the triangle
  
  BBox out = BBox();
  out.expand(mesh->positions[v1]);
  out.expand(mesh->positions[v2]);
  out.expand(mesh->positions[v3]);

  return out;
}

bool Triangle::intersect(const Ray& r) const {
  
  // TODO: implement ray-triangle intersection
  
  Vector3D p0 = mesh->positions[v1];
  Vector3D p1 = mesh->positions[v2];
  Vector3D p2 = mesh->positions[v3];
  Vector3D e1 = p1 - p0;
  Vector3D e2 = p2 - p0;
  Vector3D s = r.o - p0;
  Vector3D d = r.d;

  Vector3D exd = cross(e1, d);
  Vector3D sxe = cross(s, e2);
  double inv = dot(exd, e2);
  if (inv == 0.0) return false;
  Vector3D n = Vector3D(dot(sxe, d), dot(exd, s), dot(sxe, e1)) / inv;
  return (0 <= n[0] && n[0] <= 1 &&
          0 <= n[1] && n[1] <= 1 &&
          r.min_t < n[2] && r.max_t > n[2]);
}

bool Triangle::intersect(const Ray& r, Intersection *i) const {
  
  // TODO: 
  // implement ray-triangle intersection. When an intersection takes
  // place, the Intersection data should be updated accordingly
  
  Vector3D p0 = mesh->positions[v1];
  Vector3D p1 = mesh->positions[v2];
  Vector3D p2 = mesh->positions[v3];
  Vector3D e1 = p1 - p0;
  Vector3D e2 = p2 - p0;
  Vector3D s = r.o - p0;
  Vector3D d = r.d;

  Vector3D exd = cross(e1, d);
  Vector3D sxe = cross(e2, s);
  double inv = dot(exd, e2);
  if (inv == 0.0) return false;
  Vector3D n = Vector3D(dot(sxe, d), dot(exd, s), dot(sxe, e1)) / inv;

  if (n[0] + n[1] > 1 || n[0] < 0 || n[1] < 0) return false;
  if (n[2] > r.max_t || n[2] < r.min_t) return false;

  if (i->t > n[2]) {
    i->t = n[2];
    i->primitive = this;
    Vector3D N1 = mesh->normals[v1];
    Vector3D N2 = mesh->normals[v2];
    Vector3D N3 = mesh->normals[v3];
    Vector3D N = (1-n[1]-n[0])*N1 + n[0]*N2 + n[1]*N3;
    N.normalize();
    i->n = (dot(N, d) >= 0) ? -N : N;
    i->bsdf = this->get_bsdf();
    return true;
  }

  return false;
}

void Triangle::draw(const Color& c) const {
  glColor4f(c.r, c.g, c.b, c.a);
  glBegin(GL_TRIANGLES);
  glVertex3d(mesh->positions[v1].x,
             mesh->positions[v1].y,
             mesh->positions[v1].z);
  glVertex3d(mesh->positions[v2].x,
             mesh->positions[v2].y,
             mesh->positions[v2].z);
  glVertex3d(mesh->positions[v3].x,
             mesh->positions[v3].y,
             mesh->positions[v3].z);
  glEnd();
}

void Triangle::drawOutline(const Color& c) const {
  glColor4f(c.r, c.g, c.b, c.a);
  glBegin(GL_LINE_LOOP);
  glVertex3d(mesh->positions[v1].x,
             mesh->positions[v1].y,
             mesh->positions[v1].z);
  glVertex3d(mesh->positions[v2].x,
             mesh->positions[v2].y,
             mesh->positions[v2].z);
  glVertex3d(mesh->positions[v3].x,
             mesh->positions[v3].y,
             mesh->positions[v3].z);
  glEnd();
}



} // namespace StaticScene
} // namespace CMU462

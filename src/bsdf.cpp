#include "bsdf.h"

#include <iostream>
#include <algorithm>
#include <utility>

using std::min;
using std::max;
using std::swap;

namespace CMU462 {

void make_coord_space(Matrix3x3& o2w, const Vector3D& n) {

    Vector3D z = Vector3D(n.x, n.y, n.z);
    Vector3D h = z;
    if (fabs(h.x) <= fabs(h.y) && fabs(h.x) <= fabs(h.z)) h.x = 1.0;
    else if (fabs(h.y) <= fabs(h.x) && fabs(h.y) <= fabs(h.z)) h.y = 1.0;
    else h.z = 1.0;

    z.normalize();
    Vector3D y = cross(h, z);
    y.normalize();
    Vector3D x = cross(z, y);
    x.normalize();

    o2w[0] = x;
    o2w[1] = y;
    o2w[2] = z;
}

// Diffuse BSDF //

Spectrum DiffuseBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return albedo * (1.0 / PI);
}

Spectrum DiffuseBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  *wi = CosineWeightedHemisphereSampler3D().get_sample(pdf);
  return this->f(wo, *wi);
}

// Mirror BSDF //

Spectrum MirrorBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum MirrorBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {

  // TODO:
  // Implement MirrorBSDF

  reflect(wo, wi);
  *pdf = 1;
  return (wo[2]) ? (this->reflectance * (1.0 / wo[2])) : Spectrum();
}

// Glossy BSDF //

/*
Spectrum GlossyBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum GlossyBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  *pdf = 1.0f;
  return reflect(wo, wi, reflectance);
}
*/

// Refraction BSDF //

Spectrum RefractionBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum RefractionBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {

  // TODO:
  // Implement RefractionBSDF

  return Spectrum();
}

// Glass BSDF //

Spectrum GlassBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum GlassBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {

  // TODO:
  // Compute Fresnel coefficient and either reflect or refract based on it.


  float n_i, n_t;
  float cosi = fabs(wo[2]);
  float sini = sqrt(1 - cosi*cosi);

  if (wo[2] > 0) { // entering material
   n_i = 1.0;
   n_t = this->ior;
  } else { // exiting material 
   n_i = this->ior;
   n_t = 1.0;
  }

  float sint = sini * (n_i / n_t);
  float cost = sqrt(1 - sint*sint);

  float r_par = (n_t*cosi - n_i*cost) / (n_t*cosi + n_i*cost);
  float r_prp = (n_i*cosi - n_t*cost) / (n_i*cosi + n_t*cost);

  float Fr = ((r_par * r_par) + (r_prp * r_prp)) / 2.0;


//  Fr = 0;

  if (UniformGridSampler2D().get_sample().x < Fr) { // reflect
    reflect(wo, wi);
    *pdf = 1;
    return (wo[2]) ? (this->reflectance * (1.0 /  fabs(wo[2]))) : Spectrum();
  } else { // refract
    if (refract(wo, wi, this->ior)) {
      *pdf = 1;
      return ((n_t*n_t)/(n_i*n_i))/fabs(cosi)*(this->transmittance);
    } else {
      return Spectrum(std::numeric_limits<double>::infinity(), 0, 0);
    }
  }

  return Spectrum();
}

void BSDF::reflect(const Vector3D& wo, Vector3D* wi) {

  // TODO:
  // Implement reflection of wo about normal (0,0,1) and store result in wi.
  *wi = Vector3D(-wo[0], -wo[1], wo[2]);
}

bool BSDF::refract(const Vector3D& wo, Vector3D* wi, float ior) {

  // TODO:
  // Use Snell's Law to refract wo surface and store result ray in wi.
  // Return false if refraction does not occur due to total internal reflection
  // and true otherwise. When dot(wo,n) is positive, then wo corresponds to a
  // ray entering the surface through vacuum.

  float n_i, n_t;
  float cosi = wo[2];
  float sini = sqrt(1 - cosi*cosi);

  Vector3D tmp = Vector3D(wo[0] / sini, wo[1] / sini, 1);

  if (wo[2] > 0) { // entering material
   n_i = 1.0;
   n_t = ior;
  } else { // exiting material 
   n_i = ior;
   n_t = 1.0;
  }

  float sint = sini * (n_i / n_t);
  float cost = sqrt(1 - sint*sint);
  if (wo[2] > 0) { 
   (*wi) = -Vector3D(tmp[0] * sint, tmp[1] * sint, cost);
  } else {
   (*wi) = Vector3D(-tmp[0] * sint, -tmp[1] * sint, cost);
  }

  float ior_rat = n_i / n_t;

  if ((ior_rat*ior_rat) * (1 - cosi*cosi) > 1) return false;

  return true;

}

// Emission BSDF //

Spectrum EmissionBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum EmissionBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  *wi  = sampler.get_sample(pdf);
  return Spectrum();
}

} // namespace CMU462

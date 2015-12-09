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

Spectrum DiffuseBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf, bool& a) {
  a = false;
  *wi = CosineWeightedHemisphereSampler3D().get_sample(pdf);
  return this->f(wo, *wi);
}

// Mirror BSDF //

Spectrum MirrorBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum MirrorBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf, bool& a) {
  a = false;
  reflect(wo, wi);
  *pdf = 1;
  return (wo[2]) ? (this->reflectance * (1.0 / wo[2])) : Spectrum();
}

// Glossy BSDF //


Spectrum GlossyBSDF::f(const Vector3D& wo, const Vector3D& wi) {
//  return reflectance * (1.0 / PI);
  return Spectrum();
}

Spectrum GlossyBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf, bool& a) {
  a = false;
/*  reflect(wo, wi);
  Vector3D w_perturb = roughness*CosineWeightedHemisphereSampler3D().get_sample(pdf);
  while (w_perturb[2] < 0) w_perturb = roughness*CosineWeightedHemisphereSampler3D().get_sample(pdf);
  (*wi) = (*wi) + w_perturb;
  wi->normalize();
  *pdf = (*pdf) * wi[2];
  return reflectance; */
  Vector3D w;
  reflect(wo, &w);
  Vector3D u = cross(Vector3D(0.00424, 1.0, 0.0076), w);
  u.normalize();
  Vector3D v = cross(u,w);
//  Vector3D sp = UniformHemisphereSampler3D().get_sample();
  Vector3D sp = CosineWeightedHemisphereSampler3D().get_sample(pdf);
  *wi = roughness * (sp[0] * u + sp[1] * v) + sp[2] * w;
  if ((*wi)[2] < 0) (*wi) = -(*wi);
  wi->normalize();
  //*pdf = ((roughness) ? (1.0 / PI) : 1.0);// * wo[2];
  return reflectance * (*wi)[2];
}

// Refraction BSDF //

Spectrum RefractionBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum RefractionBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf, bool& a) {

  // TODO:
  // Implement RefractionBSDF

  return Spectrum();
}

// Glass BSDF //

Spectrum GlassBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum GlassBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf, bool& inMat) {
  
  float n_i, n_t;
  float cosi = fabs(wo[2]);
  float sini = sqrt(1 - cosi*cosi);

  if (inMat) { // exiting material
   n_i = this->ior;
   n_t = 1.0;
  } else { // entering material 
   n_i = 1.0;
   n_t = this->ior;
  }

  float sint = sini * (n_i / n_t);
  float cost = sqrt(1 - sint*sint);

  float r_par = (n_t*cosi - n_i*cost) / (n_t*cosi + n_i*cost);
  float r_prp = (n_i*cosi - n_t*cost) / (n_i*cosi + n_t*cost);

  float Fr = ((r_par * r_par) + (r_prp * r_prp)) / 2.0;

  if (UniformGridSampler2D().get_sample().x < Fr) { // reflect
    reflect(wo, wi);
    *pdf = 1;
    return (wo[2]) ? (this->reflectance * (1.0 /  fabs(wo[2]))) : Spectrum();
  } else { // refract
    if (refract(wo, wi, this->ior, inMat)) {
      *pdf = 1;
      inMat = !inMat;
      return ((n_t*n_t)/(n_i*n_i))/fabs(cosi)*(this->transmittance);
    } else {
      reflect(wo, wi);
      *pdf = 1;
      return (wo[2]) ? (this->reflectance * (1.0 / fabs(wo[2]))) : Spectrum();
    }
  }

}

void BSDF::reflect(const Vector3D& wo, Vector3D* wi) {
  *wi = Vector3D(-wo[0], -wo[1], wo[2]);
}

bool BSDF::refract(const Vector3D& wo, Vector3D* wi, float ior, bool inMat) {

  float n_i, n_t;
  float cosi = wo[2];
  float sini = sqrt(1 - cosi*cosi);

  Vector3D tmp = Vector3D(wo[0] / sini, wo[1] / sini, 1);

  if (!inMat) { // entering material
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

Spectrum EmissionBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf, bool& a) {
  a = false;
  *wi  = sampler.get_sample(pdf);
  return Spectrum();
}


// Water BSDF //
WaterBSDF::WaterBSDF(Spectrum transmittance, Spectrum reflectance, float roughness, float ior) {
  this->transmittance = transmittance;
  this->reflectance = reflectance;
  this->roughness = roughness;
  this->ior = ior;
  this->glossy = new GlossyBSDF(reflectance, roughness);
}

Spectrum WaterBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum WaterBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf, bool& inMat) {
  
  float n_i, n_t;
  float cosi = fabs(wo[2]);
  float sini = sqrt(1 - cosi*cosi);

  if (inMat) { // exiting material
   n_i = this->ior;
   n_t = 1.0;
  } else { // entering material 
   n_i = 1.0;
   n_t = this->ior;
  }

  float sint = sini * (n_i / n_t);
  float cost = sqrt(1 - sint*sint);

  float r_par = (n_t*cosi - n_i*cost) / (n_t*cosi + n_i*cost);
  float r_prp = (n_i*cosi - n_t*cost) / (n_i*cosi + n_t*cost);

  float Fr = ((r_par * r_par) + (r_prp * r_prp)) / 2.0;

  if (UniformGridSampler2D().get_sample().x < Fr) { // reflect
    return glossy->sample_f(wo, wi, pdf, inMat);
//    reflect(wo, wi);
//    *pdf = 1;
//    return (wo
  } else { // refract
    if (refract(wo, wi, this->ior, inMat)) {
      *pdf = 1;
      inMat = !inMat;
      return ((n_t*n_t)/(n_i*n_i))/fabs(cosi)*(this->transmittance);
    } else {
      reflect(wo, wi);
      *pdf = 1;
      return (wo[2]) ? (this->reflectance * (1.0 / fabs(wo[2]))) : Spectrum();
    }
  }

}



} // namespace CMU462

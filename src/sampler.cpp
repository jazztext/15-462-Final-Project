#include "sampler.h"

namespace CMU462 {

// Uniform Sampler2D Implementation //

Vector2D UniformGridSampler2D::get_sample() const {

  // TODO:
  // Implement uniform 2D grid sampler

  double m = RAND_MAX;
  return Vector2D(std::rand() / m, std::rand() / m);
/*  std::random_device rd;
  std::mt19937 gen(rd());
  return Vector2D(std::generate_canonical<double, 10>(gen),
                  std::generate_canonical<double, 10>(gen));  
  double m = RAND_MAX;
  return Vector2D(std::rand() / m, std::rand() / m); */

}

// Uniform Hemisphere Sampler3D Implementation //

Vector3D UniformHemisphereSampler3D::get_sample() const {

  // TODO:
  // Implement uniform 3D hemisphere sampler

  Vector2D r = UniformGridSampler2D().get_sample();
  double x1 = r.x, x2 = 2.0 * M_PI * r.y;
  double s = sqrt(1 - x1*x1);
  return Vector3D(s*cos(x2), s*sin(x2), x1);

}

Vector3D CosineWeightedHemisphereSampler3D::get_sample() const {
  float f;
  return get_sample(&f);
}

Vector3D CosineWeightedHemisphereSampler3D::get_sample(float *pdf) const {
  Vector2D rtheta = UniformGridSampler2D().get_sample();
  float y = sqrt(1.0f - rtheta.x);
  rtheta.x = sqrt(rtheta.x);
  rtheta.y *= M_PI*2;
  Vector2D diskPt = Vector2D(rtheta.x*cos(rtheta.y), rtheta.x*sin(rtheta.y));
  *pdf = 1.0 / M_PI;
  return Vector3D(diskPt.x, diskPt.y, y);
}


} // namespace CMU462

#include "environment_light.h"
#include <iostream>

namespace CMU462 { namespace StaticScene {

EnvironmentLight::EnvironmentLight(const HDRImageBuffer* envMap)
    : envMap(envMap) {
	    // TODO: initialize things here as needed
  this->envMap = envMap;
  std::vector<float> illums = std::vector<float>();
  float total = 0.0f;
  for (int i=0;i<envMap->h;i++) {
    for (int j=0;j<envMap->w;j++) {
      float theta = ((float)i / envMap->h) * M_PI;
      float ill = envMap->data[i*envMap->w + j].illum();
      illums.push_back(ill * sin(theta));
      total += ill * sin(theta);
    }
  }
  for (int i=0;i<illums.size();i++) illums[i] /= total;
  std::vector<float> rows = std::vector<float>();
  std::vector<float> conds = illums;
  for (int row = 0; row < envMap->h; row++) {
    rows.push_back(0.0f);   
    for (int col = 0; col < envMap->w; col++) {
      rows[row] += illums[row*envMap->w + col];
    }
    for (int col = 0; col < envMap->w; col++) {
      if (rows[row]) conds[row*envMap->w + col] /= rows[row];
      else conds[row*envMap->w + col] = 0;
    }
  }
  std::vector<float> CDF = rows;
  for (int row = 1; row < envMap->h; row++) {
    CDF[row] += CDF[row-1];
  }
  CDF[envMap->h - 1] = std::max(CDF[envMap->h - 1], 1.0f);
  std::vector<std::vector<float>> CDF_cols = std::vector<std::vector<float>>();
  for (int row = 0; row < envMap->h; row++) {
    CDF_cols.push_back(std::vector<float>(1, conds[row*envMap->w]));
    for (int col = 1; col < envMap->w; col++) {
      int i = row*envMap->w + col;
      CDF_cols[row].push_back(conds[i] + conds[i - 1]);
    }
    CDF_cols[row][envMap->w - 1] = std::max(CDF_cols[row][envMap->w - 1], 1.0f);
  }
  this->CDF_row = CDF; 
  this->CDF_col = CDF_cols;
}
  

Spectrum EnvironmentLight::sample_L(const Vector3D& p, Vector3D* wi,
                                    float* distToLight,
                                    float* pdf) const {
  // TODO: Implement

  Vector2D tmp = UniformGridSampler2D().get_sample();
  int row = std::lower_bound(this->CDF_row.begin(), this->CDF_row.end(), tmp.x) - this->CDF_row.begin();
  float theta = ((float)row / this->CDF_row.size()) * M_PI;
  std::vector<float> CDF = CDF_col[row];
  int col = std::lower_bound(CDF.begin(), CDF.end(), tmp.y) - CDF.begin();
  float phi = 2.0 * M_PI * ((float)col / CDF.size());
  *distToLight = std::numeric_limits<float>::infinity();
  float pr_row = (row) ? (CDF_row[row] - CDF_row[row-1]) : CDF_row[row];
  float pr_col = (col) ? (CDF[col] - CDF[col-1]) : CDF[col]; 
  float pr = pr_row * pr_col;
  *pdf = pr * (envMap->w * envMap->h) / (4.0 * M_PI);
  *wi = Vector3D(cos(phi)*sin(theta), sin(phi)*sin(theta), cos(theta));
  return sample_dir(Ray(Vector3D(), *wi));
}

Spectrum EnvironmentLight::sample_dir(const Ray& r) const {
  // TODO: Implement
  float cos_t = r.d[1];
  float sin_t = sqrt(1 - cos_t*cos_t);
  float theta = acos(r.d[1]);
  float sinp = (theta != 0) ? (r.d[0] / sin_t) : 1;
  float cosp = (theta != 0) ? (r.d[2] / sin_t) : 1;

  sinp = clamp(sinp, -1, 1);
  cosp = clamp(cosp, -1, 1);


  float phi = (sinp >= 0) ? acos(cosp) : (2.0 * M_PI - acos(cosp));

  float v = phi / (2.0 * M_PI);
  float u = theta / M_PI;

  float ux = round(u * envMap->h);
  float vx = round(v * envMap->w);
  float u0 = ux - 0.5f, u1 = ux + 0.5f;
  float v0 = vx - 0.5f, v1 = vx + 0.5f;
  if (v0 < 0) v0 = v1;
  if (v1 >= envMap->h - 1.0f) v1 = v0;

  float du = fabs(u*envMap->h - u0);
  float dv = fabs(v*envMap->w - v0);

  int u_0 = (int)u0 * envMap->w;
  int u_1 = (int)u1 * envMap->w;
  int v_0 = (int)v0 % envMap->w;
  int v_1 = (int)v1 % envMap->w;

  Spectrum s1 = Spectrum(), s2 = Spectrum();
  s1 += (1-du)*envMap->data[u_0 + v_0];
  s2 += (1-du)*envMap->data[u_0 + v_1];
  s1 += du*envMap->data[u_1 + v_0];
  s2 += du*envMap->data[u_1 + v_1];

  return (1-dv)*s1 + dv*s2;
}

} // namespace StaticScene
} // namespace CMU462

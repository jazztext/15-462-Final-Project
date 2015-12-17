#include "halfEdgeMesh.h"
#include "student_code.h"
#include "kdTree.h"

double Halfedge::cot(Eigen::MatrixXd &points)
{
  int iv1 = vertex()->index;
  int iv2 = twin()->vertex()->index;
  int iv3 = next()->twin()->vertex()->index;
  Eigen::Vector3d v1 = points.row(iv1);
  Eigen::Vector3d v2 = points.row(iv2);
  Eigen::Vector3d v3 = points.row(iv3);
  Eigen::Vector3d a = v1 - v3, b = v2 - v3;
  return a.dot(b) / (a.cross(b).norm());
}

double Face::area(Eigen::MatrixXd &points)
{
  if (*((void **) &halfedge()) == NULL) return 0;
  int iv1 = halfedge()->vertex()->index;
  int iv2 = halfedge()->next()->vertex()->index;
  int iv3 = halfedge()->next()->next()->vertex()->index;
  Eigen::Vector3d v1 = points.row(iv1);
  Eigen::Vector3d v2 = points.row(iv2);
  Eigen::Vector3d v3 = points.row(iv3);
  Eigen::Vector3d a = v2 - v1, b = v3 - v1;
  if (a.cross(b).norm() == 0) {
    VertexIter ve1 = halfedge()->vertex();
    VertexIter ve2 = halfedge()->next()->vertex();
    VertexIter ve3 = halfedge()->next()->next()->vertex();
  }
  return a.cross(b).norm() / 2;
}

double Face::volume(Eigen::MatrixXd &points)
{
  int iv1 = halfedge()->vertex()->index;
  int iv2 = halfedge()->next()->vertex()->index;
  int iv3 = halfedge()->next()->next()->vertex()->index;
  Eigen::Vector3d v1 = points.row(iv1);
  Eigen::Vector3d v2 = points.row(iv2);
  Eigen::Vector3d v3 = points.row(iv3);
  double result = (-v3[0]*v2[1]*v1[2] + v2[0]*v3[1]*v1[2] + v3[0]*v1[1]*v2[2] -
                    v1[0]*v3[1]*v2[2] - v2[0]*v1[1]*v3[2] + v1[0]*v2[1]*v3[2]);
  return result / 6;
}

double Vertex::dualArea(Eigen::MatrixXd &points)
{
  double total = 0;
  HalfedgeIter he = halfedge();
  do {
    total += he->face()->area(points);
    he = he->twin()->next();
  } while (he != halfedge());
  return total / 3;
}

double HalfedgeMesh::volume()
{
  double v = 0;
  for (FaceIter f = facesBegin(); f != facesEnd(); f++)
    v += f->volume(positions);
  return v;
}

double HalfedgeMesh::surfaceArea()
{
  double a = 0;
  for (FaceIter f = facesBegin(); f != facesEnd(); f++)
    a += f->area(positions);
  return a;
}

void HalfedgeMesh::correctVolume()
{
  double currentVolume = volume();
  double sa = surfaceArea();
  double h = (_V - currentVolume) / sa;
  double precision = .005;
  int count = 0;
  while (fabs(_V - currentVolume) / _V > precision) {
    for (VertexIter v = verticesBegin(); v != verticesEnd(); v++) {
      positions.row(v->index) += h * normals.row(v->index);
    }
    double newVolume = volume();
    if ((_V - currentVolume) * (_V - newVolume) < 0) {
      h = -h / 2;
    }
    currentVolume = newVolume;
  }
}

void HalfedgeMesh::indexVertices()
{
  int i = 0;
  for (VertexIter v = verticesBegin(); v!= verticesEnd(); v++) {
    v->index = i++;
  }
  positions = Eigen::MatrixXd(nVertices(), 3);
  normals = Eigen::MatrixXd(nVertices(), 3);
  vels = Eigen::VectorXd(nVertices());
  for (VertexIter v = verticesBegin(); v != verticesEnd(); v++) {
    for (i = 0; i < 3; i++) {
      positions(v->index, i) = v->position[i];
      normals(v->index, i) = v->normal[i];
      vels[v->index] = v->velocity;
    }
  }
}

Eigen::SparseMatrix<double> HalfedgeMesh::areas(bool inverted)
{
  std::vector<Eigen::Triplet<double> > coeffs;
  for (VertexIter v = verticesBegin(); v != verticesEnd(); v++) {
    double val = inverted? 1 / v->dualArea(positions) : v->dualArea(positions);
    coeffs.push_back(Eigen::Triplet<double>(v->index, v->index, val));
  }
  Eigen::SparseMatrix<double> A(nVertices(), nVertices());
  A.setFromTriplets(coeffs.begin(), coeffs.end());
  return A;
}


Eigen::SparseMatrix<double> HalfedgeMesh::laplacian(bool debug)
{
  std::vector<Eigen::Triplet<double> > coeffs;
  for (VertexIter v = verticesBegin(); v != verticesEnd(); v++) {
    HalfedgeIter he = v->halfedge();
    int i = v->index;
    double iCoeff = 0;
    do {
      int j = he->twin()->vertex()->index;
      double cotAlpha = he->cot(positions);
      double cotBeta = he->twin()->cot(positions);
      double val = (cotAlpha + cotBeta) / 2;
      coeffs.push_back(Eigen::Triplet<double>(i, j, val));
      iCoeff -= val;
      he = he->twin()->next();
    } while (he != v->halfedge());
    coeffs.push_back(Eigen::Triplet<double>(i, i, iCoeff));
  }
  Eigen::SparseMatrix<double> L(nVertices(), nVertices());
  L.setFromTriplets(coeffs.begin(), coeffs.end());
  return L;
}

void HalfedgeMesh::initWaveEquation(double dt, double c)
{
  displacements = Eigen::VectorXd::Zero(nVertices());
  for (VertexIter v = verticesBegin(); v != verticesEnd(); v++)
    v->velocity = 0;
  indexVertices();
  _L = laplacian();
  _dt = dt;
  _c = c;
}

void HalfedgeMesh::initSurfaceTension(double dt, double sigma)
{
  double dTau = dt / 60;
  initVCF(dTau, sigma);
  initWaveEquation(dt/10, 25);
}

void HalfedgeMesh::initVCF(double dt, double gamma)
{
  indexVertices();
  _V = volume();
  _dtVCF = dt;
  _gamma = gamma;
}

void HalfedgeMesh::stepVCF(bool debug)
{
  _L = laplacian(debug);
  Eigen::MatrixXd k = _L * positions;
  double kAvg = 0;
  double areaSum = 0;
  for (VertexIter v = verticesBegin(); v != verticesEnd(); v++) {
    double a = v->dualArea(positions);
    for (int i = 0; i < 3; i++) k(v->index, i) /= a;
    kAvg += a * sqrt(pow(k(v->index, 0), 2) + pow(k(v->index, 1), 2) +
                     pow(k(v->index, 2), 2));
    areaSum += a;
  }
  kAvg /= areaSum;
  for (int i = 0; i < k.rows(); i++) {
    Vector3D kn(k(i, 0), k(i, 1), k(i, 2));
    Vector3D normal = kn.unit();
    Vector3D update = _gamma * _dtVCF * (kn - kAvg*normal);
    for (int n = 0; n < 3; n++) {
      positions(i, n) += update[n];
    }
  }
}

void HalfedgeMesh::stepSurfaceTension()
{
  static int thing = 0;
  //std::cout << thing++ << " " << volume() << " " << nVertices() << "\n";
  displacements = Eigen::VectorXd::Zero(nVertices());
  Eigen::MatrixXd oldPositions(positions);
  for (int i = 0; i < 60; i++) {
    stepVCF(false);
  }

  for (VertexIter v = verticesBegin(); v != verticesEnd(); v++) {
    Eigen::VectorXd xi = oldPositions.row(v->index);
    int ind;
    double bestDist = -1;
    for (int i = 0; i < positions.rows(); i++) {
      double dist = sqrt(pow(positions(i, 0) - xi[0], 2) +
                         pow(positions(i, 1) - xi[1], 2) +
                         pow(positions(i, 2) - xi[2], 2));
      if (bestDist == -1 || dist < bestDist) {
        ind = i;
        bestDist = dist;
      }
    }
    Eigen::VectorXd xS = positions.row(ind);
    Eigen::VectorXd n = normals.row(v->index);
    double newDisplacement = (xi - xS).transpose() * n;
    displacements[v->index] = newDisplacement;
    oldPositions.row(v->index) = xi - newDisplacement * n;
  }
  positions = oldPositions;
  for (int i = 0; i < 10; i++) stepWaveEquation();
  for (VertexIter v = verticesBegin(); v != verticesEnd(); v++) {
    Eigen::Vector3d newPos = positions.row(v->index)
                             + displacements[v->index] * normals.row(v->index);
    for (int i = 0; i < 3; i++) v->position[i] = newPos[i];
    v->velocity = vels[v->index];
  }

}

void HalfedgeMesh::stepWaveEquation()
{
  Eigen::SparseMatrix<double> Area = areas(false);
  Eigen::SparseMatrix<double> AreaInv = areas(true);
  Eigen::SparseMatrix<double> A = Area - _dt * _dt * _c * _c / 4 * _L;
  Eigen::VectorXd b = Area * displacements + Area * _dt * vels + _dt*_dt*_c*_c/4*_L * displacements;
  Eigen::ConjugateGradient<Eigen::SparseMatrix<double> > cg;
  cg.compute(A);
  Eigen::VectorXd newDisplacements = cg.solve(b);
  Eigen::VectorXd newVels = vels + _dt*_c*_c/2 * AreaInv * _L * (displacements +
                                                   newDisplacements);
  //damping
  newDisplacements *= .999;
  displacements = newDisplacements;
  vels = newVels;
}

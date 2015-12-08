#include "halfedgeMesh.h"

double Halfedge::cot()
{
  Vector3D a = vertex()->position - next()->twin()->vertex()->position;
  Vector3D b = twin()->vertex()->position-next()->twin()->vertex()->position;
  return dot(a, b) / (cross(a, b).norm());
}

double Face::area()
{
  Vector3D v1 = halfedge()->vertex()->position;
  Vector3D v2 = halfedge()->next()->vertex()->position;
  Vector3D v3 = halfedge()->next()->next()->vertex()->position;
  Vector3D a = v2 - v1, b = v3 - v1;
  return cross(a, b).norm() / 2;
}

double Vertex::dualArea()
{
  double total = 0;
  HalfedgeIter he = halfedge();
  do {
    total += he->face()->area();
    he = he->twin()->next();
  } while (he != halfedge());
  return total;
}

void HalfedgeMesh::indexVertices()
{
  int i = 0;
  for (VertexIter v = verticesBegin(); v!= verticesEnd(); v++) {
    v->index = i++;
  }
}

Eigen::SparseMatrix<double> HalfedgeMesh::laplacian()
{
  std::vector<Eigen::Triplet<double> > coeffs;
  for (VertexIter v = verticesBegin(); v != verticesEnd(); v++) {
    HalfedgeIter he = v->halfedge();
    int i = v->index;
    double iCoeff = 0;
    do {
      int j = he->twin()->vertex()->index;
      double cotAlpha = he->cot();
      double cotBeta = he->twin()->cot();
      double val = (he->cot() + he->twin()->cot()) / 2;
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
  indexVertices();
  displacements = Eigen::VectorXd::Zero(nVertices());
  vels = Eigen::VectorXd::Zero(nVertices());
  for (VertexIter v = verticesBegin(); v != verticesEnd(); v++) {
    v->initialPosition = v->position;
    v->initialNormal = v->normal;
  }
  _L = laplacian();
  _dt = dt;
  _c = c;
}

void HalfedgeMesh::initVCF(double dt, double gamma)
{
  indexVertices();
  _dt = dt;
  _gamma = gamma;
  positions = Eigen::MatrixXd(nVertices(), 3);
  for (VertexIter v = verticesBegin(); v != verticesEnd(); v++) {
    positions(v->index, 0) = v->position.x;
    positions(v->index, 1) = v->position.y;
    positions(v->index, 2) = v->position.z;
  }
}

void HalfedgeMesh::stepVCF()
{
  Eigen::MatrixXd L = laplacian();
  Eigen::MatrixXd k = L * positions;
  for (VertexIter v = verticesBegin(); v != verticesEnd(); v++) {
    for (int i = 0; i < 3; i++) k(v->index, i) /= v->dualArea();
  }
  Eigen::MatrixXd dVAvg(nVertices(), 3);
  for (VertexIter v = verticesBegin(); v != verticesEnd(); v++) {
    double areaSum = v->dualArea();
    HalfedgeIter he = v->halfedge();
    do {
      VertexIter vNeighbor = he->twin()->vertex();
      areaSum += vNeighbor->dualArea();
      he = he->twin()->next();
    } while (he != v->halfedge());
    for (int i = 0; i < 3; i++)
      dVAvg(v->index, i) = v->dualArea() / areaSum * k(v->index, i);
  }
  Eigen::MatrixXd kAvg(nVertices(), 3);
  for (VertexIter v = verticesBegin(); v != verticesEnd(); v++) {
    double dVSum = sqrt(pow(dVAvg(v->index, 0), 2) +
                        pow(dVAvg(v->index, 1), 2) +
                        pow(dVAvg(v->index, 2), 2));
    Eigen::Vector3d n;
    for (int i = 0; i < 3; i++) n[i] = dVAvg(v->index, i) / dVSum;
    HalfedgeIter he = v->halfedge();
    do {
      VertexIter vNeighbor = he->twin()->vertex();
      dVSum += sqrt(pow(dVAvg(vNeighbor->index, 0), 2) +
                    pow(dVAvg(vNeighbor->index, 1), 2) +
                    pow(dVAvg(vNeighbor->index, 2), 2));
      he = he->twin()->next();
    } while (he != v->halfedge());
    for (int i = 0; i < 3; i++)
      kAvg(v->index, i) = dVSum * n[i];
  }
  positions += _gamma * _dt * (k - kAvg);
  for (VertexIter v = verticesBegin(); v != verticesEnd(); v++) {
    v->position.x = positions(v->index, 0);
    v->position.y = positions(v->index, 1);
    v->position.z = positions(v->index, 2);
  }
}

void HalfedgeMesh::stepWaveEquation()
{
  Eigen::SparseMatrix<double> I(nVertices(), nVertices());
  I.setIdentity();
  Eigen::SparseMatrix<double> A = I - _dt * _dt * _c * _c / 4 * _L;
  Eigen::VectorXd b = displacements + _dt * vels + _dt*_dt*_c*_c/4*_L * displacements;
  Eigen::ConjugateGradient<Eigen::SparseMatrix<double> > cg;
  cg.compute(A);
  Eigen::VectorXd newDisplacements = cg.solve(b);
  Eigen::VectorXd newVels = vels + _dt*_c*_c/2 * _L * (displacements +
                                                   newDisplacements);
  displacements = newDisplacements;
  vels = newVels;

  for (VertexIter v = vertices.begin(); v != vertices.end(); v++) {
    v->position = v->initialPosition + displacements[v->index]*v->initialNormal;
  }
  for (VertexIter v = vertices.begin(); v != vertices.end(); v++) {
    v->computeNormal();
  }
}

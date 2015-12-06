#include "bvh.h"

#include "CMU462/CMU462.h"
#include "static_scene/triangle.h"

#include <iostream>
#include <stack>

using namespace std;

namespace CMU462 { namespace StaticScene {
// @TODO: clean this up
unsigned int morton3D(Vector3D);
typedef float red_func(float);
typedef std::pair<unsigned int, Vector3D> morton_pt;

float aac_reduce(float x) {
  return 0.6*pow(x, 0.6);
}

/*

BVHNode *build_tree(std::vector<Primitive*> &prims, int start, int end, int leaf_prims, red_func f) {
  if (prims.size() < leaf_prims) {
    BVHNode *root = BVHNode(start, end - start, 
*/

/* Stuff for highest-level BVH function (AAC in paper) */


bool pair_cmp(morton_pt &p1, morton_pt &p2) {
  if (p1.first < p2.first) return -1;
  else if (p1.first == p2.first) return 0;
  return 1;
}
  

BVHNode *make_bvh_aac(std::vector<Primitive*> &prims, int leaf_prims, red_func f) {
  std::vector<Vector3D> centroids;
  int n = prims.size();
  centroids.resize(n);
  BBox total = BBox();
  for (int i = 0; i < n; i++) {
    BBox tmp = prims[i]->get_bbox();
    total.expand(tmp);
    centroids[i] = (tmp.max + tmp.min) / 2.0;
  }
  Vector3D span = total.max - total.min;
  Vector3D span_inv = Vector3D(1.0 / span[0], 1.0 / span[1], 1.0 / span[2]);
  Vector3D adj = Vector3D(total.min[0]*span_inv[0],total.min[1]*span_inv[1],total.min[2]*span_inv[2]);
  for (int i = 0; i < n; i++) {
    centroids[i][0] = centroids[i][0] * span_inv[0] - adj[0];
    centroids[i][1] = centroids[i][1] * span_inv[1] - adj[1];
    centroids[i][2] = centroids[i][2] * span_inv[2] - adj[2];
  }
  std::vector<morton_pt> cent_zs;
  cent_zs.resize(n);
  for (int i = 0; i < n; i++) {
    cent_zs[i].first = morton3D(centroids[i]);
    cent_zs[i].second = centroids[i];
  }

  // TODO: Maybe implement radix sort here?
  std::sort(cent_zs.begin(), cent_zs.end(), pair_cmp); 
  
  return NULL;
}

/* This code copied from NVidia devblogs */

// Expands a 10-bit integer into 30 bits
// by inserting 2 zeros after each bit.
unsigned int expandBits(unsigned int v) {
  v = (v * 0x00010001u) & 0xFF0000FFu;
  v = (v * 0x00000101u) & 0x0F00F00Fu;
  v = (v * 0x00000011u) & 0xC30C30C3u;
  v = (v * 0x00000005u) & 0x49249249u;
  return v;
}

// Calculates a 30-bit Morton code for the
// given 3D point located within the unit cube [0,1].
unsigned int morton3D(Vector3D v) {
  float x = std::min(std::max(v[0] * 1024.0, 0.0), 1023.0);
  float y = std::min(std::max(v[1] * 1024.0, 0.0), 1023.0);
  float z = std::min(std::max(v[2] * 1024.0, 0.0), 1023.0);
  unsigned int xx = expandBits((unsigned int)x);
  unsigned int yy = expandBits((unsigned int)y);
  unsigned int zz = expandBits((unsigned int)z);
  return xx * 4 + yy * 2 + zz;
}

/* End NVidia code */

double surfaceAreaCost(BVHNode *node, vector<BBox>& bboxes,
                       vector<int>& primCounts, int bin, int numBins)
{
  double sn = node->bb.surface_area();
  BBox bboxA, bboxB;
  int na = 0, nb = 0;
  for (int i = 0; i < numBins; i++) {
    if (i < bin) {
      bboxA.expand(bboxes[i]);
      na += primCounts[i];
    }
    else {
      bboxB.expand(bboxes[i]);
      nb += primCounts[i];
    }
  }
  if (na == 0 || nb == 0) return INF_D;
  return bboxA.surface_area() / sn * na + bboxB.surface_area() / sn * nb;
}

void BVHAccel::splitNode(BVHNode *node, int maxLeafSize) {
  //base case
  if (node->range <= maxLeafSize) {
    return;
  }
  //choose plane to split along
  double smallestBins = INF_D;
  int numBins = 32;
  vector<int> binAssignments;
  int plane;
  double bestCost = INF_D;
  for (int n = 0; n < 3; n++) {
    double start = INF_D, end = -INF_D;
    for (int i = 0; i < node->range; i++) {
      double p = primitives[node->start + i]->get_bbox().centroid()[n];
      if (p < start) start = p;
      if (p > end) end = p;
    }
    double binSize = (end - start) / numBins * 1.0001;
    if (binSize == 0) {
    }
    else if (binSize > 0) {
    vector<BBox> binBBoxes(numBins);
    vector<int> binPrimCounts(numBins);
    vector<int> currentBinAssignments(node->range);
    for (unsigned i = 0; i < node->range; i++) {
      BBox bbox = primitives[node->start + i]->get_bbox();
      int bin = floor((bbox.centroid()[n] - start) / binSize);
      binBBoxes[bin].expand(bbox);
      binPrimCounts[bin]++;
      currentBinAssignments[i] = bin;
    }
    for (int i = 1; i < numBins; i++) {
      double cost = surfaceAreaCost(node, binBBoxes, binPrimCounts, i, numBins);
      if (cost < bestCost) {
        plane = i;
        bestCost = cost;
        binAssignments = currentBinAssignments;
      }
    }
    if (binSize < smallestBins) smallestBins = binSize;
    }
  }
  if (bestCost == INF_D) {
    for (int i = 0; i < node->range; i++) {
      if (i < node->range / 2) binAssignments.push_back(0);
      else binAssignments.push_back(1);
      plane = 1;
    }
  }
  //perform split
  //rearrange primitives
  int i = 0, j = node->range - 1;
  BBox lbox, rbox;
  while (i < j) {
    while (binAssignments[i] < plane && i <= j)
      lbox.expand(primitives[node->start + i++]->get_bbox());
    while (binAssignments[j] >= plane && j >= i)
      rbox.expand(primitives[node->start + j--]->get_bbox());
    if (i < j) {
      Primitive *tmp = primitives[node->start + i];
      primitives[node->start + i] = primitives[node->start + j];
      primitives[node->start + j] = tmp;
      int tmpBin = binAssignments[i];
      binAssignments[i] = binAssignments[j];
      binAssignments[j] = tmpBin;
    }
  }
  //create child nodes
  node->l = new BVHNode(lbox, node->start, i);
  node->r = new BVHNode(rbox, i + node->start, node->range - i);
  splitNode(node->l, maxLeafSize);
  splitNode(node->r, maxLeafSize);
}

BVHAccel::BVHAccel(const std::vector<Primitive *> &_primitives,
                   size_t max_leaf_size) {
  this->primitives = _primitives;
  BBox bb;
  for (size_t i = 0; i < primitives.size(); i++) {
    bb.expand(primitives[i]->get_bbox());
  }
  root = new BVHNode(bb, 0, primitives.size());
  splitNode(root, max_leaf_size);
}

void destroyBVHNode(BVHNode *node) {
  if (node->l) destroyBVHNode(node->l);
  if (node->r) destroyBVHNode(node->r);
  delete node;
}

BVHAccel::~BVHAccel() {

  // TODO:
  // Implement a proper destructor for your BVH accelerator aggregate
  destroyBVHNode(this->root);

}

BBox BVHAccel::get_bbox() const {
  return root->bb;
}

bool BVHAccel::intersect(const Ray &ray) const {

  // TODO:
  // Implement ray - bvh aggregate intersection test. A ray intersects
  // with a BVH aggregate if and only if it intersects a primitive in
  // the BVH that is not an aggregate.

  Intersection i;
  return intersect(ray, &i);
}

bool BVHAccel::intersectNode(BVHNode *node, const Ray &ray, Intersection *i) const {
  bool hit = false;
  if (node->isLeaf()) {
    for (int n = node->start; n < node->start + node->range; n++) {
      if (primitives[n]->intersect(ray, i)) hit = true;
    }
    return hit;
  }
  double minTL = ray.min_t, minTR = ray.min_t;
  double maxTL = ray.max_t, maxTR = ray.max_t;
  bool hitL = node->l->bb.intersect(ray, minTL, maxTL);
  bool hitR = node->r->bb.intersect(ray, minTR, maxTR);
  BVHNode *first, *second;
  bool hitFirst, hitSecond;
  double minTS;
  if (minTL < minTR) {
    first = node->l; second = node->r; hitFirst = hitL; hitSecond = hitR; minTS = minTR;
  } else {
    first = node->r; second = node->l; hitFirst = hitR; hitSecond = hitL; minTS = minTL;
  }

  if (hitFirst && intersectNode(first, ray, i)) hit = true;
  if (hitSecond && minTS < ray.max_t && intersectNode(second, ray, i)) hit = true;
  
  return hit;
}

bool BVHAccel::intersect(const Ray &ray, Intersection *i) const {

  // TODO:
  // Implement ray - bvh aggregate intersection test. A ray intersects
  // with a BVH aggregate if and only if it intersects a primitive in
  // the BVH that is not an aggregate. When an intersection does happen.
  // You should store the non-aggregate primitive in the intersection data
  // and not the BVH aggregate itself.
  //

  return intersectNode(root, ray, i);
}

}  // namespace StaticScene
}  // namespace CMU462

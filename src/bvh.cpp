#include "bvh.h"

#include "CMU462/CMU462.h"
#include "static_scene/triangle.h"

#include <iostream>
#include <stack>

using namespace std;

namespace CMU462 { namespace StaticScene {


void split_bvhnode(BVHNode *root, size_t max_leaf_size, std::vector<Primitive*> &prims) {
  if (root->range <= max_leaf_size) return;
  int start = root->start;
  int end = start + root->range;

 
  std::vector<Vector3D> centroids;
  for (int i = start;i < end; i++) {
    BBox tmp = prims[i]->get_bbox();
    centroids.push_back((tmp.max + tmp.min) / 2.0);
  }


  int buckets = 32;

  std::vector<BBox> cuts;
  std::vector<int> cut_counts;
  cuts.resize(buckets*buckets*buckets);
  cut_counts.assign(buckets*buckets*buckets, 0);

  float best_cost = std::numeric_limits<float>::max();
  int best_i = -1, best_dim = -1;
  BBox best_l, best_r;

  for (int dim = 0; dim < 3; dim++) {
    double dd = (root->bb.max[dim] - root->bb.min[dim]) / buckets;
    for (int p = 0; p < root->range;p++) {
      int bucket_i = abs(centroids[p][dim] - root->bb.min[dim]) / dd;
      bucket_i = max(min(bucket_i, buckets-1), 0);
      cuts[dim*buckets + bucket_i].expand(prims[p + start]->get_bbox());
      cut_counts[dim*buckets + bucket_i]++;
    }
    
    for (int b = 1; b < buckets; b++) {
      BBox lower, upper;
      int lower_c = 0, upper_c = 0;
      for (int i = 0; i < b; i++) {
        lower.expand(cuts[dim*buckets + i]);
        lower_c += cut_counts[dim*buckets+i];
      }
      for (int i=b; i < buckets; i++) {
        upper.expand(cuts[dim*buckets + i]);
        upper_c += cut_counts[dim*buckets+i];
      }
      float cost = lower_c*(lower.surface_area() / root->bb.surface_area()) + 
                   upper_c*(upper.surface_area() / root->bb.surface_area());
      if (cost < best_cost && !lower.empty() && !upper.empty()) {
        best_cost = cost;
        best_i = b; best_dim = dim;
        best_l = lower; best_r = upper;
      }
    }
  }

  std::vector<Primitive*> prims_lower;
  std::vector<Primitive*> prims_upper;

  double dd = (root->bb.max[best_dim] - root->bb.min[best_dim]) / buckets;
  for (int p = 0; p < root->range;p++) {
    int bucket_i = abs(centroids[p][best_dim] - root->bb.min[best_dim]) / dd;
    bucket_i = max(min(bucket_i, buckets-1), 0);
    ((bucket_i < best_i) ? prims_lower : prims_upper).push_back(prims[start+p]);
  }

  if (prims_lower.size() == 0 || prims_upper.size() == 0) {
    prims_lower.clear();
    prims_upper.clear();
    for (int p=0;p<root->range;p++) {
      ((p < root->range/2) ? prims_lower : prims_upper).push_back(prims[start+p]);
    }
  }

  for (int i=0;i<prims_lower.size();i++) prims[start+i] = prims_lower[i];
  for (int i=0;i<prims_upper.size();i++) prims[start+i+prims_lower.size()] = prims_upper[i];

  root->l = new BVHNode(best_l, start, prims_lower.size());
  root->r = new BVHNode(best_r, start + prims_lower.size(), prims_upper.size());


  split_bvhnode(root->l, max_leaf_size, prims);
  split_bvhnode(root->r, max_leaf_size, prims);

}

BVHAccel::BVHAccel(const std::vector<Primitive *> &_primitives,
                   size_t max_leaf_size) {

  // TODO:
  // Construct a BVH from the given vector of primitives and maximum leaf
  // size configuration. The starter code build a BVH aggregate with a
  // single leaf node (which is also the root) that encloses all the
  // primitives.

  std::vector<Primitive*> prims = _primitives;

  BBox bb;
  for (size_t i = 0; i < prims.size(); ++i) {
    bb.expand(prims[i]->get_bbox());
  }

  root = new BVHNode(bb, 0, prims.size());
  split_bvhnode(root, max_leaf_size, prims);

  this->primitives = prims;

}

void destroyBVHNode(BVHNode *root) {
  if (root->l != NULL) destroyBVHNode(root->l);
  if (root->r != NULL) destroyBVHNode(root->r);
  delete root;
}

BVHAccel::~BVHAccel() {

  // TODO:
  // Implement a proper destructor for your BVH accelerator aggregate
  destroyBVHNode(this->root);

}

BBox BVHAccel::get_bbox() const {
  return root->bb;
}

bool bvhNodeIntersect(BVHNode *b, const Ray &ray, std::vector<Primitive*> prims) {
  double t1, t2;
  if (!b->bb.intersect(ray, t1,t2)) return false;
  if (b->l == NULL && b->r == NULL) {
    int end = b->start + b->range;
    bool hit = false;
    for (int i=b->start;i<end;i++) {
      hit = hit || prims[i]->intersect(ray);
    }
    return hit;
  } else {
    return bvhNodeIntersect(b->l, ray, prims) || bvhNodeIntersect(b->r, ray, prims);
  }
}

bool BVHAccel::intersect(const Ray &ray) const {

  // TODO:
  // Implement ray - bvh aggregate intersection test. A ray intersects
  // with a BVH aggregate if and only if it intersects a primitive in
  // the BVH that is not an aggregate.

  return bvhNodeIntersect(root, ray, primitives);

}


int glob_i = 0;

bool bvhIntersect(BVHNode *b, const Ray &ray, std::vector<Primitive*> prims, Intersection *i) {
  glob_i++;
  double t1 = std::numeric_limits<double>::infinity(), t2 = 0;
  if (!b->bb.intersect(ray, t1, t2)) return false;
  if (b->l == NULL && b->r == NULL) {
    int end = b->start + b->range;
    bool hit = false;
    for (int j=b->start;j<end;j++) {
      if (prims[j]->intersect(ray, i)) {
        hit = 1;
        ray.max_t = i->t;
      }
    }
    return hit;
  } else {
    double lt1 = std::numeric_limits<double>::infinity(), lt2 = 0;
    double rt1 = lt1, rt2 = lt2;
    bool hit1 = false, hit2 = false;
    bool hitL = b->l->bb.intersect(ray, lt1, lt2);
    bool hitR = b->r->bb.intersect(ray, rt1, rt2);
    if (lt1 > i->t && rt1 > i->t) return false;
    if (!(hitL || hitR)) return false;
    if (lt1 < rt1) {
      hit1 = bvhIntersect(b->l, ray, prims, i);
      if (i->t > rt1 && hitR) hit2 = bvhIntersect(b->r, ray, prims, i);
      return hit1 || hit2;
    } else {
      hit1 = bvhIntersect(b->r, ray, prims, i);
      if (i->t > lt1 && hitL) hit2 = bvhIntersect(b->l, ray, prims, i);
      return hit1 || hit2;
    }

  }
}

bool BVHAccel::intersect(const Ray &ray, Intersection *i) const {

  // TODO:
  // Implement ray - bvh aggregate intersection test. A ray intersects
  // with a BVH aggregate if and only if it intersects a primitive in
  // the BVH that is not an aggregate. When an intersection does happen.
  // You should store the non-aggregate primitive in the intersection data
  // and not the BVH aggregate itself.
  //
//  bool hit = bvhIntersect(root, ray, primitives, i);
  return bvhIntersect(root, ray, primitives, i);
}

}  // namespace StaticScene
}  // namespace CMU462

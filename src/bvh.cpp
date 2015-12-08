#include "bvh.h"

#include "CMU462/CMU462.h"
#include "static_scene/triangle.h"

#include <iostream>
#include <stack>

using namespace std;

namespace CMU462 { namespace StaticScene {


unsigned int morton3D(Vector3D);
typedef std::vector<Cluster> ClusterList;


float eps = 0.1;

/* Function must return a value between 0 and x. (lower -> quicker evaluation, lower quality BVH) */
float aacReduce(float x) {
  return pow(x, 0.5 - eps)*pow(4, eps);
}

double distance(Cluster A, Cluster B) {
  BBox tmp = BBox();
  tmp.expand(A->bb);
  tmp.expand(B->bb);
  return tmp.surface_area();
}

Cluster findBestMatch(std::list<std::pair<Cluster,Cluster>> C, Cluster target) {
  double bestD = INF_D;
  Cluster bestC = NULL;
  for (auto x = C.begin(); x != C.end(); x++) {
    if (x->first == target) continue;
    double dist = distance(x->first, target);
    if (dist < bestD) {
      bestD = dist;
      bestC = x->first;
    }
  }
  return bestC;
}


/* Combine all the clusters in clusters into n clusters */
ClusterList *BVHAccel::combineClusters(ClusterList *clusters, int n, int maxLeafSize) {
  if (clusters->size() <= n) {
    return clusters;
  }

  /* Consider using 2 lists and managing those (?) */
  std::list<std::pair<Cluster,Cluster>> C;
  C.resize(clusters->size());
  int m = clusters->size();
  std::vector<std::vector<double>> dists;
  dists.resize(m);
  for (int i = 0; i < m; i++) dists[i].resize(i);

  int i = 0;
  for (auto x = clusters->begin(); x != clusters->end(); x++) {
    int j = 0;
    for (auto y = clusters->begin(); y != x; y++) {
      dists[i][j] = distance(*x, *y);
      j++;
    }
    i++;
  }


  i = 0;
  auto cIter = C.begin();
  for (auto x = clusters->begin(); x != clusters->end(); x++) {
    int j = 0;
    double best = INF_D;
    for (auto y = clusters->begin(); y != clusters->end(); y++) {
      if (i != j) {
        double dist = dists[std::max(i,j)][std::min(i,j)];
        if (dist < best) {
          best = dist;
          *cIter = std::pair<Cluster,Cluster>(*x, *y);
        }
      }
      j++;
    }

    i++;
    cIter++;
  } 
 
 

  while (C.size() > n) {
    double best = INF_D;
    Cluster L, R;
    for (auto x = C.begin(); x != C.end(); x++) { 
      double dist = distance(x->first, x->second);
      if (dist < best) {
        best = dist;
        L = x->first; R = x->second;
      }
    }
    
    Cluster newC;
    if (L->count + R->count <= maxLeafSize) {
      LBVHLeaf *LL = (LBVHLeaf*)L;
      LBVHLeaf *RR = (LBVHLeaf*)R;
      LL->L->splice(LL->L->end(), *(RR->L));
      newC = new LBVHLeaf(LL->L);
    } else {
      newC = new LBVHParent(L, R);
    }

    auto P = [&](std::pair<Cluster,Cluster> x) { return x.first == L || x.first == R; };
    C.remove_if(P);
    C.push_back(std::pair<Cluster,Cluster>(newC, findBestMatch(C, newC)));
 
    for (auto x = C.begin(); x != C.end(); x++) {
      if (x->second == L || x->second == R) {
        x->second = findBestMatch(C, x->first);
      }
    }   
  }
 
  std::vector<Cluster> *out = new std::vector<Cluster>();
  for (auto x = C.begin(); x != C.end(); x++) {
    out->push_back(x->first);
  }
  return out; 
}

/* Returns a list of BVH trees that span primitives[start..end) */
ClusterList* BVHAccel::buildTree(std::vector<int>& M, int start, int end, int maxLeafSize, int bitPos, red_func f) {

  /* Base case */
  if (end - start <= maxLeafSize) {
    std::vector<Cluster> *C = new std::vector<Cluster>();
    for (int i = start; i < end; i++) {
      LBVHLeaf *tmp = new LBVHLeaf(new std::list<Primitive*>(1, primitives[i]));
      C->push_back(tmp);
    }
    return combineClusters(C, f(maxLeafSize), maxLeafSize);
  }
 
  /* Binary search for the split where M[i] & (1 << bitPos) becomes 1 */
  int split = (start + end) / 2;
  int mask = 1 << bitPos;
  int endTmp = end, startTmp = start;
  while (startTmp <= endTmp && bitPos >= 0) {
    if (split == start || split == end - 1) {
      split = (start + end) / 2;
      break;
    }
    else if (M[split] & mask) startTmp = split + 1;
    else endTmp = split - 1;
    split = (startTmp + endTmp) / 2;
  }


  ClusterList *pL, *pR;
  if (numThreads > 1) {
    /* Multithreaded tree building */
    std::thread tL, tR;
    bool threadL = false, threadR = false;

    auto runL = [&]{pL = buildTree(M, start, split, maxLeafSize, bitPos - 1, f); };
    auto runR = [&]{pR = buildTree(M, split, end, maxLeafSize, bitPos - 1, f); };

    int threads;

    lock.lock();
    if (runningThreads.load() < numThreads) {
      runningThreads.store(runningThreads+1);
      lock.unlock();
      threadL = true;
      tL = std::thread(runL);
    } else {
      lock.unlock();
      runL();
    }


    lock.lock();
    if (runningThreads.load() < numThreads) {
      runningThreads.store(runningThreads+1);
      lock.unlock();
      threadR = true;
      tR = std::thread(runR);
    } else {
      lock.unlock();
      runR();
    }

    /* Thread cleanup */
    if (threadL) {
      tL.join();
      lock.lock();
      runningThreads.store(runningThreads-1);
      lock.unlock();
    }
    if (threadR) {
      tR.join();
      lock.lock();
      runningThreads.store(runningThreads-1);
      lock.unlock();
    }
  } else {
    /* Single-threaded */
    pL = buildTree(M, start, split, maxLeafSize, bitPos - 1, f);
    pR = buildTree(M, split, end, maxLeafSize, bitPos - 1, f);
  }

  /* Combine lists of clusters */
  pL->insert(pL->end(), pR->begin(), pR->end());

  return combineClusters(pL, f(end - start), maxLeafSize);
}

/* LBVH Initializers */

LBVHLeaf::LBVHLeaf(std::list<Primitive *> *L) {
  this->count = L->size();
  this->L = L;
  this->bb = BBox();
  for (auto x = L->begin(); x != L->end(); x++) {
    this->bb.expand((*x)->get_bbox());
  }
  this->isLeaf = true;
}

LBVHParent::LBVHParent(LBVHNode *L, LBVHNode *R) {
  this->l = L;
  this->r = R;
  this->bb = BBox();
  this->bb.expand(L->bb);
  this->bb.expand(R->bb);
  this->isLeaf = false;
  this->count = L->count + R->count;
}



/* End LBVH Initializers */



/* Stuff for highest-level BVH function (AAC in paper) */


bool pair_cmp(MortonPt &p1, MortonPt &p2) {
  return p1.first < p2.first;
}
  

Cluster BVHAccel::makeBvhAAC(int maxLeafSize, red_func f) {
  std::vector<Vector3D> centroids;
  int n = primitives.size();
  centroids.resize(n);
  BBox total = BBox();
  for (int i = 0; i < n; i++) {
    BBox tmp = primitives[i]->get_bbox();
    total.expand(tmp);
    centroids[i] = (tmp.max + tmp.min) / 2.0;
  }
  Vector3D span = total.max - total.min;
  Vector3D span_inv = Vector3D(1.0 / span[0], 1.0 / span[1], 1.0 / span[2]);
  Vector3D adj = Vector3D(total.min[0] * span_inv[0],
                          total.min[1] * span_inv[1],
                          total.min[2] * span_inv[2]);
  for (int i = 0; i < n; i++) {
    centroids[i][0] = centroids[i][0] * span_inv[0] - adj[0];
    centroids[i][1] = centroids[i][1] * span_inv[1] - adj[1];
    centroids[i][2] = centroids[i][2] * span_inv[2] - adj[2];
  }
  std::vector<MortonPt> cent_zs;
  cent_zs.resize(n);
  for (int i = 0; i < n; i++) {
    cent_zs[i].first = morton3D(centroids[i]);
    cent_zs[i].second = primitives[i];
  }

  std::sort(cent_zs.begin(), cent_zs.end(), pair_cmp);

  std::vector<int> M;
  M.resize(cent_zs.size());
  int bitPos = 256;
  for (int i = 0; i < cent_zs.size(); i++) {
    M[i] = cent_zs[i].first;
    if (log2(M[i]) < bitPos) bitPos = log2(M[i]);
    primitives[i] = cent_zs[i].second;
  }
  
  std::vector<Cluster> *C = buildTree(M, 0, cent_zs.size(), maxLeafSize, bitPos, f);

  ClusterList *out = combineClusters(C, 1, maxLeafSize);
  return (*out)[0];
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

BVHNode *BVHAccel::rebuildBVH_single(Cluster c, int start) {
  if (c->isLeaf) {
    LBVHLeaf *L = (LBVHLeaf*)c;
    int i = start;
    for (auto x = L->L->begin(); x != L->L->end(); x++) {
      primitives[i++] = (*x);
    }
    return new BVHNode(c->bb, start, c->count);
  } else {
    LBVHParent *P = (LBVHParent*)c;
    BVHNode *root = new BVHNode(c->bb, start, c->count);
    root->l = rebuildBVH_single(P->l, start);
    root->r = rebuildBVH_single(P->r, start + P->l->count);
    return root;
  }
}

BVHNode *BVHAccel::rebuildBVH_threads(Cluster c, int start) {
  std::thread tL, tR;
  bool threadL = false, threadR = false;
  if (c->isLeaf) {
    LBVHLeaf *L = (LBVHLeaf*)c;
    int i = start;
    for (auto x = L->L->begin(); x != L->L->end(); x++) {
      primitives[i++] = (*x);
    }
    return new BVHNode(c->bb, start, c->count);
  } else {
    LBVHParent *P = (LBVHParent*)c;
    BVHNode *root = new BVHNode(c->bb, start, c->count);
    lock.lock();
    if (runningThreads.load() < numThreads) {
      runningThreads.store(runningThreads+1);
      lock.unlock();
      threadL = true;
      tL = std::thread([&]{root->l = rebuildBVH_threads(P->l, start);});
    } else {
      lock.unlock();
      root->l = rebuildBVH_threads(P->l, start);
    }
    lock.lock();
    if (runningThreads.load() < numThreads) {
      runningThreads.store(runningThreads+1);
      lock.unlock();
      threadR = true;
      tR = std::thread([&]{root->r = rebuildBVH_threads(P->r, start + P->l->count);});
    } else {
      lock.unlock();
      root->r = rebuildBVH_threads(P->r, start + P->l->count);
    }
    if (threadL) {
      tL.join();
      lock.lock();
      runningThreads.store(runningThreads-1);
      lock.unlock();
    }
    if (threadR) {
      tR.join();
      lock.lock();
      runningThreads.store(runningThreads-1);
      lock.unlock();
    }
    return root;
  }
}

BVHAccel::BVHAccel(const std::vector<Primitive *> &_primitives,
                   size_t max_leaf_size, size_t num_threads) {

  this->primitives = _primitives;
  this->numThreads = num_threads;
  this->runningThreads = 1;
  this->lock.unlock();
  if (primitives.size() < max_leaf_size) {
    BBox bb;
    for (size_t i = 0; i < primitives.size(); i++) {
      bb.expand(primitives[i]->get_bbox());
    }
    root = new BVHNode(bb, 0, primitives.size());
  } else {
    Cluster c = makeBvhAAC(max_leaf_size, aacReduce);
    if (numThreads > 1) {
      root = rebuildBVH_threads(c, 0);
    } else {
      root = rebuildBVH_single(c, 0);
    }
  }
}

void destroyBVHNode(BVHNode *node) {
  if (node->l) destroyBVHNode(node->l);
  if (node->r) destroyBVHNode(node->r);
  delete node;
}

BVHAccel::~BVHAccel() {
  destroyBVHNode(this->root);
}

BBox BVHAccel::get_bbox() const {
  return root->bb;
}

bool BVHAccel::intersect(const Ray &ray) const {
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

  return intersectNode(root, ray, i);
}

}  // namespace StaticScene
}  // namespace CMU462

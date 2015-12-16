/*
 * Student solution for CMU 15-462 Project 2 (MeshEdit)
 *
 * Implemented by Chris Kaffine on October 3rd, 2015.
 *
 */

#include "student_code.h"
#include "mutablePriorityQueue.h"

namespace CMU462
{
   VertexIter HalfedgeMesh::splitEdge( EdgeIter e0 )
   {
      //check if edge is on boundary
      if (e0->isBoundary()) return VertexIter();

      //assign relevent pointers
      HalfedgeIter he0, he1, ca, ab, bd, dc;
      //he0 begins as b->c and he1 begins as c->b, e0 is the edge bc
      he0 = e0->halfedge();
      ca = he0->next();
      ab = ca->next();
      he1 = he0->twin();
      bd = he1->next();
      dc = bd->next();
      VertexIter a, b, c, d;
      a = ab->vertex();
      b = he0->vertex();
      c = he1->vertex();
      d = dc->vertex();
      FaceIter f0, f1;
      f0 = he0->face();
      f1 = he1->face();

      //create new elements
      HalfedgeIter mc, mb, am, ma, md, dm;
      //he0 becomes b->m and he1 becoms c->m
      mc = newHalfedge();
      mb = newHalfedge();
      am = newHalfedge();
      ma = newHalfedge();
      md = newHalfedge();
      dm = newHalfedge();
      VertexIter m = newVertex();
      m->position = (c->position + b->position) / 2;
      m->velocity = (c->velocity + b->velocity) / 2;
      EdgeIter eMD, eMA, eMB; //e0 will be used as eMC
      eMD = newEdge();
      eMA = newEdge();
      eMB = newEdge();
      FaceIter f2, f3;
      f2 = newFace();
      f3 = newFace();

      //adjust neighbors of halfedges
      mc->setNeighbors(ca, he1, m, e0, f0);
      he1->setNeighbors(md, mc, c, e0, f1);
      mb->setNeighbors(bd, he0, m, eMB, f3);
      he0->setNeighbors(ma, mb, b, eMB, f2);
      am->setNeighbors(mc, ma, a, eMA, f0);
      ma->setNeighbors(ab, am, m, eMA, f2);
      md->setNeighbors(dc, dm, m, eMD, f1);
      dm->setNeighbors(mb, md, d, eMD, f3);
      ab->setNeighbors(he0, ab->twin(), ab->vertex(), ab->edge(), f2);
      bd->setNeighbors(dm, bd->twin(), bd->vertex(), bd->edge(), f3);
      dc->setNeighbors(he1, dc->twin(), dc->vertex(), dc->edge(), f1);
      ca->setNeighbors(am, ca->twin(), ca->vertex(), ca->edge(), f0);

      //adjust halfedges of other elements
      e0->halfedge() = he1;
      eMD->halfedge() = dm;
      eMA->halfedge() = am;
      eMB->halfedge() = he0;
      f0->halfedge() = ca;
      f1->halfedge() = dc;
      f2->halfedge() = ab;
      f3->halfedge() = bd;
      m->halfedge() = mc;
			return m;
	 }

   VertexIter HalfedgeMesh::collapseEdge( EdgeIter e )
   {
      //do not process boundary edges
      if (e->isBoundary()) return VertexIter();

      //assign relevant pointers
      HalfedgeIter ad, da, ac, ca, bc, cb, bd, db, cd, dc;
      cd = e->halfedge();
      dc = cd->twin();
      da = cd->next();
      ad = da->twin();
      ac = da->next();
      ca = ac->twin();
      cb = dc->next();
      bc = cb->twin();
      bd = cb->next();
      db = bd->twin();
      VertexIter a, b, c, d;
      a = ad->vertex();
      b = bc->vertex();
      c = cd->vertex();
      d = dc->vertex();
      EdgeIter eAD, eBD;
      eAD = ad->edge();
      eBD = bd->edge();
      FaceIter f0, f1;
      f0 = cd->face();
      f1 = dc->face();

      //check if there exists a three edge loop connecting c and without
      //passing through a or b. If so, the collapse will pinch the mesh,
      //making it non-manifold
      HalfedgeIter he0 = c->halfedge();
      do {
        //check all vertices connected to c
        VertexIter v = he0->twin()->vertex();
        if (v !=a && v!= b && v!=d) {
          HalfedgeIter he1 = v->halfedge();
          do {
            //check all vertices connected to current vertex, check for d
            if (he1->twin()->vertex() == d) { //found a loop, stop processing
              return VertexIter();
            }
            he1 = he1->twin()->next();
          } while (he1 != v->halfedge());
        }
        he0 = he0->twin()->next();
      } while (he0 != c->halfedge());

      //check degree of a and b, make sure it won't drop below 3
      if (a->degree() <= 3 || b->degree() <= 3) {
        return VertexIter();
      }
      //check resulting degree of m, make sure it's at least 3
      if (c->degree() + d->degree() - 4 < 3) {
        return VertexIter();
      }

      //adjust halfedges of surviving elements
      c->position = (c->position + d->position) / 2;
      c->velocity = (c->velocity + d->velocity) / 2;
      c->halfedge() = ca;
      a->halfedge() = ac;
      b->halfedge() = bc;
      ad->face()->halfedge() = ac;
      db->face()->halfedge() = cb;

      //update surviving halfedges
      ac->setNeighbors(ad->next(), ac->twin(), a, ac->edge(), ad->face());
      cb->setNeighbors(db->next(), cb->twin(), c, cb->edge(), db->face());
      //change all halfedges exiting d to have a vertex at c
      for (HalfedgeIter he = ad->next(); he != db; he = he->twin()->next())
        he->vertex() = c;
      ad->next()->next()->next() = ac;
      db->next()->next()->next() = cb;

      //perform deletions
      deleteHalfedge(ad);
      deleteHalfedge(da);
      deleteHalfedge(cd);
      deleteHalfedge(dc);
      deleteHalfedge(db);
      deleteHalfedge(bd);
      deleteVertex(d);
      deleteEdge(e);
      deleteEdge(eAD);
      deleteEdge(eBD);
      deleteFace(f0);
      deleteFace(f1);

			return c;
   }

   EdgeIter HalfedgeMesh::flipEdge( EdgeIter e0 )
   {
      //check if on boundary
      if (e0->isBoundary()) return e0;

      //assign relevant pointers
      HalfedgeIter he0, he1, ca, ab, bd, dc;
      he0 = e0->halfedge();
      ca = he0->next();
      ab = ca->next();
      he1 = he0->twin();
      bd = he1->next();
      dc = bd->next();
      VertexIter a, b, c, d;
      a = ab->vertex();
      b = he0->vertex();
      c = he1->vertex();
      d = dc->vertex();
      FaceIter f0, f1;
      f0 = he0->face();
      f1 = he1->face();

      //check if flipped version of edge already exists
      HalfedgeIter he = a->halfedge();
      do {
        if (he->twin()->vertex() == d) {
          return e0;
        }
        he = he->twin()->next();
      } while (he != a->halfedge());

      //update affected halfedges
      ca->setNeighbors(he0, ca->twin(), ca->vertex(), ca->edge(), f0);
      ab->setNeighbors(bd, ab->twin(), ab->vertex(), ab->edge(), f1);
      bd->setNeighbors(he1, bd->twin(), bd->vertex(), bd->edge(), f1);
      dc->setNeighbors(ca, dc->twin(), dc->vertex(), dc->edge(), f0);
      he0->setNeighbors(dc, he1, a, e0, f0);
      he1->setNeighbors(ab, he0, d, e0, f1);

      //update halfedges of affected elements
      c->halfedge() = ca;
      b->halfedge() = bd;
      f0->halfedge() = he0;
      f1->halfedge() = he1;

			return EdgeIter();
   }

   //returns true if either vertex on the edge has been marked as new,
   //indicating that this is the result of an edge split
   inline bool alreadySplit(EdgeIter& e)
   {
     return e->halfedge()->vertex()->isNew ||
            e->halfedge()->twin()->vertex()->isNew;
   }

   //implementation of abs, STL version seems to be causing problems
   int gabs(int x)
   {
     return (x <  0) ? -x : x;
   }

   bool isNLD(EdgeIter e) {
     Vector3D a, b, p, q; //p and q are on edge, a and b are opposite edge
     p = e->halfedge()->vertex()->position;
     q = e->halfedge()->twin()->vertex()->position;
     a = e->halfedge()->next()->twin()->vertex()->position;
     b = e->halfedge()->twin()->next()->twin()->vertex()->position;
     double cos1 = dot(p - a, q - a) / (p - a).norm() / (q - a).norm();
     double cos2 = dot(p - b, q - b) / (p - b).norm() / (q - b).norm();
     return acos(cos1) + acos(cos2) > PI;
   }

   void MeshResampler::resample( HalfedgeMesh& mesh )
   {
     bool done = true;
     do {
       done = true;
       for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
         if (e->length() < .000005) {
           if (mesh.collapseEdge(e) != VertexIter()) {
             done = false;
             break;
           }
         }
       }
     } while (!done);
     int newQueueSize = -1, prevQueueSize = -1;
     do {
       prevQueueSize = newQueueSize;
       MutablePriorityQueue<EdgeRecord> queue;
       int queueSize = 0;
       for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
         if (isNLD(e)) {
           EdgeRecord record;
           record.edge = e;
           Vector3D a, b, p, q; //p and q are on edge, a and b are opposite edge
           p = e->halfedge()->vertex()->position;
           q = e->halfedge()->twin()->vertex()->position;
           a = e->halfedge()->next()->twin()->vertex()->position;
           b = e->halfedge()->twin()->next()->twin()->vertex()->position;
           double cos1 = dot(p - a, q - a) / (p - a).norm() / (q - a).norm();
           double cos2 = dot(p - b, q - b) / (p - b).norm() / (q - b).norm();
           record.score = -(acos(cos1) + acos(cos2) - PI);
           e->record = record;
           queue.insert(record);
           queueSize++;
         }
       }
       newQueueSize = queueSize;

       while (queueSize > 0) {
         //get best edge from queue
         EdgeRecord nextEdge = queue.top();
         queue.pop();
         queueSize--;
         EdgeIter e = nextEdge.edge;
         mesh.flipEdge(e);
       }

     } while (newQueueSize != prevQueueSize);


     for (VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++)
       v->isNew = false;
     for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++)
       e->isNew = false;
     bool delauney;
     do {
       int numNLD = 0;
       delauney = true;
       std::vector<EdgeIter> issues;
       for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
         if (isNLD(e)) {
           numNLD++;
           delauney = false;
           issues.push_back(e);
         }
       }
       for (int i = 0; i < issues.size(); i++) {
         EdgeIter e = issues[i];
         if (e->isNew) mesh.flipEdge(e);
         else {
           Vector3D p, q;
           double velP, velQ;
           if (!e->halfedge()->vertex()->isNew) {
             p = e->halfedge()->vertex()->position;
             velP = e->halfedge()->vertex()->velocity;
             q = e->halfedge()->twin()->vertex()->position;
             velQ = e->halfedge()->twin()->vertex()->velocity;
           }
           else {
             p = e->halfedge()->twin()->vertex()->position;
             velP = e->halfedge()->twin()->vertex()->velocity;
             q = e->halfedge()->vertex()->position;
             velQ = e->halfedge()->vertex()->velocity;
           }
           VertexIter v = mesh.splitEdge(e);
           double d = (v->position - p).norm();
           double higherPower = pow(2, ceil(log2(d)));
           double lowerPower = pow(2, floor(log2(d)));
           double alpha = (higherPower - d < d - lowerPower) ? higherPower :
                                                               lowerPower;
           v->position = p + alpha * (v->position - p).unit();
           v->velocity = velP + (alpha / (p - q).norm()) * (velQ - velP);
           v->isNew = true;
           HalfedgeIter he = v->halfedge();
           he->edge()->isNew = false;
           he = he->twin()->next();
           he->edge()->isNew = true;
           he = he->twin()->next();
           he->edge()->isNew = false;
           he = he->twin()->next();
           he->edge()->isNew = true;
         }
       }
     } while (!delauney);
   }


   // Given an edge, the constructor for EdgeRecord finds the
   // optimal point associated with the edge's current quadric,
   // and assigns this edge a cost based on how much quadric
   // error is observed at this optimal point.
   EdgeRecord::EdgeRecord( EdgeIter& _edge )
   : edge( _edge )
   {
      //compute edge quadric
      Matrix4x4 K = edge->halfedge()->vertex()->quadric +
                    edge->halfedge()->twin()->vertex()->quadric;

      //solve for the point minimizing the quadric error
      Matrix3x3 A;
      Vector3D b;
      //construct system from K
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          A[i][j] = K[i][j];
        }
        b[i] = K[i][3];
      }
      optimalPoint = -A.inv() * b;

      //store cost associated with this point in the edge record
      Vector4D x(optimalPoint);
      x[3] = 1;
      score = dot(x, K * x);
   }

   void MeshResampler::downsample( HalfedgeMesh& mesh )
   {
     //Calculate quadric for each face
      for (FaceIter f = mesh.facesBegin(); f != mesh.facesEnd(); f++) {
        Vector3D p0, p1, p2;
        p0 = f->halfedge()->vertex()->position;
        p1 = f->halfedge()->next()->vertex()->position;
        p2 = f->halfedge()->next()->next()->vertex()->position;
        Vector3D n = cross(p1 - p0, p2 - p0);
        n.normalize();
        double d = -dot(n, p0);
        Vector4D v(n);
        v[3] = d;
        f->quadric = outer(v, v);
      }


      //get vertex quadrics by summing quadrics of connected faces
      for (VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++) {
        Matrix4x4 K;
        HalfedgeIter he = v->halfedge();
        do {
          K += he->face()->quadric;
          he = he->twin()->next();
        } while (he != v->halfedge());
        v->quadric = K;
      }

      //create EdgeRecord for each edge and place it in the priority queue
      MutablePriorityQueue<EdgeRecord> queue;
      for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
        EdgeRecord record(e);
        e->record = record;
        queue.insert(record);
      }

      //Keep collapsing the best edge until the number of edges in the mesh
      //has been reduced by a factor of 4
      double targetEdges = 4000;
      int queueSize = mesh.nEdges();
      while ((mesh.nEdges() > targetEdges) && (queueSize > 0)) {
        //get best edge from queue
        EdgeRecord nextEdge = queue.top();
        queue.pop();
        queueSize--;
        EdgeIter e = nextEdge.edge;
        //get quadric for vertex resulting from collapse
        Matrix4x4 K = e->halfedge()->vertex()->quadric +
                      e->halfedge()->twin()->vertex()->quadric;
        //remove all edges touching the current edge
        HalfedgeIter he = e->halfedge()->twin()->next();
        while (he != e->halfedge()) {
          queue.remove(he->edge()->record);
          queueSize--;
          he = he->twin()->next();
        }
        he = e->halfedge()->next();
        while (he != e->halfedge()->twin()) {
          queue.remove(he->edge()->record);
          queueSize--;
          he = he->twin()->next();
        }
        //attempt to collapse the edge
        VertexIter v = mesh.collapseEdge(e);
        if (v == VertexIter()) {
          //edge collapse failed, so add neighboring edges back into the queue
          he = e->halfedge()->twin()->next();
          while (he != e->halfedge()) {
            queue.insert(he->edge()->record);
            queueSize++;
            he = he->twin()->next();
          }
        }
        else {
          //edge collapse succeeded, so adjust new vertex and connected edges
          v->quadric = K;
          v->position = nextEdge.optimalPoint;
          he = v->halfedge();
          do {
            //update EdgeRecord and reinsert into the queue
            EdgeRecord record(he->edge());
            he->edge()->record = record;
            queue.insert(record);
            queueSize++;
            he = he->twin()->next();
          } while (he != v->halfedge());
        }
      }
   }

}

#ifndef STUDENT_CODE_H
#define STUDENT_CODE_H

#include "halfEdgeMesh.h"

using namespace std;

namespace CMU462 {

   class MeshResampler{

      public:

         MeshResampler(){};
         ~MeshResampler(){}

         void resample  ( HalfedgeMesh& mesh );
         void downsample ( HalfedgeMesh& mesh);
   };
}

#endif // STUDENT_CODE_H

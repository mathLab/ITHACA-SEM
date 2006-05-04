
#include <MultiRegions/ExpList2D.h>

namespace MultiRegions{
  ExpList2D::ExpList2D(){
  }
  
  ExpList2D::~ExpList2D(){
    vector<LocalRegions::QuadExp*>::iterator defQ;
    for(defQ = _Quad.begin(); defQ != _Quad.end(); ++defQ)
      delete defQ[0];

    vector<LocalRegions::TriExp*>::iterator defT;
    for(defT = _Tri.begin(); defT != _Tri.end(); ++defT)
      delete defT[0];
  }


} //end of namespace


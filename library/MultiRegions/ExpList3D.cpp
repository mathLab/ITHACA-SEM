
#include <MultiRegions/ExpList3D.h>

namespace MultiRegions{
  ExpList3D::ExpList3D(){
  }
  
  ExpList3D::~ExpList3D(){
    vector<LocalRegions::HexExp*>::iterator defH;
    for(defH = _Hex.begin(); defH != _Hex.end(); ++defH)
      delete defH[0];

    vector<LocalRegions::PrismExp*>::iterator defPr;
    for(defPr = _Prism.begin(); defPr != _Prism.end(); ++defPr)
      delete defPr[0];

    vector<LocalRegions::PyrExp*>::iterator defPy;
    for(defPy = _Pyr.begin(); defPy != _Pyr.end(); ++defPy)
      delete defPy[0];

    vector<LocalRegions::TetExp*>::iterator defT;
    for(defT = _Tet.begin(); defT != _Tet.end(); ++defT)
      delete defT[0];
  }


} //end of namespace


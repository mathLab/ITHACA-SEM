#ifndef H_EXPLIST3D

#include <vector>
#include <MultiRegions/MultiRegions.h>
#include <MultiRegions/ExpList.h>

#include <LocalRegions/HexExp.h>
#include <LocalRegions/PrismExp.h>
#include <LocalRegions/PyrExp.h>
#include <LocalRegions/TetExp.h>

namespace MultiRegions{

  class ExpList3D: public ExpList{
  private:
    vector<LocalRegions::HexExp*>  _Hex;
    vector<LocalRegions::PrismExp*> _Prism;
    vector<LocalRegions::PyrExp*>   _Pyr;
    vector<LocalRegions::TetExp*>   _Tet;
    
  protected:


  public:
    ExpList3D();
    ~ExpList3D();
  };


} //end of namespace

#define H_EXPLIST3D
#endif

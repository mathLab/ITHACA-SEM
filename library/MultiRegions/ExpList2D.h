#ifndef H_EXPLIST2D

#include <vector>
#include <MultiRegions/MultiRegions.h>
#include <MultiRegions/ExpList.h>

#include <LocalRegions/TriExp.h>
#include <LocalRegions/QuadExp.h>

namespace MultiRegions{

  class ExpList2D: public ExpList{
  private:
    vector<LocalRegions::QuadExp*> _Quad;
    vector<LocalRegions::TriExp*> _Tri;
    
  protected:


  public:
    ExpList2D();
    ~ExpList2D();
  };


} //end of namespace

#define H_EXPLIST2D
#endif

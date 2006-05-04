#ifndef H_DISCONTEXPLIST1D

#include <MultiRegions/MultiRegions.h>
#include <MultiRegions/ExpList1D.h>

namespace MultiRegions{

  class DisContExpList1D: public ExpList1D {
  private:
  protected:
  public:
    DisContExpList1D();
    ~DisContExpList1D();
  };
} //end of namespace

#define H_DISCONTEXPLIST1D
#endif

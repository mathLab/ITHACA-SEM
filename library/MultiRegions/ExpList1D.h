#ifndef H_EXPLIST1D

#include <vector>
#include <fstream>
#include <MultiRegions/MultiRegions.h>
#include <MultiRegions/ExpList.h>

#include <StdRegions/StdBasis.h>

#include <LocalRegions/SegExp.h>

#include <SpatialDomains/MeshGraph1D.h>

namespace MultiRegions{

  class ExpList1D: public ExpList{
  private:

  protected:
    vector<LocalRegions::SegExp *> _Seg;


  public:
    ExpList1D();

    ExpList1D(const StdRegions::BasisKey &Ba, 
	      SpatialDomains::MeshGraph1D &graph1D);

    ~ExpList1D();

    double Integral(const double *inarray);
    void   IProduct_WRT_B(const double *inarray, double *outarray);
    void   IProduct_WRT_B(ExpList1D &S1, ExpList1D &S2);
    void   IProduct_WRT_B(ExpList1D &S1, double * outarray);
    void   Deriv    (const int n, double **outarray);
    void   Deriv    (const int n, const double *inarray, double ** outarray);
    void   FwdTrans (const double *inarray);
    void   BwdTrans (double *outarray); 

    virtual void V_BwdTrans(double *outarray){
      BwdTrans(outarray);
    }

    void   GetCoords(double **coords);
    void   WriteToFile(ofstream &out);
    
    inline int get_coordim(int eid){
      ASSERTL2(eid <= _Seg.size(),"ExpList1D:get_coordim()",
		"eid is larger than number of elements");
      return _Seg[eid]->get_coordim();
    }

    double Linf (const double *sol);
    double L2   (const double *sol);

  };
  
} //end of namespace

#define H_EXPLIST1D
#endif

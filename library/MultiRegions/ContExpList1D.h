#ifndef H_CONTEXPLIST1D

#include <MultiRegions/MultiRegions.h>
#include <MultiRegions/ExpList1D.h>

#include <StdRegions/StdMatrix.h>

namespace MultiRegions{

  class ContExpList1D: public ExpList1D {
  private:
    int    _cont_ncoeffs;
    int    *_LocToContMap;

    double *_cont_coeffs;
    
    StdRegions::StdMatContainer *_Mass;

  protected:

  public:
    ContExpList1D();
    ContExpList1D(const StdRegions::BasisKey &Ba, 
		  SpatialDomains::MeshGraph1D &graph1D);
    ~ContExpList1D();

    inline int get_cont_ncoeffs(){
      return _cont_ncoeffs;
    }

    inline double *get_cont_coeffs(){
      return _cont_coeffs;
    }
    
    inline void ContToLocal(){
      ContToLocal(_cont_coeffs,_coeffs);
    }
    inline void ContToLocal(const double *cont,double *loc){
      Vmath::gathr(_ncoeffs,cont,_LocToContMap,loc);
    }

    inline void LocalToCont(){
      LocalToCont(_coeffs,_cont_coeffs);
    }

    inline void LocalToCont(const double *loc, double *cont){
      Vmath::scatr(_ncoeffs,loc,_LocToContMap,cont);
    }

    inline void Assemble(){
      Assemble(_coeffs,_cont_coeffs);
    }

    inline void Assemble(const double *loc, double *cont){
      Vmath::zero(_cont_ncoeffs,cont,1);
      Vmath::assmb(_ncoeffs,loc,_LocToContMap,cont);
    }
		       
    void IProduct_WRT_B(const double *inarray, double *outarray);

    void FwdTrans(const double *inarray);

    void BwdTrans(double *outarray);

    virtual void V_BwdTrans(double *outarray){
      BwdTrans(outarray);
    }
    void GenMassMatrix(void);
    
  };
} //end of namespace

#define H_CONTEXPLIST1D
#endif

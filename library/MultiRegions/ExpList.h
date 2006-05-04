#ifndef H_EXPLIST

#include <MultiRegions/MultiRegions.h>

namespace MultiRegions{

  class ExpList{
  private:

  protected:
    int _ncoeffs; // < 
    int _nquad;
    
    double *_coeffs;
    double *_phys;

    TransState _trans_state;
    bool       _phys_state;

  public:
     ExpList();
    ~ExpList();

    inline int get_ncoeffs(){
      return _ncoeffs;
    }

    inline int get_nquad(){
      return _nquad;
    }

    inline double *get_coeffs(){
      return _coeffs;
    }

    inline double *get_phys(){
      return _phys;
    }

    inline void set_trans_state(TransState trans_state){
      _trans_state = trans_state;
    }

    inline void set_phys_state(bool phys_state){
      _phys_state = phys_state;
    }

    //virtuals
    virtual void   V_BwdTrans (double *outarray)  = 0;
  };

} //end of namespace

#define H_EXPLIST
#endif

#include <MultiRegions/ContExpList1D.h>

namespace MultiRegions{

  ContExpList1D::ContExpList1D(){
    _Mass = (StdRegions::StdMatContainer *)NULL;
    _LocToContMap = (int *)NULL;
  }
  
  ContExpList1D::~ContExpList1D(){
    if(_cont_coeffs)
      delete[] _cont_coeffs;

    if(_LocToContMap)
      delete[] _LocToContMap;
  }

  ContExpList1D::ContExpList1D(const StdRegions::BasisKey &Ba, 
     SpatialDomains::MeshGraph1D &graph1D):ExpList1D(Ba,graph1D){
    int gid,i,j;
    int order = Ba.get_order();

    ASSERTL1((Ba.get_Btype() == StdRegions::Modified_A)
	     ||(Ba.get_Btype() == StdRegions::GLL_Lagrange),
	     "ConExpList1D::ConExpList1D",
	     "Expansion not of an boundary-interior type");

    _Mass = (StdRegions::StdMatContainer *)NULL;

    _LocToContMap = new int [_ncoeffs];
    Vmath::fill(_ncoeffs,-1,_LocToContMap,1);
    
    // set up mapping based 
    StdRegions::StdExpMap vmap;
    
    // assume all elements have the same mapping and expasion order
    _Seg[0]->MapTo(StdRegions::Forwards, vmap);
    
    // set up simple map;
    for(gid = i = 0; i < _Seg.size(); ++i,++gid)
      for(j = 0; j < 2; ++j)
	_LocToContMap[order*i+vmap[j]] = gid+j;
    ++gid;

    for(i = 0; i < _Seg.size(); ++i)
      for(j = 0; j < order; ++j)
	if(_LocToContMap[order*i+j] == -1)
	  _LocToContMap[order*i+j] = gid++;

    _cont_ncoeffs = gid;
    _cont_coeffs = new double [_cont_ncoeffs];
  }
  

  void ContExpList1D::IProduct_WRT_B(const double *inarray, double *outarray){
    ExpList1D::IProduct_WRT_B(inarray,_coeffs);
    Assemble(_coeffs,outarray);
    _trans_state = LocalCont;
  }

  void ContExpList1D::FwdTrans(const double *inarray){
    IProduct_WRT_B(inarray,_cont_coeffs);
    if(!_Mass) GenMassMatrix();
    _Mass->Solve(_cont_coeffs,1);
    _trans_state = Continuous;
    _phys_state = false;
  }

  void ContExpList1D::BwdTrans(double *outarray){

    if(_trans_state == Continuous)
      ContToLocal();

    ExpList1D::BwdTrans(outarray);
  }


  void ContExpList1D::GenMassMatrix(void){
    if(!_Mass){
      int   i,j,cnt,gid1,gid2,loc_lda;
      double *loc_mat;
      StdRegions::StdMatContainer *loc_mass;
      vector<LocalRegions::SegExp *>::iterator def;

      double *mmat = new double [_cont_ncoeffs*_cont_ncoeffs];
      Vmath::zero(_cont_ncoeffs*_cont_ncoeffs,mmat,1);
      
      _Mass = new StdRegions::StdMatContainer(mmat);
      _Mass->set_lda      (_cont_ncoeffs);
      _Mass->set_mat_form (StdRegions::Symmetric_Positive);

      // fill global matrix 
      for(cnt = 0, def = _Seg.begin(); def != _Seg.end(); ++def){
	loc_mass = (*def)->GetMassMatrix();
	loc_lda = loc_mass->get_lda();
	loc_mat = loc_mass->get_matrix();
	
	for(i = 0; i < loc_lda; ++i){
	  gid1 = _LocToContMap[i+cnt];
	  for(j = 0; j < loc_lda; ++j){
	    gid2 = _LocToContMap[j+cnt];
	    mmat[gid1*_cont_ncoeffs + gid2] += loc_mat[i*loc_lda + j];
	  }
	}
	cnt+=(*def)->get_ncoeffs();
      }
    }
  }
} //end of namespace


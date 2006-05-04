#include <MultiRegions/ExpList1D.h>

namespace MultiRegions{

  ExpList1D::ExpList1D(){
  }
  
  ExpList1D::~ExpList1D(){
    vector<LocalRegions::SegExp*>::iterator def;
    for(def = _Seg.begin(); def != _Seg.end(); ++def){
      delete def[0]->get_geom();
      delete def[0];
    }
  }

  ExpList1D::ExpList1D(const StdRegions::BasisKey &Ba, 
		       SpatialDomains::MeshGraph1D &graph1D){
    LocalRegions::SegExp *seg;
    list<SpatialDomains::SegGeom *> SegGeoms = graph1D.get_SegGeoms();
    
    _ncoeffs = SegGeoms.size()*Ba.get_order();
    _nquad   = SegGeoms.size()*Ba.get_nquad();
    
    _coeffs = new double [_ncoeffs];
    _trans_state = NotSet; 

    _phys   = new double [_nquad];
    _phys_state  = false;
    
    // make sure Geofacs are defined in MeshGraph1D
    if(graph1D.get_geofac_defined() != true)
      graph1D.GenXGeoFac();

    list<SpatialDomains::SegGeom *>::const_iterator def;
    int cnt,cnt1;
    SpatialDomains::SegGeom *geom;

    cnt = cnt1 = 0;
    for(def = SegGeoms.begin(); def != SegGeoms.end(); ++def){
      geom = new SpatialDomains::SegGeom (**def);
      seg  = new LocalRegions::SegExp(Ba,_coeffs+cnt,_phys+cnt1, geom);
      seg->SetGeoFac(seg->GenGeoFac());
      _Seg.push_back(seg);
      
      cnt  += Ba.get_order();
      cnt1 += Ba.get_nquad();
    }
  }
  
  /** \brief Integrate the physical point list \a inarray over region
      and return the value
      
      Inputs:\n
      
      - \a inarray: definition of function to be returned at quadrature point 
      of expansion. 
      
      Outputs:\n

      - returns \f$ \sum_{i=1}^{n_{el}} \int_{\Omega_i} u(\xi_1)d \xi_1 \f$ 
  */
  double ExpList1D::Integral(const double *inarray){
    vector<LocalRegions::SegExp *>::iterator def;
    int    cnt = 0;
    double sum = 0.0;

    for(def = _Seg.begin(); def != _Seg.end(); ++def){
      sum += (*def)->Integral(inarray+cnt);
      cnt += (*def)->get_pointorder(0);
    }

    return sum; 
  }
  
  void ExpList1D::IProduct_WRT_B(const double *inarray, double *outarray){
    vector<LocalRegions::SegExp *>::iterator def;
    int    cnt  = 0;
    int    cnt1 = 0;
    
    for(def = _Seg.begin(); def != _Seg.end(); ++def){
      (*def)->IProduct_WRT_B(inarray+cnt,outarray+cnt1);
      cnt  += (*def)->get_pointorder(0);
      cnt1 += (*def)->get_ncoeffs();
    }
  }

  void ExpList1D::IProduct_WRT_B(ExpList1D &S1, ExpList1D &S2){
    IProduct_WRT_B(S1.get_phys(),S2.get_coeffs());
    _trans_state = Local;
  }

  void ExpList1D::IProduct_WRT_B(ExpList1D &S1, double * outarray){
    IProduct_WRT_B( S1.get_phys(),outarray);
  }

  void ExpList1D::Deriv(const int n, double **outarray){
    Deriv(n,_phys,outarray);
  }

  void ExpList1D::Deriv(const int n,const double *inarray,double **outarray){
    vector<LocalRegions::SegExp *>::iterator def;
    int    cnt = 0;
    
    if(_phys_state == false)
      V_BwdTrans(_phys);
    
    for(def = _Seg.begin(); def != _Seg.end(); ++def){
      (*def)->Deriv(n,inarray+cnt,outarray+cnt);
      cnt  += (*def)->get_pointorder(0);
    }
  }

  void ExpList1D::FwdTrans(const double *inarray){
    vector<LocalRegions::SegExp *>::iterator def;
    int    cnt = 0;

    for(def = _Seg.begin(); def != _Seg.end(); ++def){
      (*def)->FwdTrans(inarray+cnt);
      cnt  += (*def)->get_pointorder(0);
    }

    _trans_state = Local;
  }

  void ExpList1D::BwdTrans(double *outarray){
    vector<LocalRegions::SegExp *>::iterator def;
    int    cnt = 0;

    for(def = _Seg.begin(); def != _Seg.end(); ++def){
      (*def)->BwdTrans(outarray+cnt);
      cnt  += (*def)->get_pointorder(0);
    }
    _phys_state = true;
  }

  void ExpList1D::GetCoords(double **coords){
    vector<LocalRegions::SegExp *>::iterator def;
    int    i, cnt = 0;
    double *E_coords[3];

    for(def = _Seg.begin(); def != _Seg.end(); ++def){
      for(i = 0 ; i < (*def)->get_coordim(); ++i)
	E_coords[i] = coords[i]+cnt;
      
      (*def)->GetCoords(E_coords);
      cnt  += (*def)->get_pointorder(0);
    }
  }
  
  void ExpList1D::WriteToFile(ofstream &out){
    vector<LocalRegions::SegExp *>::iterator def; 

    if(_phys_state == false)
      V_BwdTrans(_phys);

    (*_Seg.begin())->WriteToFile(out,1);
    
    for(def = ++_Seg.begin(); def != _Seg.end(); ++def)
      (*def)->WriteToFile(out,0);
  }

  double  ExpList1D::Linf(const double *sol){
    vector<LocalRegions::SegExp *>::iterator def;
    double err = 0.0;
    int    cnt = 0;

    if(_phys_state == false)
      V_BwdTrans(_phys);
    
    for(def = _Seg.begin(); def != _Seg.end(); ++def){
      err  = max(err,(*def)->Linf(sol+cnt));
      cnt  += (*def)->get_pointorder(0);
    }

    return err;
  }

  double  ExpList1D::L2(const double *sol){
    vector<LocalRegions::SegExp *>::iterator def;
    double err = 0.0,errl2;
    int    cnt = 0;
    
    if(_phys_state == false)
      V_BwdTrans(_phys);

    for(def = _Seg.begin(); def != _Seg.end(); ++def){
      errl2 = (*def)->L2(sol+cnt);
      err += errl2*errl2;
      cnt  += (*def)->get_pointorder(0);
    }

    return sqrt(err);
  }
  

} //end of namespace


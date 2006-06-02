///////////////////////////////////////////////////////////////////////////////
//
// File ExpList1D.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: Expansion list 1D definition
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/ExpList1D.h>

namespace Nektar
{
    namespace MultiRegions
    {
	
	ExpList1D::ExpList1D()
	{
	}
	
	ExpList1D::~ExpList1D()
	{
	}
	
	ExpList1D::ExpList1D(const StdRegions::BasisKey &Ba, 
			     SpatialDomains::MeshGraph1D &graph1D)
	{
	    LocalRegions::SegExpSharedPtr seg;
	    SpatialDomains::SegGeomVector SegGeoms = graph1D.GetSeggeoms();
	    
	    m_ncoeffs = SegGeoms.size()*Ba.GetBasisOrder();
	    m_npoints = SegGeoms.size()*Ba.GetPointsOrder();
	    
	    m_coeffs = new double [m_ncoeffs];
	    m_transState = eNotSet; 
	    
	    m_phys   = new double [m_npoints];
	    m_physState  = false;
	    
	    // make sure Geofacs are defined in MeshGraph1D
	    if(graph1D.GetGeofac_defined() != true)
	    {
		graph1D.GenXGeoFac();
	    }
	    
	    SpatialDomains::SegGeomVectorIter def;
	    //      SpatialDomains::SegGeomSharedPtr geom;
	    int cnt,cnt1;
	    
	    cnt = cnt1 = 0;
	    for(def = SegGeoms.begin(); def != SegGeoms.end(); ++def)
	    {
		// removed copy construction of geom
		// geom = new SpatialDomains::SegGeom (**def);
		seg.reset( new LocalRegions::SegExp(Ba,m_coeffs+cnt,m_phys+cnt1, 
						    *def));
		seg->SetGeoFac(seg->GenGeoFac());
		m_seg.push_back(seg);
		
		cnt  += Ba.GetBasisOrder();
		cnt1 += Ba.GetPointsOrder();
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
      double ExpList1D::Integral(const double *inarray)
      {
	  LocalRegions::SegExpVectorIter def;
	  int    cnt = 0;
	  double sum = 0.0;
	  
	  for(def = m_seg.begin(); def != m_seg.end(); ++def){
	      sum += (*def)->Integral(inarray+cnt);
	      cnt += (*def)->GetPointsTot();
	  }
	  
	  return sum; 
      }
      
  
      void ExpList1D::IProductWRTBase(const double *inarray, double *outarray)
      {
	  LocalRegions::SegExpVectorIter def;
	  int    cnt  = 0;
	  int    cnt1 = 0;
	  
	  for(def = m_seg.begin(); def != m_seg.end(); ++def)
	  {
	      (*def)->IProductWRTBase(inarray+cnt,outarray+cnt1);
	      cnt  += (*def)->GetPointsTot();
	      cnt1 += (*def)->GetNcoeffs();
	  }
      }
      
      void ExpList1D::IProductWRTBase(ExpList1D &S1, ExpList1D &S2)
      {
	  IProductWRTBase(S1.GetPhys(),S2.GetCoeffs());
	  m_transState = eLocal;
      }
      
      void ExpList1D::IProductWRTBase(ExpList1D &S1, double * outarray)
      {
	  IProductWRTBase( S1.GetPhys(),outarray);
      }
      
      void ExpList1D::Deriv(const int n, double **outarray)
      {
	  Deriv(n,m_phys,outarray);
      }
      
      void ExpList1D::Deriv(const int n, const double *inarray,
			    double **outarray)
      {
	  LocalRegions::SegExpVectorIter def;
	  int    cnt = 0;
	  
	  if(m_physState == false)
	  {
	      v_BwdTrans(m_phys);
	  }
	  
	  for(def = m_seg.begin(); def != m_seg.end(); ++def)
	  {
	      (*def)->Deriv(n,inarray+cnt,outarray+cnt);
	      cnt  += (*def)->GetPointsTot();
	  }
      }
      
      void ExpList1D::FwdTrans(const double *inarray)
      {
	  LocalRegions::SegExpVectorIter def;
	  int    cnt = 0;
	  
	  for(def = m_seg.begin(); def != m_seg.end(); ++def)
	  {
	      (*def)->FwdTrans(inarray+cnt);
	      cnt  += (*def)->GetPointsTot();
	  }
	
	  m_transState = eLocal;
      }
      
      void ExpList1D::BwdTrans(double *outarray)
      {
	  LocalRegions::SegExpVectorIter def;
	  int    cnt = 0;
	  
	  for(def = m_seg.begin(); def != m_seg.end(); ++def)
	  {
	      (*def)->BwdTrans(outarray+cnt);
	      cnt  += (*def)->GetPointsTot();
	  }
	  m_physState = true;
      }
      
      void ExpList1D::GetCoords(double **coords)
      {
	  LocalRegions::SegExpVectorIter def;
	  int    i, cnt = 0;
	  double *E_coords[3];
	  
	  for(def = m_seg.begin(); def != m_seg.end(); ++def)
	  {
	      for(i = 0 ; i < (*def)->GetCoordim(); ++i)
	      {
		  E_coords[i] = coords[i]+cnt;
	      }
	      
	      (*def)->GetCoords(E_coords);
	      cnt  += (*def)->GetPointsTot();
	  }
      }
      
      void ExpList1D::WriteToFile(std::ofstream &out)
      {
	  LocalRegions::SegExpVectorIter def; 
	  
	  if(m_physState == false)
	  {
	      v_BwdTrans(m_phys);
	  }
	  
	  (*m_seg.begin())->WriteToFile(out,1);
	  
	  for(def = ++m_seg.begin(); def != m_seg.end(); ++def)
	  {
	      (*def)->WriteToFile(out,0);
	  }
      }
      
      double  ExpList1D::Linf(const double *sol)
      {
	  LocalRegions::SegExpVectorIter def;
	  double err = 0.0;
	  int    cnt = 0;
	  
	  if(m_physState == false)
	  {
	      v_BwdTrans(m_phys);
	  }
	  
	  for(def = m_seg.begin(); def != m_seg.end(); ++def)
	  {
	      err  = std::max(err,(*def)->Linf(sol+cnt));
	      cnt  += (*def)->GetPointsTot();
	  }
	  
	  return err;
      }
      
      double  ExpList1D::L2(const double *sol)
      {
	  LocalRegions::SegExpVectorIter def;
	  double err = 0.0,errl2;
	  int    cnt = 0;
	  
	  if(m_physState == false)
	  {
	      v_BwdTrans(m_phys);
	  }
	  
	  for(def = m_seg.begin(); def != m_seg.end(); ++def)
	  {
	      errl2 = (*def)->L2(sol+cnt);
	      err += errl2*errl2;
	      cnt  += (*def)->GetPointsTot();
	  }
	  
	  return sqrt(err);
      }
      
    
  } //end of namespace
} //end of namespace


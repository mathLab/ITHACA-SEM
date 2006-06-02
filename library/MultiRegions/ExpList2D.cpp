///////////////////////////////////////////////////////////////////////////////
//
// File ExpList2D.cpp
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
// Description: Expansion list 2D definition
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/ExpList2D.h>

namespace Nektar
{
    namespace MultiRegions
    {
	
	ExpList2D::ExpList2D()
	{
	}
	
	ExpList2D::~ExpList2D()
	{
	}
	
	ExpList2D::ExpList2D(const StdRegions::BasisKey &TriBa, 
			     const StdRegions::BasisKey &TriBb, 
			     const StdRegions::BasisKey &QuadBa, 
			     const StdRegions::BasisKey &QuadBb, 
			     SpatialDomains::MeshGraph2D &graph2D)
	{
	    SpatialDomains::TriGeomVector TriGeoms   = graph2D.GetTrigeoms();
	    SpatialDomains::QuadGeomVector QuadGeoms = graph2D.GetQuadgeoms();
	    int cnt,cnt1;
	    int tri_ncoeffs_elmt = 0, quad_ncoeffs_elmt = 0;
	    int tri_npoints_elmt = 0, quad_npoints_elmt = 0;
	    
	    // declare memory

	    if(TriGeoms.size())
	    {		
		tri_ncoeffs_elmt = (TriBa.GetBasisOrder()+1)/2*
		    TriBa.GetBasisOrder()*(TriBb.GetBasisOrder()-
					   TriBa.GetBasisOrder());
		tri_npoints_elmt = (TriBa.GetPointsOrder()*
				    TriBb.GetPointsOrder());
		
		m_ncoeffs += TriGeoms.size()*tri_ncoeffs_elmt;
		m_npoints += TriGeoms.size()*tri_npoints_elmt;
	    }		
		
	    if(QuadGeoms.size())
	    {		
		quad_ncoeffs_elmt = QuadBa.GetBasisOrder() 
		    *QuadBb.GetBasisOrder();
		quad_npoints_elmt = QuadBa.GetPointsOrder()
		    *QuadBb.GetPointsOrder();
		
		m_ncoeffs += QuadGeoms.size()*quad_ncoeffs_elmt;
		m_npoints += QuadGeoms.size()*quad_npoints_elmt;
	    }
	    
	    m_coeffs = new double [m_ncoeffs];
	    m_transState = eNotSet; 
	    
	    m_phys   = new double [m_npoints];
	    m_physState  = false;
		
	    // declare triangles using first block of data 	    
	    if(TriGeoms.size())
	    {		
		LocalRegions::TriExpSharedPtr tri;
		SpatialDomains::TriGeomVectorIter def;

		// make sure Geofacs are defined in MeshGraph1D
		if(graph2D.GetGeofac_defined() != true)
		{
		    graph2D.GenXGeoFac();
		}
	    
		cnt = cnt1 = 0;
		for(def = TriGeoms.begin(); def != TriGeoms.end(); ++def)
		{
		    // removed copy construction of geom
		    // geom = new SpatialDomains::SegGeom (**def);
		    tri.reset(new LocalRegions::TriExp(TriBa,TriBb,
						       m_coeffs+cnt,m_phys+cnt1, 
						       *def));
		    tri->SetGeoFac(tri->GenGeoFac());
		    m_tri.push_back(tri);
		    
		    cnt  += tri_ncoeffs_elmt;
		    cnt1 += tri_npoints_elmt;
		}
	    }
	    
	    // set up quads 
	    if(QuadGeoms.size())
	    {		
		LocalRegions::QuadExpSharedPtr quad;
		SpatialDomains::QuadGeomVectorIter def;

		// make sure Geofacs are defined in MeshGraph1D
		if(graph2D.GetGeofac_defined() != true)
		{
		    graph2D.GenXGeoFac();
		}
	    
		cnt = cnt1 = 0;
		for(def = QuadGeoms.begin(); def != QuadGeoms.end(); ++def)
		{
		    // removed copy construction of geom
		    // geom = new SpatialDomains::SegGeom (**def);
		    quad.reset(new LocalRegions::QuadExp(QuadBa,QuadBb,
							m_coeffs+cnt,m_phys+cnt1, 
							*def));
		    quad->SetGeoFac(quad->GenGeoFac());
		    m_quad.push_back(quad);
		    
		    cnt  += quad_ncoeffs_elmt;
		    cnt1 += quad_npoints_elmt;
		}
	    }
	}
    
	double ExpList2D::Integral(const double *inarray)
	{
	    LocalRegions::TriExpVectorIter  Tdef;
	    LocalRegions::QuadExpVectorIter Qdef;
	    int    cnt = 0;
	    double sum = 0.0;
	    
	    for(Tdef = m_tri.begin(); Tdef != m_tri.end(); ++Tdef){
		sum += (*Tdef)->Integral(inarray+cnt);
		cnt += (*Tdef)->GetPointsTot();
	    }

	    for(Qdef = m_quad.begin(); Qdef != m_quad.end(); ++Qdef){
		sum += (*Qdef)->Integral(inarray+cnt);
		cnt += (*Qdef)->GetPointsTot();
	    }
	    
	    return sum; 
	}
	
  
	void ExpList2D::IProductWRTBase(const double *inarray, 
					double *outarray)
	{
	    LocalRegions::TriExpVectorIter  Tdef;
	    LocalRegions::QuadExpVectorIter Qdef;
	    int    cnt  = 0;
	    int    cnt1 = 0;
	    
	    for(Tdef = m_tri.begin(); Tdef != m_tri.end(); ++Tdef)
	    {
		(*Tdef)->IProductWRTBase(inarray+cnt,outarray+cnt1);
		cnt  += (*Tdef)->GetPointsTot();
		cnt1 += (*Tdef)->GetNcoeffs();
	    }

	    for(Qdef = m_quad.begin(); Qdef != m_quad.end(); ++Qdef)
	    {
		(*Qdef)->IProductWRTBase(inarray+cnt,outarray+cnt1);
		cnt  += (*Qdef)->GetPointsTot();
		cnt1 += (*Qdef)->GetNcoeffs();
	    }
	}
      
	void ExpList2D::IProductWRTBase(ExpList2D &S1, ExpList2D &S2)
	{
	    IProductWRTBase(S1.GetPhys(),S2.GetCoeffs());
	    m_transState = eLocal;
	}
      
	void ExpList2D::IProductWRTBase(ExpList2D &S1, double * outarray)
	{
	    IProductWRTBase( S1.GetPhys(),outarray);
	}
      
	void ExpList2D::Deriv(const int n, double **outarray)
	{
	    Deriv(n,m_phys,outarray);
	}
      
	void ExpList2D::Deriv(const int n, const double *inarray,
			      double **outarray)
	{
	    LocalRegions::TriExpVectorIter Tdef;
	    LocalRegions::QuadExpVectorIter Qdef;
	    int    cnt = 0;
	  
	    if(m_physState == false)
	    {
		v_BwdTrans(m_phys);
	    }
	  
	    for(Tdef = m_tri.begin(); Tdef != m_tri.end(); ++Tdef)
	    {
		(*Tdef)->Deriv(n,inarray+cnt,outarray+cnt);
		cnt  += (*Tdef)->GetPointsTot();
	    }

	    for(Qdef = m_quad.begin(); Qdef != m_quad.end(); ++Qdef)
	    {
		(*Qdef)->Deriv(n,inarray+cnt,outarray+cnt);
		cnt  += (*Qdef)->GetPointsTot();
	    }
	}
	
	void ExpList2D::FwdTrans(const double *inarray)
	{
	    LocalRegions::TriExpVectorIter  Tdef;
	    LocalRegions::QuadExpVectorIter Qdef;
	    int    cnt = 0;
	    
	    for(Tdef = m_tri.begin(); Tdef != m_tri.end(); ++Tdef)
	    {
		(*Tdef)->FwdTrans(inarray+cnt);
		cnt  += (*Tdef)->GetPointsTot();
	    }
	    
	    for(Qdef = m_quad.begin(); Qdef != m_quad.end(); ++Qdef)
	    {
		(*Qdef)->FwdTrans(inarray+cnt);
		cnt  += (*Qdef)->GetPointsTot();
	    }
	    
	    m_transState = eLocal;
	}
      
	void ExpList2D::BwdTrans(double *outarray)
	{
	    LocalRegions::TriExpVectorIter  Tdef;
	    LocalRegions::QuadExpVectorIter Qdef;
	    int    cnt = 0;
	  
	    for(Tdef = m_tri.begin(); Tdef != m_tri.end(); ++Tdef)
	    {
		(*Tdef)->BwdTrans(outarray+cnt);
		cnt  += (*Tdef)->GetPointsTot();
	    }

	    for(Qdef = m_quad.begin(); Qdef != m_quad.end(); ++Qdef)
	    {
		(*Qdef)->BwdTrans(outarray+cnt);
		cnt  += (*Qdef)->GetPointsTot();
	    }
	    
	    m_physState = true;
	}
	
	void ExpList2D::GetCoords(double **coords)
	{
	    LocalRegions::TriExpVectorIter  Tdef;
	    LocalRegions::QuadExpVectorIter Qdef;
	    int    i, cnt = 0;
	    double *E_coords[3];
	    
	    for(Tdef = m_tri.begin(); Tdef != m_tri.end(); ++Tdef)
	    {
		for(i = 0 ; i < (*Tdef)->GetCoordim(); ++i)
		{
		    E_coords[i] = coords[i]+cnt;
		}
	      
		(*Tdef)->GetCoords(E_coords);
		cnt  += (*Tdef)->GetPointsTot();
	    }

	    for(Qdef = m_quad.begin(); Qdef != m_quad.end(); ++Qdef)
	    {
		for(i = 0 ; i < (*Qdef)->GetCoordim(); ++i)
		{
		    E_coords[i] = coords[i]+cnt;
		}
		
		(*Qdef)->GetCoords(E_coords);
		cnt  += (*Qdef)->GetPointsTot();
	    }
	}
      
	void ExpList2D::WriteToFile(std::ofstream &out)
	{
	    LocalRegions::TriExpVectorIter Tdef; 
	    LocalRegions::QuadExpVectorIter Qdef; 
	    
	    if(m_physState == false)
	    {
		BwdTrans(m_phys);
	    }
	    
	    if(m_tri.size()){
		(*m_tri.begin())->WriteToFile(out,1);
		
		for(Tdef = ++m_tri.begin(); Tdef != m_tri.end(); ++Tdef)
		{
		    (*Tdef)->WriteToFile(out,0);
		}
	    }

	    if(m_quad.size()){
		(*m_quad.begin())->WriteToFile(out,1);
		
		for(Qdef = ++m_quad.begin(); Qdef != m_quad.end(); ++Qdef)
		{
		    (*Qdef)->WriteToFile(out,0);
		}
	    }
	}
      
	double  ExpList2D::Linf(const double *sol)
	{
	    LocalRegions::TriExpVectorIter  Tdef;
	    LocalRegions::QuadExpVectorIter Qdef;
	    double err = 0.0;
	    int    cnt = 0;
	    
	    if(m_physState == false)
	    {
		BwdTrans(m_phys);
	    }
	    
	    for(Tdef = m_tri.begin(); Tdef != m_tri.end(); ++Tdef)
	    {
		err  = std::max(err,(*Tdef)->Linf(sol+cnt));
		cnt  += (*Tdef)->GetPointsTot();
	    }

	    for(Qdef = m_quad.begin(); Qdef != m_quad.end(); ++Qdef)
	    {
		err  = std::max(err,(*Qdef)->Linf(sol+cnt));
		cnt  += (*Qdef)->GetPointsTot();
	    }
	    
	    return err;
	}

      
	double  ExpList2D::L2(const double *sol)
	{
	    LocalRegions::TriExpVectorIter Tdef;
	    LocalRegions::QuadExpVectorIter Qdef;
	    double err = 0.0,errl2;
	    int    cnt = 0;
	  
	    if(m_physState == false)
	    {
		BwdTrans(m_phys);
	    }
	    
	    for(Tdef = m_tri.begin(); Tdef != m_tri.end(); ++Tdef)
	    {
		errl2 = (*Tdef)->L2(sol+cnt);
		err += errl2*errl2;
		cnt  += (*Tdef)->GetPointsTot();
	    }

	    for(Qdef = m_quad.begin(); Qdef != m_quad.end(); ++Qdef)
	    {
		errl2 = (*Qdef)->L2(sol+cnt);
		err += errl2*errl2;
		cnt  += (*Qdef)->GetPointsTot();
	    }
	    
	    return sqrt(err);
	}
      
    } //end of namespace
} //end of namespace

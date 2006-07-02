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
	    
	    // determine size of local expansion and quadrature space
	    // and declare memory

	    if(TriGeoms.size())
	    {		
		tri_ncoeffs_elmt = (TriBa.GetBasisOrder()*(TriBa.GetBasisOrder()+1))/2 + TriBa.GetBasisOrder()*(TriBb.GetBasisOrder()-TriBa.GetBasisOrder());
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
		
	    
	    // make sure Geofacs are defined in MeshGraph2D
	    if(graph2D.GetGeofac_defined() != true)
	    {
		graph2D.GenXGeoFac();
	    }
	    

	    // declare triangles using first block of data 	    
	    if(TriGeoms.size())
	    {		
		LocalRegions::TriExpSharedPtr tri;
		SpatialDomains::TriGeomVectorIter def;
		StdRegions::StdExpansionVector explist;

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
		    explist.push_back(tri);
		    
		    cnt  += tri_ncoeffs_elmt;
		    cnt1 += tri_npoints_elmt;

		}
		m_exp_shapes.push_back(explist);
	    }
	    
	    // set up quads 
	    if(QuadGeoms.size())
	    {		
		LocalRegions::QuadExpSharedPtr quad;
		SpatialDomains::QuadGeomVectorIter def;
		StdRegions::StdExpansionVector explist;

		cnt = cnt1 = 0;
		for(def = QuadGeoms.begin(); def != QuadGeoms.end(); ++def)
		{
		    // removed copy construction of geom
		    // geom = new SpatialDomains::SegGeom (**def);
		    quad.reset(new LocalRegions::QuadExp(QuadBa,QuadBb,
							m_coeffs+cnt,m_phys+cnt1, 
							*def));
		    quad->SetGeoFac(quad->GenGeoFac());
		    explist.push_back(quad);
		    
		    cnt  += quad_ncoeffs_elmt;
		    cnt1 += quad_npoints_elmt;
		}
		m_exp_shapes.push_back(explist);
	    }
	}
      
    } //end of namespace
} //end of namespace

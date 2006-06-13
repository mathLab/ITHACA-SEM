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
	    
	    StdRegions::StdExpansionVector seglist;
	    SpatialDomains::SegGeomVectorIter def;
	    int cnt,cnt1;
	    
	    cnt = cnt1 = 0;
	    for(def = SegGeoms.begin(); def != SegGeoms.end(); ++def)
	    {
		seg.reset( new LocalRegions::SegExp(Ba, m_coeffs+cnt,
						    m_phys+cnt1, *def));
		seg->SetGeoFac(seg->GenGeoFac());
		seglist.push_back(seg);
		
		cnt  += Ba.GetBasisOrder();
		cnt1 += Ba.GetPointsOrder();
	    }
	    m_exp_shapes.push_back(seglist);
	}
	
	
    } //end of namespace
} //end of namespace

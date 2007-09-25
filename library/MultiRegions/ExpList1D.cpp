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
        
        ExpList1D::ExpList1D(const ExpList1D &In):
            ExpList(In)
        {
            m_coeffs = Array<OneD, NekDouble>(m_ncoeffs);
            m_phys   = Array<OneD, NekDouble>(m_npoints);
        }

        ExpList1D::ExpList1D(const LibUtilities::BasisKey &Ba,
                             const SpatialDomains::MeshGraph1D &graph1D)
        {
            int i,j;
            int nel;
            LocalRegions::SegExpSharedPtr seg;
            SpatialDomains::Composite comp;
            
            SpatialDomains::CompositeVector domain = (graph1D.GetDomain());
            
            m_ncoeffs = 0;
            m_npoints = 0;
            
            m_transState = eNotSet; 
            m_physState  = false;
            
            for(i = 0; i < domain.size(); ++i)
            {
                comp = domain[i];
                
                for(j = 0; j < comp->size(); ++j)
                {
                    SpatialDomains::SegGeomSharedPtr SegmentGeom;
                    
                    if(SegmentGeom = boost::dynamic_pointer_cast<SpatialDomains::SegGeom>((*comp)[j]))
                    {
                        seg = MemoryManager<LocalRegions::SegExp>::AllocateSharedPtr(Ba,SegmentGeom);
                        (*m_exp).push_back(seg);
                    }
                    else
                    {
                        ASSERTL0(false,"dynamic cast to a SegGeom failed");
                    }  
                }
                
                m_ncoeffs += (comp->size())*Ba.GetNumModes();
                m_npoints += (comp->size())*Ba.GetNumPoints();
            } 
            
            m_coeffs = Array<OneD, NekDouble>(m_ncoeffs);
            m_phys   = Array<OneD, NekDouble>(m_npoints);
            
        }

        ExpList1D::ExpList1D(SpatialDomains::MeshGraph1D &graph1D)
        {
            int i,j;
            int nel;
            LocalRegions::SegExpSharedPtr seg;
            SpatialDomains::Composite comp;
            
            SpatialDomains::CompositeVector domain = (graph1D.GetDomain());
            
            m_ncoeffs = 0;
            m_npoints = 0;
            
            m_transState = eNotSet; 
            m_physState  = false;
            
            for(i = 0; i < domain.size(); ++i)
            {
                comp = domain[i];
                LibUtilities::BasisKey bkey = graph1D.GetBasisKey(comp,0);
                
                for(j = 0; j < comp->size(); ++j)
                {
                    SpatialDomains::SegGeomSharedPtr SegmentGeom;
                    
                    if(SegmentGeom = boost::dynamic_pointer_cast<SpatialDomains::SegGeom>((*comp)[j]))
                    {
                        seg = MemoryManager<LocalRegions::SegExp>::AllocateSharedPtr(bkey, SegmentGeom);
                        (*m_exp).push_back(seg);
                    }
                    else
                    {
                        ASSERTL0(false,"dynamic cast to a SegGeom failed");
                    }  
                }
                
                m_ncoeffs += (comp->size())*bkey.GetNumModes();
                m_npoints += (comp->size())*bkey.GetNumPoints();
            } 
            
            m_coeffs = Array<OneD, NekDouble>(m_ncoeffs);
            m_phys   = Array<OneD, NekDouble>(m_npoints);
            
        }
    } //end of namespace
} //end of namespace

/**
* $Log: ExpList1D.cpp,v $
* Revision 1.19  2007/07/22 23:04:20  bnelson
* Backed out Nektar::ptr.
*
* Revision 1.18  2007/07/20 02:04:12  bnelson
* Replaced boost::shared_ptr with Nektar::ptr
*
* Revision 1.17  2007/07/13 16:48:47  pvos
* Another HelmHoltz update (homogeneous dir BC multi-elemental solver does work)
*
* Revision 1.16  2007/07/10 08:54:29  pvos
* Updated ContField1D constructor
*
* Revision 1.15  2007/07/06 18:39:34  pvos
* ContField1D constructor updates
*
* Revision 1.14  2007/06/05 16:36:55  pvos
* Updated Explist2D ContExpList2D and corresponding demo-codes
*
**/

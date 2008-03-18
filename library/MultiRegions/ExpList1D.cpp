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
    
        ExpList1D::ExpList1D():
            ExpList()
        {
        }
        
        ExpList1D::~ExpList1D()
        {
        }
        
        ExpList1D::ExpList1D(const ExpList1D &In):
            ExpList(In)
        {
        }

        ExpList1D::ExpList1D(const LibUtilities::BasisKey &Ba,
                             const SpatialDomains::MeshGraph1D &graph1D):
            ExpList()
        {
            int i,j;
            LocalRegions::SegExpSharedPtr seg;
            SpatialDomains::SegGeomSharedPtr SegmentGeom;

            const SpatialDomains::ExpansionVector &expansions = graph1D.GetExpansions();
            
            for(i = 0; i < expansions.size(); ++i)
            {                
                if(SegmentGeom = boost::dynamic_pointer_cast<SpatialDomains::SegGeom>(expansions[i]->m_GeomShPtr))
                {
                    seg = MemoryManager<LocalRegions::SegExp>::AllocateSharedPtr(Ba,SegmentGeom);
                    (*m_exp).push_back(seg);
                }
                else
                {
                    ASSERTL0(false,"dynamic cast to a SegGeom failed");
                }  
            
                m_ncoeffs += Ba.GetNumModes();
                m_npoints += Ba.GetNumPoints();
            } 
            
            m_coeffs = Array<OneD, NekDouble>(m_ncoeffs);
            m_phys   = Array<OneD, NekDouble>(m_npoints);            
        }

        ExpList1D::ExpList1D(SpatialDomains::MeshGraph1D &graph1D):
            ExpList()
        {
            int i;
            LocalRegions::SegExpSharedPtr seg;
            SpatialDomains::SegGeomSharedPtr SegmentGeom;

            const SpatialDomains::ExpansionVector &expansions = graph1D.GetExpansions();
            
            for(i = 0; i < expansions.size(); ++i)
            {
                LibUtilities::BasisKey bkey = graph1D.GetBasisKey(expansions[i],0);

                if(SegmentGeom = boost::dynamic_pointer_cast<SpatialDomains::SegGeom>(expansions[i]->m_GeomShPtr))
                {
                    seg = MemoryManager<LocalRegions::SegExp>::AllocateSharedPtr(bkey, SegmentGeom);
                    (*m_exp).push_back(seg);
                }
                else
                {
                    ASSERTL0(false,"dynamic cast to a SegGeom failed");
                }  

                m_ncoeffs += bkey.GetNumModes();
                m_npoints += bkey.GetNumPoints();
            }
      
            m_coeffs = Array<OneD, NekDouble>(m_ncoeffs);
            m_phys   = Array<OneD, NekDouble>(m_npoints);
            
        }

        ExpList1D::ExpList1D(const SpatialDomains::CompositeVector &domain, SpatialDomains::MeshGraph2D &graph2D):
            ExpList()
        {
            int i,j;
            SpatialDomains::Composite comp;
            SpatialDomains::SegGeomSharedPtr SegmentGeom;
            LocalRegions::SegExpSharedPtr seg;
            // SpatialDomains::ElementEdgeVectorSharedPtr edgeElement;
            //SpatialDomains::ExpansionShPtr exp;
            
            for(i = 0; i < domain.size(); ++i)
            {
                comp = domain[i];
                
                for(j = 0; j < comp->size(); ++j)
                {                    
                    if(SegmentGeom = boost::dynamic_pointer_cast<SpatialDomains::SegGeom>((*comp)[j]))
                    {
                        // edgeElement = graph2D.GetElementsFromEdge(SegmentGeom);
//                         ASSERTL1(edgeElement->size()==1,"This boundary edge is bounded by more than one element");
//                         exp = graph2D.GetExpansion((*edgeElement)[0]->m_Element);
//                         // (are we sure that flag should be set to zero (for triangles or anisotropic expansions for example?))
//                         LibUtilities::BasisKey bkey = graph2D.GetBasisKey(exp,0);

                        LibUtilities::BasisKey bkey = graph2D.GetEdgeBasisKey(SegmentGeom);

                        // Do a fix for nodal triangular expansions
                      //   if((exp->m_ExpansionType == SpatialDomains::eNodal)&&(bkey.GetBasisType()!=LibUtilities::eGLL_Lagrange))
//                         {
//                             LibUtilities::BasisKey bkey2(LibUtilities::eGLL_Lagrange,bkey.GetNumModes(),bkey.GetPointsKey());
//                             seg = MemoryManager<LocalRegions::SegExp>::AllocateSharedPtr(bkey2, SegmentGeom);
//                         }
//                         else
//                         {
//                             seg = MemoryManager<LocalRegions::SegExp>::AllocateSharedPtr(bkey, SegmentGeom);
//                         }

                        seg = MemoryManager<LocalRegions::SegExp>::AllocateSharedPtr(bkey, SegmentGeom);

                        (*m_exp).push_back(seg);  

                        m_ncoeffs += bkey.GetNumModes();
                        m_npoints += bkey.GetNumPoints();
                    }
                    else
                    {
                        ASSERTL0(false,"dynamic cast to a EdgeComp failed");
                    }  
                }                
            } 
            
            m_coeffs = Array<OneD, NekDouble>(m_ncoeffs);
            m_phys   = Array<OneD, NekDouble>(m_npoints);            
        }

    } //end of namespace
} //end of namespace

/**
* $Log: ExpList1D.cpp,v $
* Revision 1.23  2008/03/12 15:25:45  pvos
* Clean up of the code
*
* Revision 1.22  2007/12/06 22:52:30  pvos
* 2D Helmholtz solver updates
*
* Revision 1.21  2007/11/07 20:29:53  jfrazier
* Modified to use new expansion list contained in meshgraph.
*
* Revision 1.20  2007/09/25 14:25:29  pvos
* Update for helmholtz1D with different expansion orders
*
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

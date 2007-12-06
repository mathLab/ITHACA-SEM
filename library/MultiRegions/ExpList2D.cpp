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
        
        ExpList2D::ExpList2D(const ExpList2D &In):
            ExpList(In)
        {
            m_coeffs = Array<OneD, NekDouble>(m_ncoeffs);
            m_phys   = Array<OneD, NekDouble>(m_npoints);
        }
        
        ExpList2D::ExpList2D(const LibUtilities::BasisKey &TriBa, 
                             const LibUtilities::BasisKey &TriBb, 
                             const LibUtilities::BasisKey &QuadBa, 
                             const LibUtilities::BasisKey &QuadBb, 
                             const SpatialDomains::MeshGraph2D &graph2D,
                             const LibUtilities::PointsType TriNb)
        {
            int i,j;
            int nel;
            LocalRegions::TriExpSharedPtr tri;
            LocalRegions::NodalTriExpSharedPtr Ntri;
            LocalRegions::QuadExpSharedPtr quad;
            SpatialDomains::Composite comp;

            const SpatialDomains::ExpansionVector &expansions = graph2D.GetExpansions();
            
            m_ncoeffs = 0;
            m_npoints = 0;
            
            m_transState = eNotSet; 
            m_physState  = false;
            
            for(i = 0; i < expansions.size(); ++i)
            {
                SpatialDomains::TriGeomSharedPtr TriangleGeom;
                SpatialDomains::QuadGeomSharedPtr QuadrilateralGeom;
                
                if(TriangleGeom = boost::dynamic_pointer_cast<SpatialDomains::TriGeom>(expansions[i]->m_GeomShPtr))
                {
                    if(TriNb < LibUtilities::SIZE_PointsType)
                    {
                        Ntri = MemoryManager<LocalRegions::NodalTriExp>::AllocateSharedPtr(TriBa,TriBb,TriNb,TriangleGeom);
                        (*m_exp).push_back(Ntri);
                    }
                    else
                    {
                        tri = MemoryManager<LocalRegions::TriExp>::AllocateSharedPtr(TriBa,TriBb,TriangleGeom);
                        (*m_exp).push_back(tri);
                    }
                    
                    m_ncoeffs += (TriBa.GetNumModes()*(TriBa.GetNumModes()+1))/2 
                        + TriBa.GetNumModes()*(TriBb.GetNumModes()-TriBa.GetNumModes());
                    m_npoints += TriBa.GetNumPoints()*TriBb.GetNumPoints();
                }
                else if(QuadrilateralGeom = boost::dynamic_pointer_cast<SpatialDomains::QuadGeom>(expansions[i]->m_GeomShPtr))
                {
                    quad = MemoryManager<LocalRegions::QuadExp>::AllocateSharedPtr(QuadBa,QuadBb,QuadrilateralGeom);
                    (*m_exp).push_back(quad);
                    m_ncoeffs += QuadBa.GetNumModes()*QuadBb.GetNumModes();
                    m_npoints += QuadBa.GetNumPoints()*QuadBb.GetNumPoints();
                }
                else
                {
                    ASSERTL0(false,"dynamic cast to a proper Geometry2D failed");
                }  
                
            }
            
            m_coeffs = Array<OneD, NekDouble>(m_ncoeffs);
            m_phys   = Array<OneD, NekDouble>(m_npoints);
        }



        ExpList2D::ExpList2D(SpatialDomains::MeshGraph2D &graph2D)
        {
            int i,j;
            int nel;
            LocalRegions::TriExpSharedPtr tri;
            LocalRegions::NodalTriExpSharedPtr Ntri;
            LocalRegions::QuadExpSharedPtr quad;
            SpatialDomains::Composite comp;

            const SpatialDomains::ExpansionVector &expansions = graph2D.GetExpansions();            
            
            m_ncoeffs = 0;
            m_npoints = 0;
            
            m_transState = eNotSet; 
            m_physState  = false;
            
            for(i = 0; i < expansions.size(); ++i)
            {
                SpatialDomains::TriGeomSharedPtr TriangleGeom;
                SpatialDomains::QuadGeomSharedPtr QuadrilateralGeom;
                
                if(TriangleGeom = boost::dynamic_pointer_cast<SpatialDomains::TriGeom>(expansions[i]->m_GeomShPtr))
                {
                    LibUtilities::BasisKey TriBa = graph2D.GetBasisKey(expansions[i],0);
                    LibUtilities::BasisKey TriBb = graph2D.GetBasisKey(expansions[i],1);
                    
                    tri = MemoryManager<LocalRegions::TriExp>::AllocateSharedPtr(TriBa,TriBb,TriangleGeom);
                    (*m_exp).push_back(tri);
                    
                    m_ncoeffs += (TriBa.GetNumModes()*(TriBa.GetNumModes()+1))/2 
                        + TriBa.GetNumModes()*(TriBb.GetNumModes()-TriBa.GetNumModes());
                    m_npoints += TriBa.GetNumPoints()*TriBb.GetNumPoints();
                }
                else if(QuadrilateralGeom = boost::dynamic_pointer_cast<SpatialDomains::QuadGeom>(expansions[i]->m_GeomShPtr))
                {
                    LibUtilities::BasisKey QuadBa = graph2D.GetBasisKey(expansions[i],0);
                    LibUtilities::BasisKey QuadBb = graph2D.GetBasisKey(expansions[i],0);
                    
                    quad = MemoryManager<LocalRegions::QuadExp>::AllocateSharedPtr(QuadBa,QuadBb,QuadrilateralGeom);
                    (*m_exp).push_back(quad);
                    
                    m_ncoeffs += QuadBa.GetNumModes()*QuadBb.GetNumModes();
                    m_npoints += QuadBa.GetNumPoints()*QuadBb.GetNumPoints();
                }
                else
                {
                    ASSERTL0(false,"dynamic cast to a proper Geometry2D failed");
                }  
                
            }            
            m_coeffs = Array<OneD, NekDouble>(m_ncoeffs);
            m_phys   = Array<OneD, NekDouble>(m_npoints);
        }

    } //end of namespace
} //end of namespace

/**
* $Log: ExpList2D.cpp,v $
* Revision 1.14  2007/07/20 02:04:12  bnelson
* Replaced boost::shared_ptr with Nektar::ptr
*
* Revision 1.13  2007/06/07 15:54:19  pvos
* Modificications to make Demos/MultiRegions/ProjectCont2D work correctly.
* Also made corrections to various ASSERTL2 calls
*
* Revision 1.12  2007/06/05 16:36:55  pvos
* Updated Explist2D ContExpList2D and corresponding demo-codes
*
**/

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
        SpatialDomains::TriGeomVector TriGeoms   = graph2D.GetTrigeoms();
        SpatialDomains::QuadGeomVector QuadGeoms = graph2D.GetQuadgeoms();

        int tri_ncoeffs_elmt = 0;
            int quad_ncoeffs_elmt = 0;
        int tri_npoints_elmt = 0;
            int quad_npoints_elmt = 0;
        
        // determine size of local expansion and quadrature space
        // and declare memory

            m_ncoeffs = 0;
            m_npoints = 0;
            
        if(TriGeoms.size())
        {        
        tri_ncoeffs_elmt = tri_ncoeffs_elmt = (TriBa.GetNumModes()*(TriBa.GetNumModes()+1))/2 
                    + TriBa.GetNumModes()*(TriBb.GetNumModes()-TriBa.GetNumModes());
        tri_npoints_elmt = (TriBa.GetNumPoints()*
                    TriBb.GetNumPoints());
        
        m_ncoeffs += TriGeoms.size()*tri_ncoeffs_elmt;
        m_npoints += TriGeoms.size()*tri_npoints_elmt;
        }        
        
        if(QuadGeoms.size())
        {        
        quad_ncoeffs_elmt = (QuadBa.GetNumModes() 
                                     *QuadBb.GetNumModes());
        quad_npoints_elmt = (QuadBa.GetNumPoints()
                                     *QuadBb.GetNumPoints());
        
        m_ncoeffs += QuadGeoms.size()*quad_ncoeffs_elmt;
        m_npoints += QuadGeoms.size()*quad_npoints_elmt;
        }
        
        m_transState = eNotSet; 
        m_physState  = false;

            m_coeffs = Array<OneD, NekDouble>(m_ncoeffs);
            m_phys   = Array<OneD, NekDouble>(m_npoints);
        
          // declare triangles using first block of data         
        if(TriGeoms.size())
        {        
        LocalRegions::TriExpSharedPtr tri;
        LocalRegions::NodalTriExpSharedPtr Ntri;
        SpatialDomains::TriGeomVectorIter def;
        
        for(def = TriGeoms.begin(); def != TriGeoms.end(); ++def)
        {
            
            if(TriNb < LibUtilities::SIZE_PointsType)
            {

                        Ntri = MemoryManager<LocalRegions::NodalTriExp>::AllocateSharedPtr(TriBa,TriBb,TriNb,*def);
                        (*m_exp).push_back(Ntri);
            }
            else
            {
                        tri = MemoryManager<LocalRegions::TriExp>::AllocateSharedPtr(TriBa,TriBb,*def);
                        (*m_exp).push_back(tri);
            }
        }
        }
        
        // set up quads 
        if(QuadGeoms.size())
        {        
        LocalRegions::QuadExpSharedPtr quad;
        SpatialDomains::QuadGeomVectorIter def;
                
        for(def = QuadGeoms.begin(); def != QuadGeoms.end(); ++def)
        {
                    quad = MemoryManager<LocalRegions::QuadExp>::AllocateSharedPtr(QuadBa,QuadBb,*def);
                    (*m_exp).push_back(quad);
        }
        }
    }
        

//     ExpList2D::ExpList2D(const StdRegions::BasisKey  &TriBa, 
//                  const StdRegions::BasisKey  &TriBb, 
//                  const StdRegions::BasisKey  &QuadBa, 
//                  const StdRegions::BasisKey  &QuadBb, 
//                  SpatialDomains::MeshGraph2D &graph2D,
//                  SpatialDomains::Domain      &domain2D,
//                  const StdRegions::NodalBasisType TriNb )
//     {
//         SpatialDomains::TriGeomVector TriGeoms   = graph2D.GetTrigeoms();
//         SpatialDomains::QuadGeomVector QuadGeoms = graph2D.GetQuadgeoms();
//         int cnt,cnt1;
//         int tri_ncoeffs_elmt = 0, quad_ncoeffs_elmt = 0;
//         int tri_npoints_elmt = 0, quad_npoints_elmt = 0;
//         CompositeVector compVector = domain2D.GetDomain();
//         BoundaryVector boundVector = domain2D.GetBoundaries();

        
//         // determine size of local expansion and quadrature space
//         // and declare memory

//         m_npoints = m_ncoeffs = 0;


//         for(
//         if(TriGeoms.size())
//         {        
//         tri_ncoeffs_elmt = (TriBa.GetBasisOrder()*(TriBa.GetBasisOrder()+1))/2 + TriBa.GetBasisOrder()*(TriBb.GetBasisOrder()-TriBa.GetBasisOrder());
//         tri_npoints_elmt = (TriBa.GetPointsOrder()*
//                     TriBb.GetPointsOrder());
        
//         m_ncoeffs += TriGeoms.size()*tri_ncoeffs_elmt;
//         m_npoints += TriGeoms.size()*tri_npoints_elmt;
//         }        
        
//         if(QuadGeoms.size())
//         {        
//         quad_ncoeffs_elmt = QuadBa.GetBasisOrder() 
//             *QuadBb.GetBasisOrder();
//         quad_npoints_elmt = QuadBa.GetPointsOrder()
//             *QuadBb.GetPointsOrder();
        
//         m_ncoeffs += QuadGeoms.size()*quad_ncoeffs_elmt;
//         m_npoints += QuadGeoms.size()*quad_npoints_elmt;
//         }
        
//         m_coeffs = new double [m_ncoeffs];
//         m_transState = eNotSet; 
        
//         m_phys   = new double [m_npoints];
//         m_physState  = false;
        
        
//         // make sure Geofacs are defined in MeshGraph2D
//         if(graph2D.GetGeofac_defined() != true)
//         {
//         graph2D.GenXGeoFac();
//         }
        

//         // declare triangles using first block of data         
//         cnt = cnt1 = 0; // use these counts for data offsets
//         if(TriGeoms.size())
//         {        
//         LocalRegions::TriExpSharedPtr tri;
//         LocalRegions::NodalTriExpSharedPtr Ntri;
//         SpatialDomains::TriGeomVectorIter def;
//         StdRegions::StdExpansionVector explist;

//         // make sure Geofacs are defined in MeshGraph1D
//         if(graph2D.GetGeofac_defined() != true)
//         {
//             graph2D.GenXGeoFac();
//         }
        
//         for(def = TriGeoms.begin(); def != TriGeoms.end(); ++def)
//         {
//             // removed copy construction of geom
//             // geom = new SpatialDomains::SegGeom (**def);
            
//             if(TriNb < StdRegions::SIZE_NodalBasisType)
//             {
//             Ntri.reset(new LocalRegions::NodalTriExp(TriBa,TriBb,
//                          TriNb, m_coeffs+cnt,
//                          m_phys+cnt1, *def));
//             Ntri->SetGeoFac(Ntri->GenGeoFac());
//             explist.push_back(Ntri);
//             }
//             else
//             {
//             tri.reset(new LocalRegions::TriExp(TriBa,TriBb,
//                           m_coeffs+cnt,m_phys+cnt1, 
//                                *def));
//             tri->SetGeoFac(tri->GenGeoFac());
//             explist.push_back(tri);
//             }
//             cnt  += tri_ncoeffs_elmt;
//             cnt1 += tri_npoints_elmt;

//         }
//         m_exp_shapes.push_back(explist);
//         }
        
//         // set up quads 
//         if(QuadGeoms.size())
//         {        
//         LocalRegions::QuadExpSharedPtr quad;
//         SpatialDomains::QuadGeomVectorIter def;
//         StdRegions::StdExpansionVector explist;

//         for(def = QuadGeoms.begin(); def != QuadGeoms.end(); ++def)
//         {
//             // removed copy construction of geom
//             // geom = new SpatialDomains::SegGeom (**def);
//             quad.reset(new LocalRegions::QuadExp(QuadBa, QuadBb,
//                     m_coeffs+cnt, m_phys+cnt1, *def));
//             quad->SetGeoFac(quad->GenGeoFac());
//             explist.push_back(quad);
            
//             cnt  += quad_ncoeffs_elmt;
//             cnt1 += quad_npoints_elmt;
//         }
//         m_exp_shapes.push_back(explist);
//         }
//     }

      
    } //end of namespace
} //end of namespace

/**
* $Log: ExpList2D.cpp,v $
* Revision 1.13  2007/06/07 15:54:19  pvos
* Modificications to make Demos/MultiRegions/ProjectCont2D work correctly.
* Also made corrections to various ASSERTL2 calls
*
* Revision 1.12  2007/06/05 16:36:55  pvos
* Updated Explist2D ContExpList2D and corresponding demo-codes
*
**/

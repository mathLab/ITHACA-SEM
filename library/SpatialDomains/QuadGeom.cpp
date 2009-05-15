////////////////////////////////////////////////////////////////////////////////
//
//  File: QuadGeom.cpp
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific 
//  Computing and Imaging Institute, University of Utah (USA).
//
//  License for the specific language governing rights and limitations under
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description:  
//
//
////////////////////////////////////////////////////////////////////////////////
#include "pchSpatialDomains.h"

#include <SpatialDomains/QuadGeom.h>

namespace Nektar
{
    namespace SpatialDomains
    {
        /// Default constructor
        
        QuadGeom::QuadGeom()
        {
            m_GeomShapeType = eQuadrilateral;
        }

        QuadGeom::QuadGeom(int id, const int coordim):
                          Geometry2D(coordim), m_fid(id)
        {

            const LibUtilities::BasisKey B(LibUtilities::eModified_A, 2,
            LibUtilities::PointsKey(3,LibUtilities::eGaussLobattoLegendre));

            m_xmap = Array<OneD, StdRegions::StdExpansion2DSharedPtr>(m_coordim);
            
            m_fid = id;

            for(int i = 0; i < m_coordim; ++i)
            {
                m_xmap[i] = MemoryManager<StdRegions::StdQuadExp>::AllocateSharedPtr(B,B);  
            }

        }

        QuadGeom::QuadGeom(const int id, 
                           const VertexComponentSharedPtr verts[], 
                           const SegGeomSharedPtr edges[], 
                           const StdRegions::EdgeOrientation eorient[]):
            Geometry2D(verts[0]->GetCoordim()),
            m_fid(id)
        {
            m_GeomShapeType = eQuadrilateral;

            /// Copy the vert shared pointers.
            m_verts.insert(m_verts.begin(), verts, verts+QuadGeom::kNverts);

            /// Copy the edge shared pointers.
            m_edges.insert(m_edges.begin(), edges, edges+QuadGeom::kNedges);


            for (int j=0; j<kNedges; ++j)
            {
                m_eorient[j] = eorient[j];
            }

            m_coordim = verts[0]->GetCoordim();
            ASSERTL0(m_coordim > 1,
                     "Cannot call function with dim == 1");

            int order0  = max(edges[0]->GetBasis(0,0)->GetNumModes(),
                             edges[2]->GetBasis(0,0)->GetNumModes());
            int points0 = max(edges[0]->GetBasis(0,0)->GetNumPoints(),
                             edges[2]->GetBasis(0,0)->GetNumPoints());
            int order1 = max(edges[1]->GetBasis(0,0)->GetNumModes(),
                             edges[3]->GetBasis(0,0)->GetNumModes());
            int points1 = max(edges[1]->GetBasis(0,0)->GetNumPoints(),
                             edges[3]->GetBasis(0,0)->GetNumPoints());

            const LibUtilities::BasisKey B0(LibUtilities::eModified_A, order0,
                  LibUtilities::PointsKey(points0,LibUtilities::eGaussLobattoLegendre));
            const LibUtilities::BasisKey B1(LibUtilities::eModified_A, order1,
                  LibUtilities::PointsKey(points1,LibUtilities::eGaussLobattoLegendre));

            m_xmap = Array<OneD, StdRegions::StdExpansion2DSharedPtr>(m_coordim);

            for(int i = 0; i < m_coordim; ++i)
            {
                m_xmap[i] = MemoryManager<StdRegions::StdQuadExp>::AllocateSharedPtr(B0,B1);  
            }
        }
    
        QuadGeom::QuadGeom(const int id, 
                           const SegGeomSharedPtr edges[], 
                           const StdRegions::EdgeOrientation eorient[]):
            Geometry2D(edges[0]->GetVertex(0)->GetCoordim()),
            m_fid(id)
        {
            int j;
            
            m_GeomShapeType = eQuadrilateral;
            
            /// Copy the edge shared pointers.
            m_edges.insert(m_edges.begin(), edges, edges+QuadGeom::kNedges);
            
            for(j=0; j <kNedges; ++j)
            {
                if(eorient[j] == StdRegions::eForwards)
                {
                    m_verts.push_back(edges[j]->GetVertex(0));
                }
                else
                {
                    m_verts.push_back(edges[j]->GetVertex(1));
                }
            }
            
            for (j=0; j<kNedges; ++j)
            {
                m_eorient[j] = eorient[j];
            }

            m_coordim = edges[0]->GetVertex(0)->GetCoordim();
            ASSERTL0(m_coordim > 1,
                "Cannot call function with dim == 1");

            int order0  = max(edges[0]->GetBasis(0,0)->GetNumModes(),
                             edges[2]->GetBasis(0,0)->GetNumModes());
            int points0 = max(edges[0]->GetBasis(0,0)->GetNumPoints(),
                             edges[2]->GetBasis(0,0)->GetNumPoints());
            int order1 = max(edges[1]->GetBasis(0,0)->GetNumModes(),
                             edges[3]->GetBasis(0,0)->GetNumModes());
            int points1 = max(edges[1]->GetBasis(0,0)->GetNumPoints(),
                             edges[3]->GetBasis(0,0)->GetNumPoints());

            const LibUtilities::BasisKey B0(LibUtilities::eModified_A, order0,
                  LibUtilities::PointsKey(points0,LibUtilities::eGaussLobattoLegendre));
            const LibUtilities::BasisKey B1(LibUtilities::eModified_A, order1,
                  LibUtilities::PointsKey(points1,LibUtilities::eGaussLobattoLegendre));

            m_xmap = Array<OneD, StdRegions::StdExpansion2DSharedPtr>(m_coordim);

            for(int i = 0; i < m_coordim; ++i)
            {
                m_xmap[i] = MemoryManager<StdRegions::StdQuadExp>::AllocateSharedPtr(B0,B1);  
            }
        }

        QuadGeom::QuadGeom(const QuadGeom &in)
        {
            // From Geometry
            m_GeomShapeType = in.m_GeomShapeType;
            
            // From QuadFaceComponent
            m_fid = in.m_fid;
			m_ownverts = in.m_ownverts;
			std::list<CompToElmt>::const_iterator def;
            for(def = in.m_elmtmap.begin(); def != in.m_elmtmap.end(); def++)
            {
                m_elmtmap.push_back(*def);
            }

			// From QuadGeom
			m_verts = in.m_verts;
			m_edges = in.m_edges;
            for (int i = 0; i < kNedges; i++)
            {
                m_eorient[i] = in.m_eorient[i];
            }
            m_owndata = in.m_owndata;
        }

        QuadGeom::~QuadGeom()
        {
        }

        void QuadGeom::AddElmtConnected(int gvo_id, int locid)
        {
            CompToElmt ee(gvo_id,locid);
            m_elmtmap.push_back(ee);
        }

        int QuadGeom::NumElmtConnected() const
        {
            return int(m_elmtmap.size());
        }

        bool QuadGeom::IsElmtConnected(int gvo_id, int locid) const
        {
            std::list<CompToElmt>::const_iterator def;
            CompToElmt ee(gvo_id,locid);

            def = find(m_elmtmap.begin(),m_elmtmap.end(),ee);

            // Found the element connectivity object in the list
            if(def != m_elmtmap.end())
            {
                return(true);
            }

            return(false);
        }

        NekDouble QuadGeom::GetCoord(const int i, 
                                  const Array<OneD, const NekDouble> &Lcoord)
        {
            ASSERTL1(m_state == ePtsFilled,
                "Goemetry is not in physical space");

            return m_xmap[i]->PhysEvaluate(Lcoord);
        }

        // Set up GeoFac for this geometry using Coord quadrature distribution
        void QuadGeom::GenGeomFactors(const Array<OneD, const LibUtilities::BasisSharedPtr> &tbasis)
        {
            int i;
            GeomType Gtype = eRegular;

            FillGeom();

            // We will first check whether we have a regular or deformed geometry.
            // We will define regular as those cases where the Jacobian and the metric
            // terms of the derivative are constants (i.e. not coordinate dependent)

            // Check to see if expansions are linear
            // If not linear => deformed geometry
            for(i = 0; i < m_coordim; ++i)
            {
                if((m_xmap[i]->GetBasisNumModes(0) != 2)||
                   (m_xmap[i]->GetBasisNumModes(1) != 2))
                {
                    Gtype = eDeformed;
                }
            }

            // For linear expansions, the mapping from standard to local 
            // element is given by the relation:
            // x_i = 0.25 * [ ( x_i^A + x_i^B + x_i^C + x_i^D)       + 
            //                (-x_i^A + x_i^B + x_i^C - x_i^D)*xi_1  + 
            //                (-x_i^A - x_i^B + x_i^C + x_i^D)*xi_2  + 
            //                ( x_i^A - x_i^B + x_i^C - x_i^D)*xi_1*xi_2 ] 
            //
            // The jacobian of the transformation and the metric terms dxi_i/dx_j,
            // involve only terms of the form dx_i/dxi_j (both for coordim == 2 or 3).
            // Inspecting the formula above, it can be appreciated that the derivatives
            // dx_i/dxi_j will be constant, if the coefficient of the non-linear term
            // is zero.
            //
            // That is why for regular geometry, we require
            //
            //     x_i^A - x_i^B + x_i^C - x_i^D = 0
            //
            // or equivalently
            //
            //     x_i^A - x_i^B = x_i^D - x_i^C
            //
            // This corresponds to quadrilaterals which are paralellograms.
            if(Gtype == eRegular)
            {
                for(i = 0; i < m_coordim; i++)
                {
                    if( fabs( (*m_verts[0])(i) - (*m_verts[1])(i) + 
                              (*m_verts[2])(i) - (*m_verts[3])(i) ) > NekConstants::kNekZeroTol )
                    {
                        Gtype = eDeformed;
                        break;
                    }
                }
            }

            m_geomfactors = MemoryManager<GeomFactors>::AllocateSharedPtr(Gtype, m_coordim, m_xmap, tbasis);
        }

        /** \brief put all quadrature information into edge structure 
        and backward transform 

        Note verts and edges are listed according to anticlockwise
        convention but points in _coeffs have to be in array format from
        left to right.

        */

        void QuadGeom::FillGeom()
        {            
            // check to see if geometry structure is already filled
            if(m_state != ePtsFilled)
            {
                int i,j,k;
                int nEdgeCoeffs;

                Array<OneD, unsigned int> mapArray;
                Array<OneD, int>          signArray;

                for(i = 0; i < kNedges; i++)
                {
                    m_edges[i]->FillGeom();
                    m_xmap[0]->GetEdgeToElementMap(i,m_eorient[i],
                                                   mapArray,signArray);

                    nEdgeCoeffs = (*m_edges[i])[0]->GetNcoeffs();
                    
                    for(j = 0 ; j < m_coordim; j++)
                    {
                        for(k = 0; k < nEdgeCoeffs; k++)
                        {
                            (m_xmap[j]->UpdateCoeffs())[mapArray[k]] 
                                = signArray[k]*((*m_edges[i])[j]->GetCoeffs())[k];
                        }
                    }
                }
                
                for(i = 0; i < m_coordim; ++i)
                {
                    m_xmap[i]->BwdTrans(m_xmap[i]->GetCoeffs(),
                                        m_xmap[i]->UpdatePhys());
                }
                
                m_state = ePtsFilled;
            }
        }
        
        void QuadGeom::GetLocCoords(const Array<OneD, const NekDouble> &coords,
                                    Array<OneD,NekDouble> &Lcoords)
        {
            int i;
            
            FillGeom();

            // calculate local coordinate for coord 
            if(GetGtype() == eRegular)
            { // can assume it is right angled rectangle
                NekDouble len0 = 0.0 ;
                NekDouble len1 = 0.0;
                NekDouble xi0 = 0.0;
                NekDouble xi1 = 0.0;
                Array<OneD, const NekDouble> pts;
                int nq0, nq1;

                // get points;
                //find end points
                for(i = 0; i < m_coordim; ++i)
                {
                    nq0 = m_xmap[i]->GetNumPoints(0);
                    nq1 = m_xmap[i]->GetNumPoints(1);

                    pts = m_xmap[i]->GetPhys();

                    // use projection to side 1 to determine xi_1 coordinate based on length
                    len0 += (pts[nq0-1]-pts[0])*(pts[nq0-1]-pts[0]);    
                    xi0  += (coords[i] -pts[0])*(pts[nq0-1]-pts[0]);

                    // use projection to side 4 to determine xi_2 coordinate based on length
                    len1 += (pts[nq0*(nq1-1)]-pts[0])*(pts[nq0*(nq1-1)]-pts[0]); 
                    xi1  += (coords[i] -pts[0])*(pts[nq0*(nq1-1)]-pts[0]);    
                }

                Lcoords[0] =  2*xi0/len0-1.0;
                Lcoords[1] =  2*xi1/len1-1.0;

            }
            else
            {
                NEKERROR(ErrorUtil::efatal,
                         "inverse mapping must be set up to use this call");
            }
        }

        //TODO: implement eight different case of face orientation
       StdRegions::FaceOrientation QuadGeom::GetFaceOrientation(const QuadGeom &face1,
                                                                const QuadGeom &face2)
       {
            StdRegions::FaceOrientation returnval = StdRegions::eDir1FwdDir1_Dir2FwdDir2;

             // TODO : implement 
//             eDir1FwdDir1_Dir2BwdDir2 
//             eDir1BwdDir1_Dir2FwdDir2 
//             eDir1BwdDir1_Dir2BwdDir2 
//             eDir1FwdDir2_Dir2FwdDir1 
//             eDir1FwdDir2_Dir2BwdDir1 
//             eDir1BwdDir2_Dir2FwdDir1 
//             eDir1BwdDir2_Dir2BwdDir1 
   
            return returnval;
        }

    }; //end of namespace
}; //end of namespace

//
// $Log: QuadGeom.cpp,v $
// Revision 1.24  2009/01/21 16:59:03  pvos
// Added additional geometric factors to improve efficiency
//
// Revision 1.23  2008/12/18 14:08:58  pvos
// NekConstants update
//
// Revision 1.22  2008/09/12 11:26:19  pvos
// Updates for mappings in 3D
//
// Revision 1.21  2008/09/09 14:26:22  sherwin
// Updates for deformed curved quads
//
// Revision 1.20  2008/08/14 22:11:03  sherwin
// Mods for HDG update
//
// Revision 1.19  2008/06/30 19:35:22  ehan
// Fixed infinity recursive-loop error.
//
// Revision 1.18  2008/06/14 01:23:07  ehan
// Implemented constructor and FillGeom().
//
// Revision 1.17  2008/06/11 21:34:42  delisi
// Removed TriFaceComponent, QuadFaceComponent, and EdgeComponent.
//
// Revision 1.16  2008/05/29 19:01:45  delisi
// Renamed eQuad to eQuadrilateral.
//
// Revision 1.15  2008/05/28 21:52:27  jfrazier
// Added GeomShapeType initialization for the different shapes.
//
// Revision 1.14  2008/05/07 16:05:37  pvos
// Mapping + Manager updates
//
// Revision 1.13  2008/04/06 06:00:38  bnelson
// Changed ConstArray to Array<const>
//
// Revision 1.12  2008/01/21 19:58:14  sherwin
// Updated so that QuadGeom and TriGeom have SegGeoms instead of EdgeComponents
//
// Revision 1.11  2007/07/20 02:15:09  bnelson
// Replaced boost::shared_ptr with Nektar::ptr
//
// Revision 1.10  2007/06/06 15:15:21  pvos
// Some minor updates for 2D routines
//
// Revision 1.9  2007/06/06 11:29:31  pvos
// Changed ErrorUtil::Error into NEKERROR (modifications in ErrorUtil.hpp caused compiler errors)
//
// Revision 1.8  2007/05/28 21:48:42  sherwin
// Update for 2D functionality
//
// Revision 1.7  2006/08/05 19:03:47  sherwin
// Update to make the multiregions 2D expansion in connected regions work
//
// Revision 1.6  2006/07/02 17:16:17  sherwin
//
// Modifications to make MultiRegions work for a connected domain in 2D (Tris)
//
// Revision 1.5  2006/06/02 18:48:40  sherwin
// Modifications to make ProjectLoc2D run bit there are bus errors for order > 3
//
// Revision 1.4  2006/06/01 14:15:30  sherwin
// Added typdef of boost wrappers and made GeoFac a boost shared pointer.
//
// Revision 1.3  2006/05/16 22:28:31  sherwin
// Updates to add in FaceComponent call to constructors
//
// Revision 1.2  2006/05/07 11:26:38  sherwin
// Modifications to get the demo LocalRegions::Project2D to compile
//
// Revision 1.1  2006/05/04 18:59:03  kirby
// *** empty log message ***
//
// Revision 1.23  2006/05/02 21:21:11  sherwin
// Corrected libraries to compile new version of spatialdomains and demo Graph1D
//
// Revision 1.22  2006/04/11 23:18:11  jfrazier
// Completed MeshGraph2D for tri's and quads.  Not thoroughly tested.
//
// Revision 1.21  2006/04/09 02:08:35  jfrazier
// Added precompiled header.
//
// Revision 1.20  2006/03/25 00:58:29  jfrazier
// Many changes dealing with fundamental structure and reading/writing.
//
// Revision 1.19  2006/03/13 19:47:54  sherwin
//
// Fixed bug related to constructor of GeoFac and also makde arguments to GeoFac all consts
//
// Revision 1.18  2006/03/12 14:20:43  sherwin
//
// First compiling version of SpatialDomains and associated modifications
//
// Revision 1.17  2006/03/12 11:06:40  sherwin
//
// First complete copy of code standard code but still not compile tested
//
// Revision 1.16  2006/02/26 21:19:43  bnelson
// Fixed a variety of compiler errors caused by updates to the coding standard.
//
// Revision 1.15  2006/02/19 01:37:34  jfrazier
// Initial attempt at bringing into conformance with the coding standard.  Still more work to be done.  Has not been compiled.
//
//

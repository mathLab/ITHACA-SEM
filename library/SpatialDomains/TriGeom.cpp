////////////////////////////////////////////////////////////////////////////////
//
//  File: TriGeom.cpp
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

#include <SpatialDomains/TriGeom.h>
#include <StdRegions/StdNodalTriExp.h>
#include <LibUtilities/Foundations/ManagerAccess.h>

namespace Nektar
{
    namespace SpatialDomains
    {
        /**
         *
         */
        TriGeom::TriGeom()
        {
            m_geomShapeType = eTriangle;
        }


        /**
         *
         */
        TriGeom::TriGeom(int id, const int coordim):
                                   Geometry2D(coordim), m_fid(id)
        {
            const LibUtilities::BasisKey B0(LibUtilities::eModified_A, 2,
                    LibUtilities::PointsKey(3,LibUtilities::eGaussLobattoLegendre));
            const LibUtilities::BasisKey B1(LibUtilities::eModified_B, 2,
                    LibUtilities::PointsKey(3,LibUtilities::eGaussRadauMAlpha1Beta0));

            m_xmap = Array<OneD, StdRegions::StdExpansion2DSharedPtr>(coordim);

            m_fid = id;

            for(int i = 0; i < m_coordim; ++i)
            {
                m_xmap[i] = MemoryManager<StdRegions::StdTriExp>::AllocateSharedPtr(B0,B1);
            }

        }


        /**
         *
         */
        TriGeom::TriGeom(const int id,
                const VertexComponentSharedPtr verts[],
                const SegGeomSharedPtr edges[],
                const StdRegions::Orientation eorient[]):
                Geometry2D(verts[0]->GetCoordim()),
                m_fid(id)
        {
            m_geomShapeType = eTriangle;

            /// Copy the vert shared pointers.
            m_verts.insert(m_verts.begin(), verts, verts+TriGeom::kNverts);

            /// Copy the edge shared pointers.
            m_edges.insert(m_edges.begin(), edges, edges+TriGeom::kNedges);

            for (int j=0; j<kNedges; ++j)
            {
                m_eorient[j] = eorient[j];
            }

            m_coordim = verts[0]->GetCoordim();
            ASSERTL0(m_coordim > 1,
                    "Cannot call function with dim == 1");

            int order0  = edges[0]->GetBasis(0,0)->GetNumModes();
            int points0 = edges[0]->GetBasis(0,0)->GetNumPoints();
            int order1  = max(order0,
                    max(edges[1]->GetBasis(0,0)->GetNumModes(),
                            edges[2]->GetBasis(0,0)->GetNumModes()));
            int points1 = max(points0,
                    max(edges[1]->GetBasis(0,0)->GetNumPoints(),
                            edges[2]->GetBasis(0,0)->GetNumPoints()));


            const LibUtilities::BasisKey B0(LibUtilities::eModified_A, order0,
                    LibUtilities::PointsKey(points0,LibUtilities::eGaussLobattoLegendre));
            const LibUtilities::BasisKey B1(LibUtilities::eModified_B, order1,
                    LibUtilities::PointsKey(points1,LibUtilities::eGaussRadauMAlpha1Beta0));

            m_xmap = Array<OneD, StdRegions::StdExpansion2DSharedPtr>(m_coordim);

            for(int i = 0; i < m_coordim; ++i)
            {
                m_xmap[i] = MemoryManager<StdRegions::StdTriExp>::AllocateSharedPtr(B0,B1);
            }
        }


        /**
         *
         */
        TriGeom::TriGeom(const int id, const SegGeomSharedPtr edges[],
                const StdRegions::Orientation eorient[]):
                Geometry2D(edges[0]->GetVertex(0)->GetCoordim()),
                m_fid(id)
        {
            m_geomShapeType = eTriangle;

            /// Copy the edge shared pointers.
            m_edges.insert(m_edges.begin(), edges, edges+TriGeom::kNedges);

            for(int j=0; j <kNedges; ++j)
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

            for (int j=0; j<kNedges; ++j)
            {
                m_eorient[j] = eorient[j];
            }

            m_coordim = edges[0]->GetVertex(0)->GetCoordim();
            ASSERTL0(m_coordim > 1,"Cannot call function with dim == 1");

            int order0  = edges[0]->GetBasis(0,0)->GetNumModes();
            int points0 = edges[0]->GetBasis(0,0)->GetNumPoints();
            int order1  = max(order0,
                    max(edges[1]->GetBasis(0,0)->GetNumModes(),
                            edges[2]->GetBasis(0,0)->GetNumModes()));
            int points1 = max(points0,
                    max(edges[1]->GetBasis(0,0)->GetNumPoints(),
                            edges[2]->GetBasis(0,0)->GetNumPoints()));

            const LibUtilities::BasisKey B0(LibUtilities::eModified_A, order0,
                    LibUtilities::PointsKey(points0,LibUtilities::eGaussLobattoLegendre));
            const LibUtilities::BasisKey B1(LibUtilities::eModified_B, order1,
                    LibUtilities::PointsKey(points1,LibUtilities::eGaussRadauMAlpha1Beta0));

            m_xmap = Array<OneD, StdRegions::StdExpansion2DSharedPtr>(m_coordim);

            for(int i = 0; i < m_coordim; ++i)
            {
                m_xmap[i] = MemoryManager<StdRegions::StdTriExp>::AllocateSharedPtr(B0,B1);
            }
        }


        /**
         *
         */
        TriGeom::TriGeom(const int id,
                const SegGeomSharedPtr edges[],
                const StdRegions::Orientation eorient[],
                const CurveSharedPtr &curve) :
                Geometry2D(edges[0]->GetVertex(0)->GetCoordim()),
                m_fid(id)
        {
            m_geomShapeType = eTriangle;

            /// Copy the edge shared pointers.
            m_edges.insert(m_edges.begin(), edges, edges+TriGeom::kNedges);

            for(int j=0; j <kNedges; ++j)
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

            for (int j=0; j<kNedges; ++j)
            {
                m_eorient[j] = eorient[j];
            }

            m_coordim = edges[0]->GetVertex(0)->GetCoordim();
            ASSERTL0(m_coordim > 1,"Cannot call function with dim == 1");

            int order0  = edges[0]->GetBasis(0,0)->GetNumModes();
            int points0 = edges[0]->GetBasis(0,0)->GetNumPoints();
            int order1  = max(order0,
                    max(edges[1]->GetBasis(0,0)->GetNumModes(),
                            edges[2]->GetBasis(0,0)->GetNumModes()));
            int points1 = max(points0,
                    max(edges[1]->GetBasis(0,0)->GetNumPoints(),
                            edges[2]->GetBasis(0,0)->GetNumPoints()));

            const LibUtilities::BasisKey B0(LibUtilities::eModified_A, order0,
                    LibUtilities::PointsKey(points0,LibUtilities::eGaussLobattoLegendre));
            const LibUtilities::BasisKey B1(LibUtilities::eModified_B, order1,
                    LibUtilities::PointsKey(points1,LibUtilities::eGaussRadauMAlpha1Beta0));

            m_xmap = Array<OneD, StdRegions::StdExpansion2DSharedPtr>(m_coordim);

            for(int i = 0; i < m_coordim; ++i)
            {
                m_xmap[i] = MemoryManager<StdRegions::StdTriExp>::AllocateSharedPtr(B0,B1);

                int pdim = LibUtilities::PointsManager()[LibUtilities::PointsKey(2, curve->m_ptype)]->GetPointsDim();

                // Deal with 2D points type separately (e.g. electrostatic or
                // Fekete points).
                if (pdim == 2)
                {
                    int N = curve->m_points.size();
                    int nEdgePts = (-1+(int)sqrt(static_cast<double>(8*N+1)))/2;
                    
                    ASSERTL0(nEdgePts*(nEdgePts+1)/2 == N,
                             "NUMPOINTS must be a triangle number for 2D basis.");
                    
                    for (int j = 0; j < kNedges; ++j)
                    {
                        ASSERTL0(edges[j]->GetXmap(i)->GetNcoeffs() == nEdgePts,
                                 "Number of edge points does not correspond "
                                 "to number of face points.");
                    }
                    
                    // Create a StdNodalTriExp.
                    const LibUtilities::PointsKey P0(
                        nEdgePts, LibUtilities::eGaussLobattoLegendre);
                    const LibUtilities::PointsKey P1(
                        nEdgePts, LibUtilities::eGaussRadauMAlpha1Beta0);
                    const LibUtilities::BasisKey  T0(
                        LibUtilities::eOrtho_A, nEdgePts, P0);
                    const LibUtilities::BasisKey  T1(
                        LibUtilities::eOrtho_B, nEdgePts, P1);
                    
                    StdRegions::StdNodalTriExpSharedPtr t =
                        MemoryManager<StdRegions::StdNodalTriExp>::AllocateSharedPtr(
                            T0,T1,curve->m_ptype);
                    
                    for (int j = 0; j < N; ++j)
                    {
                        t->UpdatePhys()[j] = (curve->m_points[j]->GetPtr())[i];
                    }
                    
                    Array<OneD, NekDouble> tmp(nEdgePts*nEdgePts);
                    t->BwdTrans(t->GetPhys(), tmp);
                    // Interpolate points to standard region.
                    LibUtilities::Interp2D(P0, P1, tmp,
                                           B0.GetPointsKey(),B1.GetPointsKey(),
                                           m_xmap[i]->UpdatePhys());
                    
                    // Forwards transform to get coefficient space.
                    m_xmap[i]->FwdTrans(m_xmap[i]->GetPhys(), m_xmap[i]->UpdateCoeffs());
                }
                else if (pdim == 1)
                {
                    int npts = curve->m_points.size();
                    int nEdgePts = (int)sqrt(static_cast<double>(npts));
                    Array<OneD,NekDouble> tmp(npts);
                    LibUtilities::PointsKey curveKey(nEdgePts, curve->m_ptype);
                    
                    // Sanity checks:
                    // - Curved faces should have square number of points;
                    // - Each edge should have sqrt(npts) points.
                    ASSERTL0(nEdgePts*nEdgePts == npts,
                             "NUMPOINTS should be a square number");
                
                    for (int j = 0; j < kNedges; ++j)
                    {
                        ASSERTL0(edges[j]->GetXmap(i)->GetNcoeffs() == nEdgePts,
                                 "Number of edge points does not correspond "
                                 "to number of face points.");
                    }
                    
                    for (int j = 0; j < npts; ++j)
                    {
                        tmp[j] = (curve->m_points[j]->GetPtr())[i];
                    }
                    
                    // Interpolate curve points to standard triangle points.
                    LibUtilities::Interp2D(curveKey,curveKey,tmp,
                                           B0.GetPointsKey(),B1.GetPointsKey(),
                                           m_xmap[i]->UpdatePhys());
                    
                    // Forwards transform to get coefficient space.
                    m_xmap[i]->FwdTrans(m_xmap[i]->GetPhys(),m_xmap[i]->UpdateCoeffs());
                }
                else
                {
                    ASSERTL0(false, "Only 1D/2D points distributions supported.");
                }
            }
        }


        /**
         *
         */
        TriGeom::TriGeom(const TriGeom &in)
        {
            // From Geomtry
            m_geomShapeType = in.m_geomShapeType;

            // From TriFaceComponent
            m_fid = in.m_fid;
            m_ownVerts = in.m_ownVerts;
            std::list<CompToElmt>::const_iterator def;
            for(def = in.m_elmtMap.begin(); def != in.m_elmtMap.end(); def++)
            {
                m_elmtMap.push_back(*def);
            }

            // From TriGeom
            m_verts = in.m_verts;
            m_edges = in.m_edges;
            for (int i = 0; i < kNedges; i++)
            {
                m_eorient[i] = in.m_eorient[i];
            }
            m_ownData = in.m_ownData;
        }


        /**
         *
         */
        TriGeom::~TriGeom()
        {
        }


        /**
         * Given local collapsed coordinate Lcoord return the value of physical
         * coordinate in direction i
         */
        NekDouble TriGeom::GetCoord(const int i,
                const Array<OneD, const NekDouble> &Lcoord)
        {
            ASSERTL1(m_state == ePtsFilled,
                    "Goemetry is not in physical space");
            return m_xmap[i]->PhysEvaluate(Lcoord);
        }


        /**
         * TODO: implement eight different case of face orientation
         */
        StdRegions::Orientation TriGeom::GetFaceOrientation(
                const TriGeom &face1,
                const TriGeom &face2)
        {
            StdRegions::Orientation returnval = StdRegions::eDir1FwdDir1_Dir2FwdDir2;
            ASSERTL0(false,"this function is not yet implemented.");
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


        /**
         *
         */
        void TriGeom::v_AddElmtConnected(int gvo_id, int locid)
        {
            CompToElmt ee(gvo_id,locid);
            m_elmtMap.push_back(ee);
        }


        /**
         *
         */
        int  TriGeom::v_NumElmtConnected() const
        {
            return int(m_elmtMap.size());
        }


        /**
         *
         */
        bool TriGeom::v_IsElmtConnected(int gvo_id, int locid) const
        {
            std::list<CompToElmt>::const_iterator def;
            CompToElmt ee(gvo_id,locid);

            def = find(m_elmtMap.begin(),m_elmtMap.end(),ee);

            // Found the element connectivity object in the list
            return (def != m_elmtMap.end());
        }


        /**
         *
         */
        int  TriGeom::v_GetFid() const
        {
            return m_fid;
        }


        /**
         *
         */
        int  TriGeom::v_GetCoordim() const
        {
            return m_coordim;
        }


        /**
         *
         */
        const LibUtilities::BasisSharedPtr TriGeom::v_GetBasis(const int i, const int j)
        {
            return m_xmap[i]->GetBasis(j);
        }


        /**
         *
         */
        const LibUtilities::BasisSharedPtr TriGeom::v_GetEdgeBasis(const int i, const int j)
        {
            ASSERTL1(j <=2,"edge is out of range");
            if(j == 0)
            {
                return m_xmap[i]->GetBasis(0);
            }
            else
            {
                return m_xmap[i]->GetBasis(1);
            }
        }


        /**
         *
         */
        Array<OneD,NekDouble> &TriGeom::v_UpdatePhys(const int i)
        {
            return m_xmap[i]->UpdatePhys();
        }


        /**
         *
         */
        NekDouble TriGeom::v_GetCoord(const int i, const Array<OneD,const NekDouble> &Lcoord)
        {
            return GetCoord(i,Lcoord);
        }


        /**
         * Set up GeoFac for this geometry using Coord quadrature distribution
         */
        void TriGeom::v_GenGeomFactors(const Array<OneD, const LibUtilities::BasisSharedPtr> &tbasis)
        {
            GeomType Gtype = eRegular;

            TriGeom::v_FillGeom();

            // check to see if expansions are linear
            for(int i = 0; i < m_coordim; ++i)
            {
                if((m_xmap[i]->GetBasisNumModes(0) != 2)||
                        (m_xmap[i]->GetBasisNumModes(1) != 2))
                {
                    Gtype = eDeformed;
                }
            }

            m_geomFactors = MemoryManager<GeomFactors2D>::AllocateSharedPtr(Gtype, m_coordim, m_xmap, tbasis);
        }


        /**
         *
         */
        void TriGeom::v_SetOwnData()
        {
            m_ownData = true;
        }


        /**
         * Note verts and edges are listed according to anticlockwise
         * convention but points in _coeffs have to be in array format from
         * left to right.
         */
        void TriGeom::v_FillGeom()
        {
            // check to see if geometry structure is already filled
            if(m_state != ePtsFilled)
            {
                int i,j,k;
                int nEdgeCoeffs = m_xmap[0]->GetEdgeNcoeffs(0);

                Array<OneD, unsigned int> mapArray (nEdgeCoeffs);
                Array<OneD, int>          signArray(nEdgeCoeffs);

                for(i = 0; i < kNedges; i++)
                {
                    m_edges[i]->FillGeom();
                    m_xmap[0]->GetEdgeToElementMap(i,m_eorient[i],mapArray,signArray);

                    //nEdgeCoeffs = m_xmap[0]->GetEdgeNcoeffs(i);
                    nEdgeCoeffs = (*m_edges[i])[0]->GetNcoeffs();

                    for(j = 0 ; j < m_coordim; j++)
                    {
                        for(k = 0; k < nEdgeCoeffs; k++)
                        {
                            (m_xmap[j]->UpdateCoeffs())[ mapArray[k] ]
                                                         = signArray[k]*
                                                         ((*m_edges[i])[j]->GetCoeffs())[k];
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


        /**
         *
         */
        void TriGeom::v_GetLocCoords(const Array<OneD,const NekDouble> &coords, Array<OneD,NekDouble> &Lcoords)
        {
            TriGeom::v_FillGeom();

            // calculate local coordinate for coord
            if(GetGtype() == eRegular)
            { // can assume it is right angled rectangle
                VertexComponent dv1, dv2, norm, orth1, orth2;
                VertexComponent *xin;

                switch(m_coordim)
                {
                case 2:
                    xin = new VertexComponent (m_coordim,0,coords[0], coords[1], 0.0);
                    break;
                case 3:
                    xin = new VertexComponent (m_coordim,0,coords[0], coords[1], coords[2]);
                    break;
                }

                dv1.Sub(*m_verts[1],*m_verts[0]);
                dv2.Sub(*m_verts[2],*m_verts[0]);

                norm.Mult(dv1,dv2);

                orth1.Mult(norm,dv1);
                orth2.Mult(norm,dv2);

                xin[0] *= 2.0;
                xin[0] -= *m_verts[1];
                xin[0] -= *m_verts[2];

                Lcoords[0] = xin->dot(orth2)/dv1.dot(orth2);
                Lcoords[1] = xin->dot(orth1)/dv2.dot(orth1);

                delete xin;
            }
            else
            {
                NEKERROR(ErrorUtil::efatal,
                        "inverse mapping must be set up to use this call");
            }
        }


        /**
         *
         */
        int TriGeom::v_GetEid(int i) const
        {
            ASSERTL2((i >=0) && (i <= 2),"Edge id must be between 0 and 2");
            return m_edges[i]->GetEid();
        }


        /**
         *
         */
        int TriGeom::v_GetVid(int i) const
        {
            ASSERTL2((i >=0) && (i <= 2),"Vertex id must be between 0 and 2");
            return m_verts[i]->GetVid();
        }


        /**
         *
         */
        const VertexComponentSharedPtr TriGeom::v_GetVertex(int i) const
        {
            ASSERTL2((i >=0) && (i <= 2),"Vertex id must be between 0 and 2");
            return m_verts[i];
        }


        /**
         *
         */
        const Geometry1DSharedPtr TriGeom::v_GetEdge(int i) const
        {
            ASSERTL2((i >=0) && (i <= 2),"Edge id must be between 0 and 3");
            return m_edges[i];
        }


        /**
         *
         */
        StdRegions::Orientation TriGeom::v_GetEorient(const int i) const
        {
            ASSERTL2((i >=0) && (i <= 2),"Edge id must be between 0 and 2");
            return m_eorient[i];
        }


        /**
         *
         */
        StdRegions::Orientation TriGeom::v_GetCartesianEorient(const int i) const
        {
            ASSERTL2((i >=0) && (i <= 3),"Edge id must be between 0 and 3");
            if(i < 2)
            {
                return m_eorient[i];
            }
            else
            {
                if(m_eorient[i] == StdRegions::eForwards)
                {
                    return StdRegions::eBackwards;
                }
                else
                {
                    return StdRegions::eForwards;
                }
            }
        }


        /**
         *
         */
        int TriGeom::v_WhichEdge(SegGeomSharedPtr edge)
        {
            int returnval = -1;

            SegGeomVector::iterator edgeIter;
            int i;

            for (i=0,edgeIter = m_edges.begin(); edgeIter != m_edges.end(); ++edgeIter,++i)
            {
                if (*edgeIter == edge)
                {
                    returnval = i;
                    break;
                }
            }

            return returnval;
        }


        /**
         *
         */
        int TriGeom::v_GetNumVerts() const
        {
            return kNverts;
        }


        /**
         *
         */
        int TriGeom::v_GetNumEdges() const
        {
            return kNedges;
        }


        /**
         *
         */
        bool TriGeom::v_ContainsPoint(const Array<OneD, const NekDouble> &gloCoord, NekDouble tol)
        {
            ASSERTL1(gloCoord.num_elements() >= 2,
                    "Two dimensional geometry expects at least two coordinates.");

            Array<OneD,NekDouble> stdCoord(GetCoordim(),0.0);
            GetLocCoords(gloCoord, stdCoord);
            if (stdCoord[0] >= -(1+tol) && stdCoord[1] >= -(1+tol)
                    && stdCoord[0] + stdCoord[1] <= tol)
            {
                return true;
            }
            return false;
        }
    }; //end of namespace
}; //end of namespace

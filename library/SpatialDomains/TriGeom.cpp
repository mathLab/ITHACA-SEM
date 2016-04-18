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

#include <SpatialDomains/TriGeom.h>
#include <StdRegions/StdNodalTriExp.h>
#include <LibUtilities/Foundations/ManagerAccess.h>
#include <LibUtilities/Foundations/Interp.h>

#include <SpatialDomains/SegGeom.h>
#include <SpatialDomains/GeomFactors.h>

using namespace std;

namespace Nektar
{
    namespace SpatialDomains
    {
        /**
         *
         */
        TriGeom::TriGeom()
        {
            m_shapeType = LibUtilities::eTriangle;
        }

        /**
         *
         */
        TriGeom::TriGeom(int id, const int coordim):
                                   Geometry2D(coordim), m_fid(id)
        {
            m_globalID = m_fid;
            SetUpXmap();
            SetUpCoeffs(m_xmap->GetNcoeffs());
        }

        /**
         *
         */
        TriGeom::TriGeom(const int id,
                const PointGeomSharedPtr verts[],
                const SegGeomSharedPtr edges[],
                const StdRegions::Orientation eorient[]):
                Geometry2D(verts[0]->GetCoordim()),
                m_fid(id)
        {
            m_globalID = m_fid;
            m_shapeType = LibUtilities::eTriangle;

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

            SetUpXmap();
            SetUpCoeffs(m_xmap->GetNcoeffs());
        }


        /**
         *
         */
        TriGeom::TriGeom(const int id, const SegGeomSharedPtr edges[],
                const StdRegions::Orientation eorient[]):
                Geometry2D(edges[0]->GetVertex(0)->GetCoordim()),
                m_fid(id)
        {
            m_globalID = m_fid;
            m_shapeType = LibUtilities::eTriangle;

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

            SetUpXmap();
            SetUpCoeffs(m_xmap->GetNcoeffs());
        }


        /**
         *
         */
        TriGeom::TriGeom(
            const int                      id,
            const SegGeomSharedPtr         edges[],
            const StdRegions::Orientation  eorient[],
            const CurveSharedPtr          &curve)
            : Geometry2D(edges[0]->GetVertex(0)->GetCoordim()),
              m_fid(id),
              m_curve(curve)
        {
            m_globalID = m_fid;
            m_shapeType =  LibUtilities::eTriangle;

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

            SetUpXmap();
            SetUpCoeffs(m_xmap->GetNcoeffs());
        }

        /**
         *
         */
        TriGeom::TriGeom(const TriGeom &in)
        {
            // From Geomtry
            m_shapeType = in.m_shapeType;

            // From TriFaceComponent
            m_fid = in.m_fid;
            m_globalID = m_fid;
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
                "Geometry is not in physical space");

            Array<OneD, NekDouble> tmp(m_xmap->GetTotPoints());
            m_xmap->BwdTrans(m_coeffs[i], tmp);

            return m_xmap->PhysEvaluate(Lcoord, tmp);
        }

        StdRegions::Orientation TriGeom::GetFaceOrientation(
            const TriGeom &face1,
            const TriGeom &face2)
        {
            return GetFaceOrientation(face1.m_verts, face2.m_verts);
        }

        StdRegions::Orientation TriGeom::GetFaceOrientation(
            const PointGeomVector &face1,
            const PointGeomVector &face2)
        {
            int i, j, vmap[3] = {-1,-1,-1};
            NekDouble x, y, z, x1, y1, z1, cx = 0.0, cy = 0.0, cz = 0.0;

            // For periodic faces, we calculate the vector between the centre
            // points of the two faces. (For connected faces this will be
            // zero). We can then use this to determine alignment later in the
            // algorithm.
            for (i = 0; i < 3; ++i)
            {
                cx += (*face2[i])(0) - (*face1[i])(0);
                cy += (*face2[i])(1) - (*face1[i])(1);
                cz += (*face2[i])(2) - (*face1[i])(2);
            }
            cx /= 3;
            cy /= 3;
            cz /= 3;

            // Now construct a mapping which takes us from the vertices of one
            // face to the other. That is, vertex j of face2 corresponds to
            // vertex vmap[j] of face1.
            for (i = 0; i < 3; ++i)
            {
                x = (*face1[i])(0);
                y = (*face1[i])(1);
                z = (*face1[i])(2);
                for (j = 0; j < 3; ++j)
                {
                    x1 = (*face2[j])(0)-cx;
                    y1 = (*face2[j])(1)-cy;
                    z1 = (*face2[j])(2)-cz;
                    if (sqrt((x1-x)*(x1-x)+(y1-y)*(y1-y)+(z1-z)*(z1-z)) < 1e-8)
                    {
                        vmap[j] = i;
                        break;
                    }
                }
            }

            if (vmap[1] == (vmap[0]+1) % 3)
            {
                switch (vmap[0])
                {
                    case 0:
                        return StdRegions::eDir1FwdDir1_Dir2FwdDir2;
                        break;
                    case 1:
                        return StdRegions::eDir1FwdDir2_Dir2BwdDir1;
                        break;
                    case 2:
                        return StdRegions::eDir1BwdDir1_Dir2BwdDir2;
                        break;
                }
            }
            else
            {
                switch (vmap[0])
                {
                    case 0:
                        return StdRegions::eDir1FwdDir2_Dir2FwdDir1;
                        break;
                    case 1:
                        return StdRegions::eDir1BwdDir1_Dir2FwdDir2;
                        break;
                    case 2:
                        return StdRegions::eDir1BwdDir2_Dir2BwdDir1;
                        break;
                }
            }

            ASSERTL0(false, "Unable to determine triangle orientation");
            return StdRegions::eNoOrientation;
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
        const LibUtilities::BasisSharedPtr TriGeom::v_GetBasis(const int i)
        {
            return m_xmap->GetBasis(i);
        }


        /**
         *
         */
        const LibUtilities::BasisSharedPtr TriGeom::v_GetEdgeBasis(const int i)
        {
            ASSERTL1(i <= 2, "edge is out of range");

            if(i == 0)
            {
                return m_xmap->GetBasis(0);
            }
            else
            {
                return m_xmap->GetBasis(1);
            }
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
        void TriGeom::v_GenGeomFactors()
        {
            if (m_geomFactorsState != ePtsFilled)
            {
                GeomType Gtype = eRegular;

                TriGeom::v_FillGeom();

                // check to see if expansions are linear
                for(int i = 0; i < m_coordim; ++i)
                {
                    if(m_xmap->GetBasisNumModes(0) != 2 ||
                       m_xmap->GetBasisNumModes(1) != 2)
                    {
                        Gtype = eDeformed;
                    }
                }

                m_geomFactors = MemoryManager<GeomFactors>::AllocateSharedPtr(
                    Gtype, m_coordim, m_xmap, m_coeffs);

                m_geomFactorsState = ePtsFilled;
            }
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
                int nEdgeCoeffs = m_xmap->GetEdgeNcoeffs(0);

                if (m_curve)
                {
                    int pdim = LibUtilities::PointsManager()[
                        LibUtilities::PointsKey(2, m_curve->m_ptype)]
                        ->GetPointsDim();

                    // Deal with 2D points type separately
                    // (e.g. electrostatic or Fekete points) to 1D tensor
                    // product.
                    if (pdim == 2)
                    {
                        int N = m_curve->m_points.size();
                        int nEdgePts = (
                            -1+(int)sqrt(static_cast<NekDouble>(8*N+1)))/2;

                        ASSERTL0(nEdgePts*(nEdgePts+1)/2 == N,
                                 "NUMPOINTS should be a triangle number for"
                                 " triangle "
                                 + boost::lexical_cast<string>(m_globalID));

                        for (i = 0; i < kNedges; ++i)
                        {
                            ASSERTL0(
                                m_edges[i]->GetXmap()->GetNcoeffs() == nEdgePts,
                                "Number of edge points does not correspond "
                                "to number of face points in triangle "
                                + boost::lexical_cast<string>(m_globalID));
                        }

                        const LibUtilities::PointsKey P0(
                            nEdgePts, LibUtilities::eGaussLobattoLegendre);
                        const LibUtilities::PointsKey P1(
                            nEdgePts, LibUtilities::eGaussRadauMAlpha1Beta0);
                        const LibUtilities::BasisKey  T0(
                            LibUtilities::eOrtho_A, nEdgePts, P0);
                        const LibUtilities::BasisKey  T1(
                            LibUtilities::eOrtho_B, nEdgePts, P1);
                        Array<OneD, NekDouble> phys(
                            max(nEdgePts*nEdgePts, m_xmap->GetTotPoints()));
                        Array<OneD, NekDouble> tmp(nEdgePts*nEdgePts);

                        for(i = 0; i < m_coordim; ++i)
                        {
                            // Create a StdNodalTriExp.
                            StdRegions::StdNodalTriExpSharedPtr t =
                                MemoryManager<StdRegions::StdNodalTriExp>
                                ::AllocateSharedPtr(T0, T1, m_curve->m_ptype);

                            for (j = 0; j < N; ++j)
                            {
                                phys[j] = (m_curve->m_points[j]->GetPtr())[i];
                            }

                            t->BwdTrans(phys, tmp);

                            // Interpolate points to standard region.
                            LibUtilities::Interp2D(
                                P0, P1, tmp,
                                m_xmap->GetBasis(0)->GetPointsKey(),
                                m_xmap->GetBasis(1)->GetPointsKey(),
                                phys);

                            // Forwards transform to get coefficient space.
                            m_xmap->FwdTrans(phys, m_coeffs[i]);
                        }
                    }
                    else if (pdim == 1)
                    {
                        int npts = m_curve->m_points.size();
                        int nEdgePts = (int)sqrt(static_cast<NekDouble>(npts));
                        Array<OneD, NekDouble> tmp (npts);
                        Array<OneD, NekDouble> phys(m_xmap->GetTotPoints());
                        LibUtilities::PointsKey curveKey(
                            nEdgePts, m_curve->m_ptype);

                        // Sanity checks:
                        // - Curved faces should have square number of points;
                        // - Each edge should have sqrt(npts) points.
                        ASSERTL0(nEdgePts * nEdgePts == npts,
                                 "NUMPOINTS should be a square number for"
                                 " triangle "
                                 + boost::lexical_cast<string>(m_globalID));

                        for (i = 0; i < kNedges; ++i)
                        {
                            ASSERTL0(
                                m_edges[i]->GetXmap()->GetNcoeffs() == nEdgePts,
                                "Number of edge points does not correspond to "
                                "number of face points in triangle "
                                + boost::lexical_cast<string>(m_globalID));
                        }

                        for (i = 0; i < m_coordim; ++i)
                        {
                            for (j = 0; j < npts; ++j)
                            {
                                tmp[j] = (m_curve->m_points[j]->GetPtr())[i];
                            }

                            // Interpolate curve points to standard triangle
                            // points.
                            LibUtilities::Interp2D(
                                curveKey, curveKey, tmp,
                                m_xmap->GetBasis(0)->GetPointsKey(),
                                m_xmap->GetBasis(1)->GetPointsKey(),
                                phys);

                            // Forwards transform to get coefficient space.
                            m_xmap->FwdTrans(phys, m_coeffs[i]);
                        }
                    }
                    else
                    {
                        ASSERTL0(false, "Only 1D/2D points distributions "
                                        "supported.");
                    }
                }

                Array<OneD, unsigned int> mapArray (nEdgeCoeffs);
                Array<OneD, int>          signArray(nEdgeCoeffs);

                for(i = 0; i < kNedges; i++)
                {
                    m_edges[i]->FillGeom();
                    m_xmap->GetEdgeToElementMap(i,m_eorient[i],mapArray,signArray);

                    nEdgeCoeffs = m_edges[i]->GetXmap()->GetNcoeffs();

                    for(j = 0 ; j < m_coordim; j++)
                    {
                        for(k = 0; k < nEdgeCoeffs; k++)
                        {
                            m_coeffs[j][mapArray[k]] =
                                signArray[k] * m_edges[i]->GetCoeffs(j)[k];
                        }
                    }
                }

                m_state = ePtsFilled;
            }
        }


        /**
         *
         */
        NekDouble TriGeom::v_GetLocCoords(const Array<OneD,const NekDouble> &coords,
                                          Array<OneD,NekDouble> &Lcoords)
        {
            NekDouble resid = 0.0;
            TriGeom::v_FillGeom();

            // calculate local coordinate for coord
            if(GetMetricInfo()->GetGtype() == eRegular)
            {
                NekDouble coords2 = (m_coordim == 3)? coords[2]: 0.0;
                PointGeom dv1, dv2, norm, orth1, orth2;
                PointGeom xin(m_coordim,0,coords[0],coords[1],coords2);

                // Calculate edge vectors from 0-1 and 0-2 edges.
                dv1.Sub(*m_verts[1],*m_verts[0]);
                dv2.Sub(*m_verts[2],*m_verts[0]);

                // Obtain normal to plane in which dv1 and dv2 lie
                norm.Mult(dv1,dv2);

                // Obtain vector which are proportional to normal of dv1 and dv2.
                orth1.Mult(norm,dv1);
                orth2.Mult(norm,dv2);

                // Start with vector of desired points minus vertex_0
                xin -= *m_verts[0];

                // Calculate length using L/|dv1| = (x-v0).n1/(dv1.n1) for coordiante 1
                // Then rescale to [-1,1].
                Lcoords[0] = xin.dot(orth2)/dv1.dot(orth2);
                Lcoords[0] = 2*Lcoords[0]-1;
                Lcoords[1] = xin.dot(orth1)/dv2.dot(orth1);
                Lcoords[1] = 2*Lcoords[1]-1;
            }
            else
            {
                // Determine nearest point of coords  to values in m_xmap
                int npts = m_xmap->GetTotPoints();
                Array<OneD, NekDouble> ptsx(npts), ptsy(npts);
                Array<OneD, NekDouble> tmpx(npts), tmpy(npts);

                m_xmap->BwdTrans(m_coeffs[0], ptsx);
                m_xmap->BwdTrans(m_coeffs[1], ptsy);

                const Array<OneD, const NekDouble> za = m_xmap->GetPoints(0);
                const Array<OneD, const NekDouble> zb = m_xmap->GetPoints(1);

                //guess the first local coords based on nearest point
                Vmath::Sadd(npts, -coords[0], ptsx,1,tmpx,1);
                Vmath::Sadd(npts, -coords[1], ptsy,1,tmpy,1);
                Vmath::Vmul (npts, tmpx,1,tmpx,1,tmpx,1);
                Vmath::Vvtvp(npts, tmpy,1,tmpy,1,tmpx,1,tmpx,1);

                int min_i = Vmath::Imin(npts,tmpx,1);

                Lcoords[0] = za[min_i%za.num_elements()];
                Lcoords[1] = zb[min_i/za.num_elements()];

                // recover cartesian coordinate from collapsed coordinate.
                Lcoords[0] = (1.0+Lcoords[0])*(1.0-Lcoords[1])/2 -1.0;

                // Perform newton iteration to find local coordinates
                NewtonIterationForLocCoord(coords, ptsx, ptsy, Lcoords,resid);
            }
            return resid;
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
        PointGeomSharedPtr TriGeom::v_GetVertex(int i) const
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
         * @brief Determines if a point specified in global coordinates is
         * located within this tetrahedral geometry.
         */
        bool TriGeom::v_ContainsPoint(
            const Array<OneD, const NekDouble> &gloCoord, NekDouble tol)
        {
            Array<OneD,NekDouble> locCoord(GetCoordim(),0.0);
            return v_ContainsPoint(gloCoord,locCoord,tol);

        }

        /**
         *
         */
        bool TriGeom::v_ContainsPoint(const Array<OneD, const NekDouble> &gloCoord,
                                      Array<OneD, NekDouble> &stdCoord,
                                      NekDouble tol)
        {
            NekDouble resid;
            return v_ContainsPoint(gloCoord,stdCoord,tol,resid);
        }

        bool TriGeom::v_ContainsPoint(const Array<OneD, const NekDouble> &gloCoord,
                                      Array<OneD, NekDouble> &stdCoord,
                                      NekDouble tol,
                                      NekDouble &resid)
        {
            ASSERTL1(gloCoord.num_elements() >= 2,
                    "Two dimensional geometry expects at least two coordinates.");

            resid = GetLocCoords(gloCoord, stdCoord);
            if (stdCoord[0] >= -(1+tol) && stdCoord[1] >= -(1+tol)
                    && stdCoord[0] + stdCoord[1] <= tol)
            {
                return true;
            }
            return false;
        }

        void TriGeom::v_Reset(
            CurveMap &curvedEdges,
            CurveMap &curvedFaces)
        {
            Geometry::v_Reset(curvedEdges, curvedFaces);
            CurveMap::iterator it = curvedFaces.find(m_globalID);

            if (it != curvedFaces.end())
            {
                m_curve = it->second;
            }

            for (int i = 0; i < 3; ++i)
            {
                m_edges[i]->Reset(curvedEdges, curvedFaces);
            }

            SetUpXmap();
            SetUpCoeffs(m_xmap->GetNcoeffs());
        }

        void TriGeom::SetUpXmap()
        {
            int order0  = m_edges[0]->GetBasis(0)->GetNumModes();
            int order1  = max(order0,
                              max(m_edges[1]->GetBasis(0)->GetNumModes(),
                                  m_edges[2]->GetBasis(0)->GetNumModes()));

            const LibUtilities::BasisKey B0(
                LibUtilities::eModified_A, order0, LibUtilities::PointsKey(
                    order0+1, LibUtilities::eGaussLobattoLegendre));
            const LibUtilities::BasisKey B1(
                LibUtilities::eModified_B, order1, LibUtilities::PointsKey(
                    order1, LibUtilities::eGaussRadauMAlpha1Beta0));

            m_xmap = MemoryManager<StdRegions::StdTriExp>
                ::AllocateSharedPtr(B0,B1);
        }
    }; //end of namespace
}; //end of namespace

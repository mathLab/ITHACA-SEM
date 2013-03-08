////////////////////////////////////////////////////////////////////////////////
//
//  File: HexGeom.cpp
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
//  Description: Hexahedral geometry definition.
//
////////////////////////////////////////////////////////////////////////////////

#include <SpatialDomains/HexGeom.h>
#include <SpatialDomains/Geometry1D.h>
#include <StdRegions/StdHexExp.h>
#include <SpatialDomains/SegGeom.h>
#include <SpatialDomains/GeomFactors3D.h>

namespace Nektar
{
    namespace SpatialDomains
    {
        HexGeom::HexGeom()
        {
            m_shapeType = LibUtilities::eHexahedron;
        }

        HexGeom::HexGeom(const QuadGeomSharedPtr faces[]):
            Geometry3D(faces[0]->GetEdge(0)->GetVertex(0)->GetCoordim())
        {
            m_shapeType = LibUtilities::eHexahedron;

            /// Copy the face shared pointers
            m_faces.insert(m_faces.begin(), faces, faces+HexGeom::kNfaces);

            /// Set up orientation vectors with correct amount of elements.
            m_eorient.resize(kNedges);
            m_forient.resize(kNfaces);

            SetUpLocalEdges();
            SetUpLocalVertices();
            SetUpEdgeOrientation();
            SetUpFaceOrientation();

            /// Determine necessary order for standard region. This can almost
            /// certainly be simplified but works for now!
            vector<int> tmp1, tmp2;

            if (m_forient[0] < 9)
            {
                tmp1.push_back(faces[0]->GetXmap(0)->GetEdgeNcoeffs  (0));
                tmp1.push_back(faces[0]->GetXmap(0)->GetEdgeNcoeffs  (2));
                tmp2.push_back(faces[0]->GetXmap(0)->GetEdgeNumPoints(0));
                tmp2.push_back(faces[0]->GetXmap(0)->GetEdgeNumPoints(2));
            }
            else
            {
                tmp1.push_back(faces[0]->GetXmap(0)->GetEdgeNcoeffs  (1));
                tmp1.push_back(faces[0]->GetXmap(0)->GetEdgeNcoeffs  (3));
                tmp2.push_back(faces[0]->GetXmap(0)->GetEdgeNumPoints(1));
                tmp2.push_back(faces[0]->GetXmap(0)->GetEdgeNumPoints(3));
            }
            
            if (m_forient[5] < 9)
            {
                tmp1.push_back(faces[5]->GetXmap(0)->GetEdgeNcoeffs  (0));
                tmp1.push_back(faces[5]->GetXmap(0)->GetEdgeNcoeffs  (2));
                tmp2.push_back(faces[5]->GetXmap(0)->GetEdgeNumPoints(0));
                tmp2.push_back(faces[5]->GetXmap(0)->GetEdgeNumPoints(2));
            }
            else
            {
                tmp1.push_back(faces[5]->GetXmap(0)->GetEdgeNcoeffs  (1));
                tmp1.push_back(faces[5]->GetXmap(0)->GetEdgeNcoeffs  (3));
                tmp2.push_back(faces[5]->GetXmap(0)->GetEdgeNumPoints(1));
                tmp2.push_back(faces[5]->GetXmap(0)->GetEdgeNumPoints(3));
            }
            
            int order0  = *max_element(tmp1.begin(), tmp1.end());
            int points0 = *max_element(tmp2.begin(), tmp2.end());
            
            tmp1.clear();
            tmp2.clear();
            
            if (m_forient[0] < 9)
            {
                tmp1.push_back(faces[0]->GetXmap(0)->GetEdgeNcoeffs  (1));
                tmp1.push_back(faces[0]->GetXmap(0)->GetEdgeNcoeffs  (3));
                tmp2.push_back(faces[0]->GetXmap(0)->GetEdgeNumPoints(1));
                tmp2.push_back(faces[0]->GetXmap(0)->GetEdgeNumPoints(3));
            }
            else
            {
                tmp1.push_back(faces[0]->GetXmap(0)->GetEdgeNcoeffs  (0));
                tmp1.push_back(faces[0]->GetXmap(0)->GetEdgeNcoeffs  (2));
                tmp2.push_back(faces[0]->GetXmap(0)->GetEdgeNumPoints(0));
                tmp2.push_back(faces[0]->GetXmap(0)->GetEdgeNumPoints(2));
            }

            if (m_forient[5] < 9)
            {
                tmp1.push_back(faces[5]->GetXmap(0)->GetEdgeNcoeffs  (1));
                tmp1.push_back(faces[5]->GetXmap(0)->GetEdgeNcoeffs  (3));
                tmp2.push_back(faces[5]->GetXmap(0)->GetEdgeNumPoints(1));
                tmp2.push_back(faces[5]->GetXmap(0)->GetEdgeNumPoints(3));
            }
            else
            {
                tmp1.push_back(faces[5]->GetXmap(0)->GetEdgeNcoeffs  (0));
                tmp1.push_back(faces[5]->GetXmap(0)->GetEdgeNcoeffs  (2));
                tmp2.push_back(faces[5]->GetXmap(0)->GetEdgeNumPoints(0));
                tmp2.push_back(faces[5]->GetXmap(0)->GetEdgeNumPoints(2));
            }
            
            int order1  = *max_element(tmp1.begin(), tmp1.end());
            int points1 = *max_element(tmp2.begin(), tmp2.end());
            
            tmp1.clear();
            tmp2.clear();

            if (m_forient[1] < 9)
            {
                tmp1.push_back(faces[1]->GetXmap(0)->GetEdgeNcoeffs  (1));
                tmp1.push_back(faces[1]->GetXmap(0)->GetEdgeNcoeffs  (3));
                tmp2.push_back(faces[1]->GetXmap(0)->GetEdgeNumPoints(1));
                tmp2.push_back(faces[1]->GetXmap(0)->GetEdgeNumPoints(3));
            }
            else
            {
                tmp1.push_back(faces[1]->GetXmap(0)->GetEdgeNcoeffs  (0));
                tmp1.push_back(faces[1]->GetXmap(0)->GetEdgeNcoeffs  (2));
                tmp2.push_back(faces[1]->GetXmap(0)->GetEdgeNumPoints(0));
                tmp2.push_back(faces[1]->GetXmap(0)->GetEdgeNumPoints(2));
            }
            
            if (m_forient[3] < 9)
            {
                tmp1.push_back(faces[3]->GetXmap(0)->GetEdgeNcoeffs  (1));
                tmp1.push_back(faces[3]->GetXmap(0)->GetEdgeNcoeffs  (3));
                tmp2.push_back(faces[3]->GetXmap(0)->GetEdgeNumPoints(1));
                tmp2.push_back(faces[3]->GetXmap(0)->GetEdgeNumPoints(3));
            }
            else
            {
                tmp1.push_back(faces[3]->GetXmap(0)->GetEdgeNcoeffs  (0));
                tmp1.push_back(faces[3]->GetXmap(0)->GetEdgeNcoeffs  (2));
                tmp2.push_back(faces[3]->GetXmap(0)->GetEdgeNumPoints(0));
                tmp2.push_back(faces[3]->GetXmap(0)->GetEdgeNumPoints(2));
            }

            int order2  = *max_element(tmp1.begin(), tmp1.end());
            int points2 = *max_element(tmp2.begin(), tmp2.end());

            const LibUtilities::BasisKey A(
                LibUtilities::eModified_A, order0,
                LibUtilities::PointsKey(points0,LibUtilities::eGaussLobattoLegendre));
            const LibUtilities::BasisKey B(
                LibUtilities::eModified_A, order1,
                LibUtilities::PointsKey(points1,LibUtilities::eGaussLobattoLegendre));
            const LibUtilities::BasisKey C(
                LibUtilities::eModified_A, order2,
                LibUtilities::PointsKey(points2,LibUtilities::eGaussLobattoLegendre));

            m_xmap = Array<OneD, StdRegions::StdExpansion3DSharedPtr>(m_coordim);

            for(int i = 0; i < m_coordim; ++i)
            {
                m_xmap[i] = MemoryManager<StdRegions::StdHexExp>::AllocateSharedPtr(A,B,C);
            }
        }

        /*
        HexGeom::HexGeom(const QuadGeomSharedPtr faces[], 
                         const Array<OneD, StdRegions::StdExpansion3DSharedPtr> &xMap) :
            Geometry3D(faces[0]->GetEdge(0)->GetVertex(0)->GetCoordim())
        {
            HexGeom::HexGeom(faces);
            
            for(int i = 0; i < xMap.num_elements(); ++i)
            {
                m_xmap[i] = xMap[i];
            }
        }
        */

        HexGeom::~HexGeom()
        {
            
        }

        void HexGeom::v_GenGeomFactors(
                const Array<OneD, const LibUtilities::BasisSharedPtr> &tbasis)
        {
            if (m_geomFactorsState != ePtsFilled)
            {
                int i,f;
                GeomType Gtype = eRegular;

                v_FillGeom();

                // check to see if expansions are linear
                for(i = 0; i < m_coordim; ++i)
                {
                    if (m_xmap[i]->GetBasisNumModes(0) != 2 ||
                        m_xmap[i]->GetBasisNumModes(1) != 2 ||
                        m_xmap[i]->GetBasisNumModes(2) != 2 )
                    {
                        Gtype = eDeformed;
                    }
                }

                // check to see if all faces are parallelograms
                if(Gtype == eRegular)
                {
                    const unsigned int faceVerts[kNfaces][QuadGeom::kNverts] =
                        { {0,1,2,3} ,
                          {0,1,5,4} ,
                          {1,2,6,5} ,
                          {3,2,6,7} ,
                          {0,3,7,4} ,
                          {4,5,6,7} };

                    for(f = 0; f < kNfaces; f++)
                    {
                        // Ensure each face is a parallelogram? Check this.
                        for (i = 0; i < m_coordim; i++)
                        {
                            if( fabs( (*m_verts[ faceVerts[f][0] ])(i) -
                                      (*m_verts[ faceVerts[f][1] ])(i) +
                                      (*m_verts[ faceVerts[f][2] ])(i) -
                                      (*m_verts[ faceVerts[f][3] ])(i) )
                                    > NekConstants::kNekZeroTol )
                            {
                                Gtype = eDeformed;
                                break;
                            }
                        }

                        if (Gtype == eDeformed)
                        {
                            break;
                        }
                    }
                }

                m_geomFactors = MemoryManager<GeomFactors3D>::AllocateSharedPtr(
                    Gtype, m_coordim, m_xmap, tbasis, true);

                m_geomFactorsState = ePtsFilled;
            }
        }

        void HexGeom::v_GetLocCoords(
            const Array<OneD, const NekDouble> &coords, 
                  Array<OneD,       NekDouble> &Lcoords)
        {
            int i;

            v_FillGeom();

            // calculate local coordinate for coord
            if(GetGtype() == eRegular)
            {   
                NekDouble len0 = 0.0 ;
                NekDouble len1 = 0.0;
                NekDouble len2 = 0.0;
                NekDouble xi0 = 0.0;
                NekDouble xi1 = 0.0;
                NekDouble xi2 = 0.0;
                Array<OneD, const NekDouble> pts;
                int nq0, nq1, nq2;

                // get points;
                //find end points
                for(i = 0; i < m_coordim; ++i)
                {
                    nq0 = m_xmap[i]->GetNumPoints(0);
                    nq1 = m_xmap[i]->GetNumPoints(1);
                    nq2 = m_xmap[i]->GetNumPoints(2);

                    pts = m_xmap[i]->GetPhys();

                    // use projection to side 1 to determine xi_1 coordinate based on length
                    len0 += (pts[nq0-1]-pts[0])*(pts[nq0-1]-pts[0]);
                    xi0  += (coords[i] -pts[0])*(pts[nq0-1]-pts[0]);

                    // use projection to side 4 to determine xi_2 coordinate based on length
                    len1 += (pts[nq0*(nq1-1)]-pts[0])*(pts[nq0*(nq1-1)]-pts[0]);
                    xi1  += (coords[i] -pts[0])*(pts[nq0*(nq1-1)]-pts[0]);

                    // use projection to side 4 to determine xi_3 coordinate based on length
                    len2 += (pts[nq0*nq1*(nq2-1)]-pts[0])*(pts[nq0*nq1*(nq2-1)]-pts[0]);
                    xi2  += (coords[i] -pts[0])*(pts[nq0*nq1*(nq2-1)]-pts[0]);
                }

                Lcoords[0] =  2*xi0/len0-1.0;
                Lcoords[1] =  2*xi1/len1-1.0;
                Lcoords[2] =  2*xi2/len2-1.0;
            }
            else
            {
                // Determine nearest point of coords  to values in m_xmap
                Array<OneD, NekDouble> ptsx = m_xmap[0]->GetPhys();
                Array<OneD, NekDouble> ptsy = m_xmap[1]->GetPhys();
                Array<OneD, NekDouble> ptsz = m_xmap[2]->GetPhys();
                int npts = ptsx.num_elements();
                Array<OneD, NekDouble> tmp1(npts), tmp2(npts);
                const Array<OneD, const NekDouble> za = m_xmap[0]->GetPoints(0);
                const Array<OneD, const NekDouble> zb = m_xmap[0]->GetPoints(1);
                const Array<OneD, const NekDouble> zc = m_xmap[0]->GetPoints(2);
                
                //guess the first local coords based on nearest point
                Vmath::Sadd(npts, -coords[0], ptsx,1,tmp1,1);
                Vmath::Vmul (npts, tmp1,1,tmp1,1,tmp1,1);
                Vmath::Sadd(npts, -coords[1], ptsy,1,tmp2,1);
                Vmath::Vvtvp(npts, tmp2,1,tmp2,1,tmp1,1,tmp1,1);
                Vmath::Sadd(npts, -coords[2], ptsz,1,tmp2,1);
                Vmath::Vvtvp(npts, tmp2,1,tmp2,1,tmp1,1,tmp1,1);
                          
                int min_i = Vmath::Imin(npts,tmp1,1);
                
                // Get Local coordinates
                int qa = za.num_elements(), qb = zb.num_elements();
                Lcoords[2] = zc[min_i/(qa*qb)];
                min_i = min_i%(qa*qb);
                Lcoords[1] = zb[min_i/qa];
                Lcoords[0] = za[min_i%qa];

                // Perform newton iteration to find local coordinates 
                NewtonIterationForLocCoord(coords,Lcoords);
            }
        }

        bool HexGeom::v_ContainsPoint(
            const Array<OneD, const NekDouble> &gloCoord, NekDouble tol)
        {
            ASSERTL1(gloCoord.num_elements() == 3,
                     "Three dimensional geometry expects three coordinates.");

            Array<OneD,NekDouble> stdCoord(GetCoordim(),0.0);
            v_GetLocCoords(gloCoord, stdCoord);
            if (stdCoord[0] >= -(1+tol) && stdCoord[0] <= 1+tol
                && stdCoord[1] >= -(1+tol) && stdCoord[1] <= 1+tol
                && stdCoord[2] >= -(1+tol) && stdCoord[2] <= 1+tol)
            {
                return true;
            }
            return false;
        }

        int HexGeom::v_GetNumVerts() const
        {
            return 8;
        }
        
        int HexGeom::v_GetNumEdges() const
        {
            return 12;
        }

        int HexGeom::v_GetNumFaces() const
        {
            return 6;
        }

        int HexGeom::v_GetVertexEdgeMap(const int i, const int j) const
	{
	    const unsigned int VertexEdgeConnectivity[][3] = {
	        {0,3,4},{0,1,5},{1,2,6},{2,3,7},
                {4,8,11},{5,8,9},{6,9,10},{7,10,11}};

	    return VertexEdgeConnectivity[i][j];
	}

        int HexGeom::v_GetVertexFaceMap(const int i, const int j) const
	{
	    const unsigned int VertexFaceConnectivity[][3] = {
	        {0,1,4},{0,1,2},{0,2,3},{0,3,4},
	        {1,4,5},{1,2,5},{2,3,5},{3,4,5}};

	    return VertexFaceConnectivity[i][j];
	}

        int HexGeom::v_GetEdgeFaceMap(const int i, const int j) const
	{
	    const unsigned int EdgeFaceConnectivity[][2] = {
                {0,1},{0,2},{0,3},{0,4},{1,4},{1,2},{2,3},{3,4},
	        {1,5},{2,5},{3,5},{4,5}};

	    return EdgeFaceConnectivity[i][j];
	}

        void HexGeom::SetUpLocalEdges()
        {
            // find edge 0
            int i,j;
            unsigned int check;

            SegGeomSharedPtr edge;

            // First set up the 4 bottom edges
            int f;
            for(f = 1; f < 5 ; f++)
            {
                check = 0;
                for(i = 0; i < 4; i++)
                {
                    for(j = 0; j < 4; j++)
                    {
                        if( (m_faces[0])->GetEid(i) == (m_faces[f])->GetEid(j) )
                        {
                            edge = boost::dynamic_pointer_cast<SegGeom>((m_faces[0])->GetEdge(i));
                            m_edges.push_back(edge);
                            check++;
                        }
                    }
                }

                if( check < 1 )
                {
                    std::ostringstream errstrm;
                    errstrm << "Connected faces do not share an edge. Faces ";
                    errstrm << (m_faces[0])->GetFid() << ", " << (m_faces[f])->GetFid();
                    ASSERTL0(false, errstrm.str());
                }
                else if( check > 1)
                {
                    std::ostringstream errstrm;
                    errstrm << "Connected faces share more than one edge. Faces ";
                    errstrm << (m_faces[0])->GetFid() << ", " << (m_faces[f])->GetFid();
                    ASSERTL0(false, errstrm.str());
                }
            }

            // Then, set up the 4 vertical edges
            check = 0;
            for(i = 0; i < 4; i++)
            {
                for(j = 0; j < 4; j++)
                {
                    if( (m_faces[1])->GetEid(i) == (m_faces[4])->GetEid(j) )
                    {
                        edge = boost::dynamic_pointer_cast<SegGeom>((m_faces[1])->GetEdge(i));
                        m_edges.push_back(edge);
                        check++;
                    }
                }
            }
            if( check < 1 )
            {
                std::ostringstream errstrm;
                errstrm << "Connected faces do not share an edge. Faces ";
                errstrm << (m_faces[1])->GetFid() << ", " << (m_faces[4])->GetFid();
                ASSERTL0(false, errstrm.str());
            }
            else if( check > 1)
            {
                std::ostringstream errstrm;
                errstrm << "Connected faces share more than one edge. Faces ";
                errstrm << (m_faces[1])->GetFid() << ", " << (m_faces[4])->GetFid();
                ASSERTL0(false, errstrm.str());
            }
            for(f = 1; f < 4 ; f++)
            {
                check = 0;
                for(i = 0; i < 4; i++)
                {
                    for(j = 0; j < 4; j++)
                    {
                        if( (m_faces[f])->GetEid(i) == (m_faces[f+1])->GetEid(j) )
                        {
                            edge = boost::dynamic_pointer_cast<SegGeom>((m_faces[f])->GetEdge(i));
                            m_edges.push_back(edge);
                            check++;
                        }
                    }
                }

                if( check < 1 )
                {
                    std::ostringstream errstrm;
                    errstrm << "Connected faces do not share an edge. Faces ";
                    errstrm << (m_faces[f])->GetFid() << ", " << (m_faces[f+1])->GetFid();
                    ASSERTL0(false, errstrm.str());
                }
                else if( check > 1)
                {
                    std::ostringstream errstrm;
                    errstrm << "Connected faces share more than one edge. Faces ";
                    errstrm << (m_faces[f])->GetFid() << ", " << (m_faces[f+1])->GetFid();
                    ASSERTL0(false, errstrm.str());
                }
            }

            // Finally, set up the 4 top vertices
            for(f = 1; f < 5 ; f++)
            {
                check = 0;
                for(i = 0; i < 4; i++)
                {
                    for(j = 0; j < 4; j++)
                    {
                        if( (m_faces[5])->GetEid(i) == (m_faces[f])->GetEid(j) )
                        {
                            edge = boost::dynamic_pointer_cast<SegGeom>((m_faces[5])->GetEdge(i));
                            m_edges.push_back(edge);
                            check++;
                        }
                    }
                }

                if( check < 1 )
                {
                    std::ostringstream errstrm;
                    errstrm << "Connected faces do not share an edge. Faces ";
                    errstrm << (m_faces[5])->GetFid() << ", " << (m_faces[f])->GetFid();
                    ASSERTL0(false, errstrm.str());
                }
                else if( check > 1)
                {
                    std::ostringstream errstrm;
                    errstrm << "Connected faces share more than one edge. Faces ";
                    errstrm << (m_faces[5])->GetFid() << ", " << (m_faces[f])->GetFid();
                    ASSERTL0(false, errstrm.str());
                }
            }
        }

        void HexGeom::SetUpLocalVertices()
        {
            // Set up the first 2 vertices (i.e. vertex 0,1)
            if( ( m_edges[0]->GetVid(0) == m_edges[1]->GetVid(0) ) ||
                ( m_edges[0]->GetVid(0) == m_edges[1]->GetVid(1) ) )
            {
                m_verts.push_back(m_edges[0]->GetVertex(1));
                m_verts.push_back(m_edges[0]->GetVertex(0));
            }
            else if( ( m_edges[0]->GetVid(1) == m_edges[1]->GetVid(0) ) ||
                     ( m_edges[0]->GetVid(1) == m_edges[1]->GetVid(1) ) )
            {
                m_verts.push_back(m_edges[0]->GetVertex(0));
                m_verts.push_back(m_edges[0]->GetVertex(1));
            }
            else
            {
                std::ostringstream errstrm;
                errstrm << "Connected edges do not share a vertex. Edges ";
                errstrm << m_edges[0]->GetEid() << ", " << m_edges[1]->GetEid();
                ASSERTL0(false, errstrm.str());
            }

            // set up the other bottom vertices (i.e. vertex 2,3)
            int i;
            for(i = 1; i < 3; i++)
            {
                if( m_edges[i]->GetVid(0) == m_verts[i]->GetVid() )
                {
                    m_verts.push_back(m_edges[i]->GetVertex(1));
                }
                else if( m_edges[i]->GetVid(1) == m_verts[i]->GetVid() )
                {
                    m_verts.push_back(m_edges[i]->GetVertex(0));
                }
                else
                {
                    std::ostringstream errstrm;
                    errstrm << "Connected edges do not share a vertex. Edges ";
                    errstrm << m_edges[i]->GetEid() << ", " << m_edges[i-1]->GetEid();
                    ASSERTL0(false, errstrm.str());
                }
            }

            // set up top vertices
            // First, set up vertices 4,5
            if( ( m_edges[8]->GetVid(0) == m_edges[9]->GetVid(0) ) ||
                ( m_edges[8]->GetVid(0) == m_edges[9]->GetVid(1) ) )
            {
                m_verts.push_back(m_edges[8]->GetVertex(1));
                m_verts.push_back(m_edges[8]->GetVertex(0));
            }
            else if( ( m_edges[8]->GetVid(1) == m_edges[9]->GetVid(0) ) ||
                     ( m_edges[8]->GetVid(1) == m_edges[9]->GetVid(1) ) )
            {
                m_verts.push_back(m_edges[8]->GetVertex(0));
                m_verts.push_back(m_edges[8]->GetVertex(1));
            }
            else
            {
                std::ostringstream errstrm;
                errstrm << "Connected edges do not share a vertex. Edges ";
                errstrm << m_edges[8]->GetEid() << ", " << m_edges[9]->GetEid();
                ASSERTL0(false, errstrm.str());
            }

            // set up the other top vertices (i.e. vertex 6,7)
            for(i = 9; i < 11; i++)
            {
                if( m_edges[i]->GetVid(0) == m_verts[i-4]->GetVid() )
                {
                    m_verts.push_back(m_edges[i]->GetVertex(1));
                }
                else if( m_edges[i]->GetVid(1) == m_verts[i-4]->GetVid() )
                {
                    m_verts.push_back(m_edges[i]->GetVertex(0));
                }
                else
                {
                    std::ostringstream errstrm;
                    errstrm << "Connected edges do not share a vertex. Edges ";
                    errstrm << m_edges[i]->GetEid() << ", " << m_edges[i-1]->GetEid();
                    ASSERTL0(false, errstrm.str());
                }
            }
        }

        void HexGeom::SetUpFaceOrientation()
        {
            int f,i;

            // These arrays represent the vector of the A and B
            // coordinate of the local elemental coordinate system
            // where A corresponds with the coordinate direction xi_i
            // with the lowest index i (for that particular face)
            // Coordinate 'B' then corresponds to the other local
            // coordinate (i.e. with the highest index)
            Array<OneD,NekDouble> elementAaxis(m_coordim);
            Array<OneD,NekDouble> elementBaxis(m_coordim);

            // These arrays correspond to the local coordinate
            // system of the face itself (i.e. the Geometry2D)
            // faceAaxis correspond to the xi_0 axis
            // faceBaxis correspond to the xi_1 axis
            Array<OneD,NekDouble> faceAaxis(m_coordim);
            Array<OneD,NekDouble> faceBaxis(m_coordim);

            // This is the base vertex of the face (i.e. the Geometry2D)
            // This corresponds to thevertex with local ID 0 of the
            // Geometry2D
            unsigned int baseVertex;

            // The lenght of the vectors above
            NekDouble elementAaxis_length;
            NekDouble elementBaxis_length;
            NekDouble faceAaxis_length;
            NekDouble faceBaxis_length;

            // This 2D array holds the local id's of all the vertices
            // for every face. For every face, they are ordered in such
            // a way that the implementation below allows a unified approach
            // for all faces.
            const unsigned int faceVerts[kNfaces][QuadGeom::kNverts] =
                { {0,1,2,3} ,
                  {0,1,5,4} ,
                  {1,2,6,5} ,
                  {3,2,6,7} ,
                  {0,3,7,4} ,
                  {4,5,6,7} };

            NekDouble dotproduct1 = 0.0;
            NekDouble dotproduct2 = 0.0;

            unsigned int orientation;

            // Loop over all the faces to set up the orientation
            for(f = 0; f < kNqfaces + kNtfaces; f++)
            {
                // initialisation
                elementAaxis_length = 0.0;
                elementBaxis_length = 0.0;
                faceAaxis_length = 0.0;
                faceBaxis_length = 0.0;

                dotproduct1 = 0.0;
                dotproduct2 = 0.0;

                baseVertex = m_faces[f]->GetVid(0);

                // We are going to construct the vectors representing the A and B axis
                // of every face. These vectors will be constructed as a vector-representation
                // of the edges of the face. However, for both coordinate directions, we can
                // represent the vectors by two different edges. That's why we need to make sure that
                // we pick the edge to which the baseVertex of the Geometry2D-representation of the face
                // belongs...
                if( baseVertex == m_verts[ faceVerts[f][0] ]->GetVid() )
                {
                    for(i = 0; i < m_coordim; i++)
                    {
                        elementAaxis[i] = (*m_verts[ faceVerts[f][1] ])[i] - (*m_verts[ faceVerts[f][0] ])[i];
                        elementBaxis[i] = (*m_verts[ faceVerts[f][3] ])[i] - (*m_verts[ faceVerts[f][0] ])[i];
                    }
                }
                else if( baseVertex == m_verts[ faceVerts[f][1] ]->GetVid() )
                {
                    for(i = 0; i < m_coordim; i++)
                    {
                        elementAaxis[i] = (*m_verts[ faceVerts[f][1] ])[i] - (*m_verts[ faceVerts[f][0] ])[i];
                        elementBaxis[i] = (*m_verts[ faceVerts[f][2] ])[i] - (*m_verts[ faceVerts[f][1] ])[i];
                    }
                }
                else if( baseVertex == m_verts[ faceVerts[f][2] ]->GetVid() )
                {
                    for(i = 0; i < m_coordim; i++)
                    {
                        elementAaxis[i] = (*m_verts[ faceVerts[f][2] ])[i] - (*m_verts[ faceVerts[f][3] ])[i];
                        elementBaxis[i] = (*m_verts[ faceVerts[f][2] ])[i] - (*m_verts[ faceVerts[f][1] ])[i];
                    }
                }
                else if( baseVertex == m_verts[ faceVerts[f][3] ]->GetVid() )
                {
                    for(i = 0; i < m_coordim; i++)
                    {
                        elementAaxis[i] = (*m_verts[ faceVerts[f][2] ])[i] - (*m_verts[ faceVerts[f][3] ])[i];
                        elementBaxis[i] = (*m_verts[ faceVerts[f][3] ])[i] - (*m_verts[ faceVerts[f][0] ])[i];
                    }
                }
                else
                {
                    ASSERTL0(false, "Could not find matching vertex for the face");
                }

                // Now, construct the edge-vectors of the local coordinates of
                // the Geometry2D-representation of the face
                for(i = 0; i < m_coordim; i++)
                {
                    faceAaxis[i] = (*m_faces[f]->GetVertex(1))[i] - (*m_faces[f]->GetVertex(0))[i];
                    faceBaxis[i] = (*m_faces[f]->GetVertex(3))[i] - (*m_faces[f]->GetVertex(0))[i];

                    elementAaxis_length += pow(elementAaxis[i],2);
                    elementBaxis_length += pow(elementBaxis[i],2);
                    faceAaxis_length += pow(faceAaxis[i],2);
                    faceBaxis_length += pow(faceBaxis[i],2);
                }

                elementAaxis_length = sqrt(elementAaxis_length);
                elementBaxis_length = sqrt(elementBaxis_length);
                faceAaxis_length = sqrt(faceAaxis_length);
                faceBaxis_length = sqrt(faceBaxis_length);

                // Calculate the inner product of both the A-axis
                // (i.e. Elemental A axis and face A axis)
                for(i = 0 ; i < m_coordim; i++)
                {
                    dotproduct1 += elementAaxis[i]*faceAaxis[i];
                }

                orientation = 0;
                // if the innerproduct is equal to the (absolute value of the ) products of the lengths
                // of both vectors, then, the coordinate systems will NOT be transposed
                if( fabs(elementAaxis_length*faceAaxis_length - fabs(dotproduct1)) < NekConstants::kNekZeroTol )
                {
                    // if the inner product is negative, both A-axis point
                    // in reverse direction
                    if(dotproduct1 < 0.0)
                    {
                        orientation += 2;
                    }

                    // calculate the inner product of both B-axis
                    for(i = 0 ; i < m_coordim; i++)
                    {
                        dotproduct2 += elementBaxis[i]*faceBaxis[i];
                    }

                    // check that both these axis are indeed parallel
                    ASSERTL1(fabs(elementBaxis_length*faceBaxis_length - fabs(dotproduct2)) <
                             NekConstants::kNekZeroTol,
                             "These vectors should be parallel");

                    // if the inner product is negative, both B-axis point
                    // in reverse direction
                    if( dotproduct2 < 0.0 )
                    {
                        orientation++;
                    }
                }
                // The coordinate systems are transposed
                else
                {
                    orientation = 4;

                    // Calculate the inner product between the elemental A-axis
                    // and the B-axis of the face (which are now the corresponding axis)
                    dotproduct1 = 0.0;
                    for(i = 0 ; i < m_coordim; i++)
                    {
                        dotproduct1 += elementAaxis[i]*faceBaxis[i];
                    }

                    // check that both these axis are indeed parallel
                    ASSERTL1(fabs(elementAaxis_length*faceBaxis_length - fabs(dotproduct1)) <
                             NekConstants::kNekZeroTol,
                             "These vectors should be parallel");

                    // if the result is negative, both axis point in reverse
                    // directions
                    if(dotproduct1 < 0.0)
                    {
                        orientation += 2;
                    }

                    // Do the same for the other two corresponding axis
                    dotproduct2 = 0.0;
                    for(i = 0 ; i < m_coordim; i++)
                    {
                        dotproduct2 += elementBaxis[i]*faceAaxis[i];
                    }

                    // check that both these axis are indeed parallel
                    ASSERTL1(fabs(elementBaxis_length*faceAaxis_length - fabs(dotproduct2)) <
                             NekConstants::kNekZeroTol,
                             "These vectors should be parallel");

                    if( dotproduct2 < 0.0 )
                    {
                        orientation++;
                    }
                }
				
				orientation = orientation + 5;
                // Fill the m_forient array
                m_forient[f] = (StdRegions::Orientation) orientation;
            }
        }

        void HexGeom::SetUpEdgeOrientation()
        {

            // This 2D array holds the local id's of all the vertices
            // for every edge. For every edge, they are ordered to what we
            // define as being Forwards
            const unsigned int edgeVerts[kNedges][2] =
                { {0,1} ,
                  {1,2} ,
                  {2,3} ,
                  {3,0} ,
                  {0,4} ,
                  {1,5} ,
                  {2,6} ,
                  {3,7} ,
                  {4,5} ,
                  {5,6} ,
                  {6,7} ,
                  {7,4} };

            int i;
            for(i = 0; i < kNedges; i++)
            {
                if( m_edges[i]->GetVid(0) == m_verts[ edgeVerts[i][0] ]->GetVid() )
                {
                    m_eorient[i] = StdRegions::eForwards;
                }
                else if( m_edges[i]->GetVid(0) == m_verts[ edgeVerts[i][1] ]->GetVid() )
                {
                    m_eorient[i] = StdRegions::eBackwards;
                }
                else
                {
                    ASSERTL0(false,"Could not find matching vertex for the edge");
                }
            }
        }

    }; //end of namespace
}; //end of namespace

////////////////////////////////////////////////////////////////////////////////
//
//  File: PrismGeom.cpp
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
//  Description: Prismatic geometry definition.
//
////////////////////////////////////////////////////////////////////////////////


#include <SpatialDomains/PrismGeom.h>
#include <SpatialDomains/Geometry1D.h>
#include <SpatialDomains/Geometry2D.h>
#include <StdRegions/StdPrismExp.h>
#include <SpatialDomains/SegGeom.h>
#include <SpatialDomains/QuadGeom.h>
#include <SpatialDomains/MeshComponents.h>
#include <SpatialDomains/GeomFactors3D.h>


namespace Nektar
{
    namespace SpatialDomains
    {
        
        PrismGeom::PrismGeom()
        {
            m_geomShapeType = ePrism;
        }

        PrismGeom::PrismGeom(const Geometry2DSharedPtr faces[]):
            Geometry3D(faces[0]->GetEdge(0)->GetVertex(0)->GetCoordim())
        {
            m_geomShapeType = ePrism;
            
            /// Copy the face shared pointers.
            m_faces.insert(m_faces.begin(), faces, faces+PrismGeom::kNfaces);

            /// Set up orientation vectors with correct amount of elements.
            m_eorient.resize(kNedges);
            m_forient.resize(kNfaces);
            
            /// Set up local objects.
            SetUpLocalEdges();
            SetUpLocalVertices();
            SetUpEdgeOrientation();
            SetUpFaceOrientation();
            
            /// Determine necessary order for standard region.
            vector<int> tmp;
            
            int order0, points0, order1, points1;
            
            if (m_forient[0] < 9)
            {
                tmp.push_back(faces[0]->GetXmap(0)->GetEdgeNcoeffs(0));
                tmp.push_back(faces[0]->GetXmap(0)->GetEdgeNcoeffs(2));
                order0 = *max_element(tmp.begin(), tmp.end());

                tmp.clear();
                tmp.push_back(faces[0]->GetXmap(0)->GetEdgeNumPoints(0));
                tmp.push_back(faces[0]->GetXmap(0)->GetEdgeNumPoints(2));
                points0 = *max_element(tmp.begin(), tmp.end());
            }
            else
            {
                tmp.push_back(faces[0]->GetXmap(0)->GetEdgeNcoeffs(1));
                tmp.push_back(faces[0]->GetXmap(0)->GetEdgeNcoeffs(3));
                order0 = *max_element(tmp.begin(), tmp.end());

                tmp.clear();
                tmp.push_back(faces[0]->GetXmap(0)->GetEdgeNumPoints(1));
                tmp.push_back(faces[0]->GetXmap(0)->GetEdgeNumPoints(3));
                points0 = *max_element(tmp.begin(), tmp.end());
            }
            
            if (m_forient[0] < 9)
            {
                tmp.clear();
                tmp.push_back(faces[0]->GetXmap(0)->GetEdgeNcoeffs(1));
                tmp.push_back(faces[0]->GetXmap(0)->GetEdgeNcoeffs(3));
                tmp.push_back(faces[2]->GetXmap(0)->GetEdgeNcoeffs(2));
                order1 = *max_element(tmp.begin(), tmp.end());
                
                tmp.clear();
                tmp.push_back(faces[0]->GetXmap(0)->GetEdgeNumPoints(1));
                tmp.push_back(faces[0]->GetXmap(0)->GetEdgeNumPoints(3));
                tmp.push_back(faces[2]->GetXmap(0)->GetEdgeNumPoints(2));
                points1 = *max_element(tmp.begin(), tmp.end());
            }
            else
            {
                tmp.clear();
                tmp.push_back(faces[0]->GetXmap(0)->GetEdgeNcoeffs(0));
                tmp.push_back(faces[0]->GetXmap(0)->GetEdgeNcoeffs(2));
                tmp.push_back(faces[2]->GetXmap(0)->GetEdgeNcoeffs(2));
                order1 = *max_element(tmp.begin(), tmp.end());
                
                tmp.clear();
                tmp.push_back(faces[0]->GetXmap(0)->GetEdgeNumPoints(0));
                tmp.push_back(faces[0]->GetXmap(0)->GetEdgeNumPoints(2));
                tmp.push_back(faces[2]->GetXmap(0)->GetEdgeNumPoints(2));
                points1 = *max_element(tmp.begin(), tmp.end());
            }
            
            tmp.clear();
            tmp.push_back(order0);
            tmp.push_back(order1);
            tmp.push_back(faces[1]->GetXmap(0)->GetEdgeNcoeffs(1));
            tmp.push_back(faces[1]->GetXmap(0)->GetEdgeNcoeffs(2));
            tmp.push_back(faces[3]->GetXmap(0)->GetEdgeNcoeffs(1));
            tmp.push_back(faces[3]->GetXmap(0)->GetEdgeNcoeffs(2));
            int order2 = *max_element(tmp.begin(), tmp.end());
            
            tmp.clear();
            tmp.push_back(points0);
            tmp.push_back(points1);
            tmp.push_back(faces[1]->GetXmap(0)->GetEdgeNumPoints(1));
            tmp.push_back(faces[1]->GetXmap(0)->GetEdgeNumPoints(2));
            tmp.push_back(faces[3]->GetXmap(0)->GetEdgeNumPoints(1));
            tmp.push_back(faces[3]->GetXmap(0)->GetEdgeNumPoints(2));
            tmp.push_back(faces[1]->GetEdge(1)->GetBasis(0,0)->GetNumPoints());
            tmp.push_back(faces[1]->GetEdge(2)->GetBasis(0,0)->GetNumPoints());
            tmp.push_back(faces[3]->GetEdge(1)->GetBasis(0,0)->GetNumPoints());
            tmp.push_back(faces[3]->GetEdge(2)->GetBasis(0,0)->GetNumPoints());
            int points2 = *max_element(tmp.begin(), tmp.end());
            
            const LibUtilities::BasisKey A(
                LibUtilities::eModified_A, order0,
                LibUtilities::PointsKey(points0,LibUtilities::eGaussLobattoLegendre));
            const LibUtilities::BasisKey B(
                LibUtilities::eModified_A, order1,
                LibUtilities::PointsKey(points1,LibUtilities::eGaussLobattoLegendre));
            const LibUtilities::BasisKey C(
                LibUtilities::eModified_B, order2,
                LibUtilities::PointsKey(points2,LibUtilities::eGaussRadauMAlpha1Beta0));

            m_xmap = Array<OneD, StdRegions::StdExpansion3DSharedPtr>(m_coordim);

            for(int i = 0; i < m_coordim; ++i)
            {
                m_xmap[i] = MemoryManager<StdRegions::StdPrismExp>::AllocateSharedPtr(A,B,C);
            }
        }

        PrismGeom::~PrismGeom()
        {
        }
        
        int PrismGeom::v_GetNumVerts() const
        {
            return 6;
        }
        
        int PrismGeom::v_GetNumEdges() const
        {
            return 9;
        }

        
        /**
         * @brief Determines if a point specified in global coordinates is
         * located within this tetrahedral geometry.
         */
        bool PrismGeom::v_ContainsPoint(
            const Array<OneD, const NekDouble> &gloCoord, NekDouble tol)
        {
            // Validation checks
            ASSERTL1(gloCoord.num_elements() == 3,
                     "Three dimensional geometry expects three coordinates.");
            
            // Convert to the local (eta) coordinates.
            Array<OneD,NekDouble> locCoord(GetCoordim(),0.0);
            v_GetLocCoords(gloCoord, locCoord);
            
            // Check local coordinate is within [-1,1]^3 bounds.
            if (locCoord[0] >= -(1+tol) && locCoord[1] >= -(1+tol) &&
                locCoord[2] >= -(1+tol) && locCoord[1] <=  (1+tol) &&
                locCoord[0] + locCoord[2] <= tol)
            {
                return true;
            }
            
            return false;
        }

        void PrismGeom::v_GenGeomFactors(
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

                // check to see if all quadrilateral faces are parallelograms
                if(Gtype == eRegular)
                {
                    // Vertex ids of quad faces
                    const unsigned int faceVerts[3][4] =
                        { {0,1,2,3} ,
                          {1,2,5,4} ,
                          {0,3,5,4} };

                    for(f = 0; f < 3; f++)
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
                    Gtype, m_coordim, m_xmap, tbasis);

                m_geomFactorsState = ePtsFilled;
            }
        }


        void PrismGeom::v_GetLocCoords(
            const Array<OneD, const NekDouble> &coords, 
                  Array<OneD,       NekDouble> &Lcoords)
        {
            // calculate local coordinate for coord
            if(GetGtype() == eRegular)
            {
                // Point inside tetrahedron
                VertexComponent r(m_coordim, 0, coords[0], coords[1], coords[2]);

                // Edges
                VertexComponent er0, e10, e30, e40;
                er0.Sub(r,*m_verts[0]);
                e10.Sub(*m_verts[1],*m_verts[0]);
                e30.Sub(*m_verts[3],*m_verts[0]);
                e40.Sub(*m_verts[4],*m_verts[0]);


                // Cross products (Normal times area)
                VertexComponent cp1030, cp3040, cp4010;
                cp1030.Mult(e10,e30);
                cp3040.Mult(e30,e40);
                cp4010.Mult(e40,e10);


                // Barycentric coordinates (relative volume)
                NekDouble V = e40.dot(cp1030); // Prism Volume = {(e40)dot(e10)x(e30)}/2
                NekDouble beta  = er0.dot(cp3040) / (2.0*V); // volume1 = {(er0)dot(e30)x(e40)}/4
                NekDouble gamma = er0.dot(cp4010) / (3.0*V); // volume2 = {(er0)dot(e40)x(e10)}/6
                NekDouble delta = er0.dot(cp1030) / (2.0*V); // volume3 = {(er0)dot(e10)x(e30)}/4

                // Make Prism bigger
                Lcoords[0] = 2.0*beta  - 1.0;
                Lcoords[1] = 2.0*gamma - 1.0;
                Lcoords[2] = 2.0*delta - 1.0;

            }
            else
            {
                v_FillGeom();
            
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
                
                // Get collapsed coordinate
                int qa = za.num_elements(), qb = zb.num_elements();
                Lcoords[2] = zc[min_i/(qa*qb)];
                min_i = min_i%(qa*qb);
                Lcoords[1] = zb[min_i/qa];
                Lcoords[0] = za[min_i%qa];

                // recover cartesian coordinate from collapsed coordinate. 
                Lcoords[0] = (1.0+Lcoords[0])*(1.0-Lcoords[2])/2 -1.0;            
                Lcoords[1] = (1.0+Lcoords[0])*(1.0-Lcoords[2])/2 -1.0;


                // Perform newton iteration to find local coordinates 
                NewtonIterationForLocCoord(coords,Lcoords);
            }
        }

        void PrismGeom::SetUpLocalEdges(){
            // find edge 0
            int i,j;
            unsigned int check;

            SegGeomSharedPtr edge;

            // First set up the 4 bottom edges
            int f;  //  Connected face index
            for (f = 1; f < 5 ; f++)
            {
                int nEdges = m_faces[f]->GetNumEdges();
                check = 0;
                for (i = 0; i < 4; i++)
                {
                    for (j = 0; j < nEdges; j++)
                    {
                        if (m_faces[0]->GetEid(i) == m_faces[f]->GetEid(j))
                        {
                            edge = boost::dynamic_pointer_cast<SegGeom>((m_faces[0])->GetEdge(i));
                            m_edges.push_back(edge);
                            check++;
                        }
                    }
                }

                if (check < 1)
                {
                    std::ostringstream errstrm;
                    errstrm << "Connected faces do not share an edge. Faces ";
                    errstrm << (m_faces[0])->GetFid() << ", " << (m_faces[f])->GetFid();
                    ASSERTL0(false, errstrm.str());
                }
                else if (check > 1)
                {
                    std::ostringstream errstrm;
                    errstrm << "Connected faces share more than one edge. Faces ";
                    errstrm << (m_faces[0])->GetFid() << ", " << (m_faces[f])->GetFid();
                    ASSERTL0(false, errstrm.str());
                }
            }

            // Then, set up the 4 vertical edges
            check = 0;
            for(i = 0; i < 3; i++) //Set up the vertical edge :face(1) and face(4)
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
            // Set up vertical edges: face(1) through face(4)
            for(f = 1; f < 4 ; f++)
            {
                check = 0;
                for(i = 0; i < m_faces[f]->GetNumEdges(); i++)
                {
                    for(j = 0; j < m_faces[f+1]->GetNumEdges(); j++)
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

            // Finally, set up the 1 top edge
            check = 0;
            for(i = 0; i < 4; i++)
            {
                for(j = 0; j < 4; j++)
                {
                    if( (m_faces[2])->GetEid(i) == (m_faces[4])->GetEid(j) )
                    {
                        edge = boost::dynamic_pointer_cast<SegGeom>((m_faces[2])->GetEdge(i));
                        m_edges.push_back(edge);
                        check++;
                    }
                }
            }

            if( check < 1 )
            {
                std::ostringstream errstrm;
                errstrm << "Connected faces do not share an edge. Faces ";
                errstrm << (m_faces[1])->GetFid() << ", " << (m_faces[3])->GetFid();
                ASSERTL0(false, errstrm.str());
            }
            else if( check > 1)
            {
                std::ostringstream errstrm;
                errstrm << "Connected faces share more than one edge. Faces ";
                errstrm << (m_faces[1])->GetFid() << ", " << (m_faces[3])->GetFid();
                ASSERTL0(false, errstrm.str());
            }
        };


        void PrismGeom::SetUpLocalVertices(){

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
            for(int i = 1; i < 3; i++)
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
            if( (m_edges[8]->GetVid(0) == m_edges[4]->GetVid(0)) ||
                (m_edges[8]->GetVid(0) == m_edges[4]->GetVid(1))  )
            {
                m_verts.push_back(m_edges[8]->GetVertex(0));
                m_verts.push_back(m_edges[8]->GetVertex(1));
            }
            else if( (m_edges[8]->GetVid(1) == m_edges[4]->GetVid(0)) ||
                     (m_edges[8]->GetVid(1) == m_edges[4]->GetVid(1))  )
            {
                m_verts.push_back(m_edges[8]->GetVertex(1));
                m_verts.push_back(m_edges[8]->GetVertex(0));
            }
            else
            {
                std::ostringstream errstrm;
                errstrm << "Connected edges do not share a vertex. Edges ";
                errstrm << m_edges[8]->GetEid();
                ASSERTL0(false, errstrm.str());
            }
        };

        void PrismGeom::SetUpEdgeOrientation(){

            // This 2D array holds the local id's of all the vertices
            // for every edge. For every edge, they are ordered to what we
            // define as being Forwards
            const unsigned int edgeVerts[kNedges][2] =
                { {0,1} ,
                  {1,2} ,
                  {3,2} ,
                  {0,3} ,
                  {0,4} ,
                  {1,4} ,
                  {2,5} ,
                  {3,5} ,
                  {4,5} };


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


        };

        void PrismGeom::SetUpFaceOrientation(){
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
                  {0,1,4,0},  // This is triangle requires only three vertices
                  {1,2,5,4} ,
                  {3,2,5,0},  // This is triangle requires only three vertices
                  {0,3,5,4} ,};

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

                // Compute the length of edges on a base-face
                if( f==1  ||  f==3 ) {   // Face is a Triangle
                    for(i = 0; i < m_coordim; i++)
                    {
                        elementAaxis[i] = (*m_verts[ faceVerts[f][1] ])[i] - (*m_verts[ faceVerts[f][0] ])[i];
                        elementBaxis[i] = (*m_verts[ faceVerts[f][2] ])[i] - (*m_verts[ faceVerts[f][0] ])[i];
                    }
                }
                else { // Face is a Quad
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
                }
                // Now, construct the edge-vectors of the local coordinates of
                // the Geometry2D-representation of the face
                for(i = 0; i < m_coordim; i++)
                {
                    int v = m_faces[f]->GetNumVerts()-1;
                    faceAaxis[i] = (*m_faces[f]->GetVertex(1))[i] - (*m_faces[f]->GetVertex(0))[i];
                    faceBaxis[i] = (*m_faces[f]->GetVertex(v))[i] - (*m_faces[f]->GetVertex(0))[i];

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

//                     // check that both these axis are indeed parallel
//                     ASSERTL1(fabs(elementBaxis_length*faceBaxis_length - fabs(dotproduct2)) <
//                              StdRegions::NekConstants::kEvaluateTol,
//                              "These vectors should be parallel");

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
                    if (fabs(elementAaxis_length*faceBaxis_length
                            - fabs(dotproduct1)) > NekConstants::kNekZeroTol)
                    {
                        cout << "Warning: Prism axes not parallel" << endl;
                    }

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
                    if (fabs(elementBaxis_length*faceAaxis_length
                            - fabs(dotproduct2)) > NekConstants::kNekZeroTol)
                    {
                        cout << "Warning: Prism axes not parallel" << endl;
                    }

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
    }; //end of namespace
}; //end of namespace

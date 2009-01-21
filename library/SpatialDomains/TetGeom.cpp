////////////////////////////////////////////////////////////////////////////////
//
//  File:  $Source: /usr/sci/projects/Nektar/cvs/Nektar++/library/SpatialDomains/TetGeom.cpp,v $
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

#include <SpatialDomains/TetGeom.h>



namespace Nektar
{
    namespace SpatialDomains
    {

        TetGeom::TetGeom()
        {
            m_GeomShapeType = eTetrahedron;
        }

        TetGeom::TetGeom(const TriGeomSharedPtr faces[]):
                 Geometry3D(faces[0]->GetEdge(0)->GetVertex(0)->GetCoordim())
        {
            m_GeomShapeType = eTetrahedron;
            
            /// Copy the face shared pointers
            m_faces.insert(m_faces.begin(), faces, faces+TetGeom::kNfaces);

            SetUpLocalEdges();
            SetUpLocalVertices();
            SetUpEdgeOrientation();
            SetUpFaceOrientation();

            // BasisKey (const BasisType btype, const int nummodes, const PointsKey pkey)
            //PointsKey (const int &numpoints, const PointsType &pointstype)
            const LibUtilities::BasisKey A(LibUtilities::eModified_A, 2,
                                           LibUtilities::PointsKey(3,LibUtilities::eGaussLobattoLegendre));
            const LibUtilities::BasisKey B(LibUtilities::eModified_B, 2,
                                           LibUtilities::PointsKey(3,LibUtilities::eGaussRadauMAlpha1Beta0));
            const LibUtilities::BasisKey C(LibUtilities::eModified_C, 2,
                                           LibUtilities::PointsKey(3,LibUtilities::eGaussRadauMAlpha2Beta0));
            
            m_xmap = Array<OneD, StdRegions::StdExpansion3DSharedPtr>(m_coordim);

            for(int i = 0; i < m_coordim; ++i)
            {
                m_xmap[i] = MemoryManager<StdRegions::StdTetExp>::AllocateSharedPtr(A,B,C);
            }
        }

        TetGeom::TetGeom (const TriGeomSharedPtr faces[], const StdRegions::FaceOrientation forient[])
        {
            m_GeomShapeType = eTetrahedron;

            /// Copy the face shared pointers
            m_faces.insert(m_faces.begin(), faces, faces+TetGeom::kNfaces);

            for (int j=0; j<kNfaces; ++j)
            {
               m_forient[j] = forient[j];
            }

            m_coordim = faces[0]->GetEdge(0)->GetVertex(0)->GetCoordim();
            ASSERTL0(m_coordim > 2,"Cannot call function with dim == 2");
        }
        
         TetGeom::TetGeom(const VertexComponentSharedPtr verts[], const SegGeomSharedPtr edges[],
                          const TriGeomSharedPtr faces[], const StdRegions::EdgeOrientation eorient[],
                          const StdRegions::FaceOrientation forient[])
         {
			m_GeomShapeType = eTetrahedron;
 
            /// Copy the vert shared pointers.
            m_verts.insert(m_verts.begin(), verts, verts+TetGeom::kNverts);

            /// Copy the edge shared pointers.
            m_edges.insert(m_edges.begin(), edges, edges+TetGeom::kNedges);

            /// Copy the face shared pointers
            m_faces.insert(m_faces.begin(), faces, faces+TetGeom::kNfaces);

            for (int i=0; i<kNedges; ++i)
            {
                m_eorient[i] = eorient[i];
            }

            for (int j=0; j<kNfaces; ++j)
            {
               m_forient[j] = forient[j];
            }

            m_coordim = verts[0]->GetCoordim();
            ASSERTL0(m_coordim > 2,"Cannot call function with dim == 2");
        }
       
        TetGeom::TetGeom(const SegGeomSharedPtr edges[], const StdRegions::EdgeOrientation eorient[])          
        {
            m_GeomShapeType = eTetrahedron;

            /// Copy the edge shared pointers.
            m_edges.insert(m_edges.begin(), edges, edges+TetGeom::kNedges);

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
            ASSERTL0(m_coordim > 2,"Cannot call function with dim == 2");
        }

        TetGeom::~TetGeom()
        {
        }

        void TetGeom::SetUpLocalEdges(){
        
            // find edge 0
            int i,j;
            unsigned int check;

            SegGeomSharedPtr edge;

            // First set up the 3 bottom edges
            int faceConnected;
            for(faceConnected = 1; faceConnected < 4 ; faceConnected++)
            {
                check = 0;
                for(i = 0; i < 3; i++)
                {
                    for(j = 0; j < 3; j++)
                    {
                        if( (m_faces[0])->GetEid(i) == (m_faces[faceConnected])->GetEid(j) )
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
                    errstrm << (m_faces[0])->GetFid() << ", " << (m_faces[faceConnected])->GetFid();
                    ASSERTL0(false, errstrm.str());
                }
                else if( check > 1)
                {
                    std::ostringstream errstrm;
                    errstrm << "Connected faces share more than one egde. Faces ";
                    errstrm << (m_faces[0])->GetFid() << ", " << (m_faces[faceConnected])->GetFid();
                    ASSERTL0(false, errstrm.str());
                }                
            }           
            
            // Then, set up the 3 vertical edges
            check = 0;
            for(i = 0; i < 3; i++) //Set up the vertical edge :face(1) and face(3)
            {
                for(j = 0; j < 3; j++)
                {
                    if( (m_faces[1])->GetEid(i) == (m_faces[3])->GetEid(j) )
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
                errstrm << (m_faces[1])->GetFid() << ", " << (m_faces[3])->GetFid();
                ASSERTL0(false, errstrm.str());
            }
            else if( check > 1)
            {
                std::ostringstream errstrm;
                errstrm << "Connected faces share more than one egde. Faces ";
                errstrm << (m_faces[1])->GetFid() << ", " << (m_faces[3])->GetFid();
                ASSERTL0(false, errstrm.str());
            }
            // Set up vertical edges: face(1) through face(3)
            for(faceConnected = 1; faceConnected < 3 ; faceConnected++)
            {
                check = 0;
                for(i = 0; i < 3; i++)
                {
                    for(j = 0; j < 3; j++)
                    {
                        if( (m_faces[faceConnected])->GetEid(i) == (m_faces[faceConnected+1])->GetEid(j) )
                        {
                            edge = boost::dynamic_pointer_cast<SegGeom>((m_faces[faceConnected])->GetEdge(i));
                            m_edges.push_back(edge);
                            check++;
                        }
                    }
                }
                
                if( check < 1 )
                {
                    std::ostringstream errstrm;
                    errstrm << "Connected faces do not share an edge. Faces ";
                    errstrm << (m_faces[faceConnected])->GetFid() << ", " << (m_faces[faceConnected+1])->GetFid();
                    ASSERTL0(false, errstrm.str());
                }
                else if( check > 1)
                {
                    std::ostringstream errstrm;
                    errstrm << "Connected faces share more than one egde. Faces ";
                    errstrm << (m_faces[faceConnected])->GetFid() << ", " << (m_faces[faceConnected+1])->GetFid();
                    ASSERTL0(false, errstrm.str());
                }                
            }

        };

        
        void TetGeom::SetUpLocalVertices(){
        
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

            // set up the other bottom vertices (i.e. vertex 2)
            for(int i = 1; i < 2; i++)
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
            // First, set up top vertice 3 TODO


            
        };
        void TetGeom::SetUpEdgeOrientation(){
        
            // This 2D array holds the local id's of all the vertices
            // for every edge. For every edge, they are ordered to what we 
            // define as being Forwards
            const unsigned int edgeVerts[kNedges][2] = 
                { {0,1},
                  {1,2},
                  {0,2},
                  {0,3},
                  {1,3},
                  {2,3} };

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
        
        void TetGeom::SetUpFaceOrientation(){
        
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
            const unsigned int faceVerts[kNfaces][TriGeom::kNverts] =  
                { {0,1,2} ,
                  {0,1,3} ,
                  {1,2,3} ,
                  {0,2,3} };

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
                if( baseVertex == m_verts[ faceVerts[f][0] ]->GetVid() )
                {
                    for(i = 0; i < m_coordim; i++)
                    {                    
                        elementAaxis[i] = (*m_verts[ faceVerts[f][1] ])[i] - (*m_verts[ faceVerts[f][0] ])[i];
                        elementBaxis[i] = (*m_verts[ faceVerts[f][2] ])[i] - (*m_verts[ faceVerts[f][0] ])[i];
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
                        elementAaxis[i] = (*m_verts[ faceVerts[f][1] ])[i] - (*m_verts[ faceVerts[f][2] ])[i];
                        elementBaxis[i] = (*m_verts[ faceVerts[f][2] ])[i] - (*m_verts[ faceVerts[f][0] ])[i];
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
                    faceBaxis[i] = (*m_faces[f]->GetVertex(2))[i] - (*m_faces[f]->GetVertex(0))[i];
                
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

                // Fill the m_forient array
                m_forient[f] = (StdRegions::FaceOrientation) orientation;
            }


        };

        void TetGeom::AddElmtConnected(int gvo_id, int locid)
        {
            CompToElmt ee(gvo_id,locid);
            m_elmtmap.push_back(ee);
        }


        int TetGeom::NumElmtConnected() const
        {
            return int(m_elmtmap.size());
        }


        bool TetGeom::IsElmtConnected(int gvo_id, int locid) const
        {
            std::list<CompToElmt>::const_iterator def;
            CompToElmt ee(gvo_id,locid);

            def = find(m_elmtmap.begin(),m_elmtmap.end(),ee);

            // Found the element connectivity object in the list
            return (def != m_elmtmap.end());
        }

        /** given local collapsed coordinate Lcoord return the value of
        physical coordinate in direction i **/

        NekDouble TetGeom::GetCoord(const int i, 
                                          const Array<OneD, const NekDouble> &Lcoord)
        {
            ASSERTL1(m_state == ePtsFilled,
                "Goemetry is not in physical space");

            return m_xmap[i]->PhysEvaluate(Lcoord);
        }

        // Set up GeoFac for this geometry using Coord quadrature distribution

        void TetGeom::GenGeomFactors(const Array<OneD, const LibUtilities::BasisSharedPtr> &tbasis)
        {
            GeomType Gtype = eRegular;
            
            FillGeom();

            // check to see if expansions are linear
            for(int i = 0; i < m_coordim; ++i)
            {
                if((m_xmap[i]->GetBasisNumModes(0) != 2)||
                   (m_xmap[i]->GetBasisNumModes(1) != 2)||
                   (m_xmap[i]->GetBasisNumModes(2) != 2))
                {
                    Gtype = eDeformed;
                }
            }
            m_geomfactors = MemoryManager<GeomFactors>::AllocateSharedPtr(Gtype, m_coordim, m_xmap, tbasis);
        }

          /** \brief put all quadrature information into edge structure 
        and backward transform 

        Note verts, edges, and faces are listed according to anticlockwise
        convention but points in _coeffs have to be in array format from
        left to right.

        */

        void TetGeom::FillGeom()
        {
    // check to see if geometry structure is already filled
            if(m_state != ePtsFilled)
            {
                int i,j,k;
                int nFaceCoeffs = m_xmap[0]->GetFaceNcoeffs(0); 

                Array<OneD, unsigned int> mapArray (nFaceCoeffs);
                Array<OneD, int>    signArray(nFaceCoeffs);
                NekDouble sign;

                for(i = 0; i < kNfaces; i++)
                {
                    m_faces[i]->FillGeom();
                    m_xmap[0]->GetFaceToElementMap(i,m_forient[i],mapArray,signArray); 
                    
                    nFaceCoeffs = m_xmap[0]->GetFaceNcoeffs(i);

                    for(j = 0 ; j < m_coordim; j++)
                    {
                        for(k = 0; k < nFaceCoeffs; k++)
                        {
//                             const Array<OneD, const NekDouble> & coeffs = (*m_faces[i])[j]->GetCoeffs(); //TODO fix this
//                             double v = signArray[k]* coeffs[k];
//                             (m_xmap[j]->UpdateCoeffs())[ mapArray[k] ] = v;

                        }
                    }
                }
                
                for(i = 0; i < m_coordim; ++i)
                {
                    m_xmap[i]->BwdTrans(m_xmap[i]->GetCoeffs(), m_xmap[i]->UpdatePhys());
                }
                
                m_state = ePtsFilled;
            }

        }


       void TetGeom::GetLocCoords(const Array<OneD, const NekDouble>& coords, Array<OneD, NekDouble>& Lcoords)
       {
            FillGeom();

            // calculate local coordinate for coord 
            if(GetGtype() == eRegular)
            {   // Based on Spen's book, page 99

                // Point inside tetrahedron
                VertexComponent r(m_coordim, 0, coords[0], coords[1], coords[2]);

                // Edges
                VertexComponent er0, e10, e20, e30;
                er0.Sub(r,*m_verts[0]);
                e10.Sub(*m_verts[1],*m_verts[0]);
                e20.Sub(*m_verts[2],*m_verts[0]);
                e30.Sub(*m_verts[3],*m_verts[0]);


                // Cross products (Normal times area)
                VertexComponent cp1020, cp2030, cp3010;
                cp1020.Mult(e10,e20);
                cp2030.Mult(e20,e30);
                cp3010.Mult(e30,e10);


                // Barycentric coordinates (relative volume)
                NekDouble V = e30.dot(cp1020); // Tet Volume = {(e30)dot(e10)x(e20)}/6
                NekDouble beta  = er0.dot(cp2030) / V; // volume1 = {(er0)dot(e20)x(e30)}/6
                NekDouble gamma = er0.dot(cp3010) / V; // volume1 = {(er0)dot(e30)x(e10)}/6
                NekDouble delta = er0.dot(cp1020) / V; // volume1 = {(er0)dot(e10)x(e20)}/6

                // Make tet bigger
                Lcoords[0] = 2.0*beta  - 1.0;
                Lcoords[1] = 2.0*gamma - 1.0;
                Lcoords[2] = 2.0*delta - 1.0;
                
            }
            else
            {
                NEKERROR(ErrorUtil::efatal,
                    "inverse mapping must be set up to use this call");
            }

       }



        
    }; //end of namespace
}; //end of namespace

//
// $Log: TetGeom.cpp,v $
// Revision 1.16  2008/12/18 14:08:59  pvos
// NekConstants update
//
// Revision 1.15  2008/11/17 08:58:42  ehan
// Added necessary mapping routines for Tet
//
// Revision 1.14  2008/06/30 19:35:44  ehan
// Fixed infinity recursive-loop error.
//
// Revision 1.13  2008/06/18 19:28:27  ehan
// Added some comments.
//
// Revision 1.12  2008/06/14 01:23:50  ehan
// Implemented constructor and FillGeom().
//
// Revision 1.11  2008/06/12 21:22:55  delisi
// Added method stubs for GenGeomFactors, FillGeom, and GetLocCoords.
//
// Revision 1.10  2008/06/11 16:10:12  delisi
// Added the 3D reader.
//
// Revision 1.9  2008/05/29 19:02:55  delisi
// Renamed eTet to eTetrahedron.
//
// Revision 1.8  2008/05/28 21:52:27  jfrazier
// Added GeomShapeType initialization for the different shapes.
//
// Revision 1.7  2008/05/12 17:30:18  ehan
// Added virtual functions
//
// Revision 1.6  2008/04/06 06:00:38  bnelson
// Changed ConstArray to Array<const>
//
// Revision 1.5  2008/03/25 08:42:16  ehan
// Added constructor, GetLocCoords(), and FillGeom().
//
// Revision 1.4  2008/02/08 23:05:52  jfrazier
// More work on 3D components.
//
// Revision 1.3  2008/02/05 00:43:10  ehan
// Included geometry3D and meshgraphics inorder to prevent compile error.
//
// Revision 1.2  2008/02/03 05:05:16  jfrazier
// Initial checkin of 3D components.
//
// Revision 1.1  2006/05/04 18:59:04  kirby
// *** empty log message ***
//
// Revision 1.10  2006/04/09 02:08:36  jfrazier
// Added precompiled header.
//
// Revision 1.9  2006/03/12 11:06:40  sherwin
//
// First complete copy of code standard code but still not compile tested
//
// Revision 1.8  2006/02/19 01:37:34  jfrazier
// Initial attempt at bringing into conformance with the coding standard.  Still more work to be done.  Has not been compiled.
//
//

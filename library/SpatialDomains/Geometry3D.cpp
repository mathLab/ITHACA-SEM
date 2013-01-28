////////////////////////////////////////////////////////////////////////////////
//
//  File:  Geometry3D.cpp
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
//  Description: 3D geometry information.
//
////////////////////////////////////////////////////////////////////////////////

#include <SpatialDomains/Geometry3D.h>

#include <SpatialDomains/Geometry2D.h>
#include <SpatialDomains/GeomFactors3D.h>
#include <SpatialDomains/MeshComponents.h>
#include <SpatialDomains/SegGeom.h>     // for SegGeomSharedPtr

namespace Nektar
{
  namespace SpatialDomains
  {
      Geometry3D::Geometry3D()
      {
      }

      Geometry3D::Geometry3D(const int coordim):
          Geometry(coordim)
      {
          ASSERTL0(m_coordim > 2,
                   "Coordinate dimension should be at least 3 for a 3D geometry.");
      }

      Geometry3D::~Geometry3D()
      {
      }


      //---------------------------------------
      // Helper functions
      //---------------------------------------

      /**
       * @brief Return the ID of edge i in this element.
       */
      int Geometry3D::GetEid(int i) const
      {
          return v_GetEid(i);
      }

      /**
       * @brief Return face i in this element.
       */
      Geometry2DSharedPtr Geometry3D::GetFace(int i)
      {
          return v_GetFace(i);
      }

      /**
       * @brief Return the orientation of face i in this element.
       */
      StdRegions::Orientation Geometry3D::GetFaceOrient(const int i) const
      {
          return v_GetFaceOrient(i);
      }

      /**
       * @brief Return the ID of face i in this element.
       */
      int Geometry3D::GetFid(int i) const
      {
          return v_GetFid(i);
      }


      //---------------------------------------
      // 3D Geometry Methods
      //---------------------------------------
      void Geometry3D::NewtonIterationForLocCoord
      (const Array<OneD, const NekDouble> &coords, 
                          Array<OneD,NekDouble> &Lcoords)
      {
          
          Array<OneD, NekDouble> ptsx = m_xmap[0]->GetPhys();
          Array<OneD, NekDouble> ptsy = m_xmap[1]->GetPhys();
          Array<OneD, NekDouble> ptsz = m_xmap[2]->GetPhys();
          NekDouble xmap,ymap,zmap, F1,F2, F3;
          NekDouble der1_x, der2_x, der3_x, der1_y, der2_y, der3_y,
              der1_z, der2_z, der3_z; 
          const Array<TwoD, const NekDouble> &gmat = GetGmat();
          
          // Unfortunately need the points in an Array to interpolate
          Array<OneD, NekDouble> D1Dx(ptsx.num_elements(),&gmat[0][0]);
          Array<OneD, NekDouble> D1Dy(ptsx.num_elements(),&gmat[1][0]);
          Array<OneD, NekDouble> D1Dz(ptsx.num_elements(),&gmat[2][0]);
          Array<OneD, NekDouble> D2Dx(ptsx.num_elements(),&gmat[3][0]);
          Array<OneD, NekDouble> D2Dy(ptsx.num_elements(),&gmat[4][0]);
          Array<OneD, NekDouble> D2Dz(ptsx.num_elements(),&gmat[5][0]);
          Array<OneD, NekDouble> D3Dx(ptsx.num_elements(),&gmat[6][0]);
          Array<OneD, NekDouble> D3Dy(ptsx.num_elements(),&gmat[7][0]);
          Array<OneD, NekDouble> D3Dz(ptsx.num_elements(),&gmat[8][0]);
          
          int cnt=0; 
          int MaxIterations = 40;
          NekDouble Tol = 1e-12; // This is error*error; 
          
          F1 = F2 = F3 = 2000; // Starting value of Function
          
          while(cnt++ < MaxIterations)
          {
              //calculate the global point `corresponding to Lcoords
              xmap = m_xmap[0]->PhysEvaluate(Lcoords, ptsx);
              ymap = m_xmap[1]->PhysEvaluate(Lcoords, ptsy);
              zmap = m_xmap[2]->PhysEvaluate(Lcoords, ptsz);
              
              F1 = coords[0] - xmap;
              F2 = coords[1] - ymap;
              F3 = coords[2] - zmap;
              
              // stopping criterion
              if(F1*F1 + F2*F2 + F3*F3 < Tol)
              {
                  break;
              }
              
              //Interpolate derivative metric at Lcoords
              der1_x = m_xmap[0]->PhysEvaluate(Lcoords, D1Dx);
              der2_x = m_xmap[0]->PhysEvaluate(Lcoords, D1Dy);
              der3_z = m_xmap[0]->PhysEvaluate(Lcoords, D1Dz);
              der1_y = m_xmap[1]->PhysEvaluate(Lcoords, D2Dx);
              der2_y = m_xmap[1]->PhysEvaluate(Lcoords, D2Dy);                  
              der3_y = m_xmap[1]->PhysEvaluate(Lcoords, D2Dz);          
              der1_z = m_xmap[2]->PhysEvaluate(Lcoords, D3Dx);
              der2_z = m_xmap[2]->PhysEvaluate(Lcoords, D3Dy);                  
              der3_z = m_xmap[2]->PhysEvaluate(Lcoords, D3Dz);          
              
              Lcoords[0] = Lcoords[0] + der1_x*(coords[0]-xmap) + 
                  der1_y*(coords[1]-ymap) +  der1_z*(coords[2]-zmap);
              Lcoords[1] = Lcoords[1] + der2_x*(coords[0]-xmap) + 
                  der2_y*(coords[1]-ymap) +  der2_z*(coords[2]-zmap);
              Lcoords[2] = Lcoords[2] + der3_x*(coords[0]-xmap) + 
                  der3_y*(coords[1]-ymap) +  der3_z*(coords[2]-zmap);
          }
          
          if(cnt >= 40)
          {
              Lcoords[0] = Lcoords[1] = Lcoords[2] = 2.0;    
          }                        
      }
      
      
      /** 
       * @brief Put all quadrature information into face/edge structure and
       * backward transform.
       * 
       * Note verts, edges, and faces are listed according to anticlockwise
       * convention but points in _coeffs have to be in array format from left
       * to right.
       */
      void Geometry3D::v_FillGeom()
      {
          if (m_state == ePtsFilled)
              return;

          int i,j,k;

          for(i = 0; i < m_forient.size(); i++)
          {
              m_faces[i]->FillGeom();

              int nFaceCoeffs = (*m_faces[i])[0]->GetNcoeffs();
              Array<OneD, unsigned int> mapArray (nFaceCoeffs);
              Array<OneD,          int> signArray(nFaceCoeffs);
              
              if (m_forient[i] < 9)
              {
                  m_xmap[0]->GetFaceToElementMap(
                      i,m_forient[i],mapArray,signArray,
                      m_faces[i]->GetXmap(0)->GetEdgeNcoeffs(0),
                      m_faces[i]->GetXmap(0)->GetEdgeNcoeffs(1));
              }
              else
              {
                  m_xmap[0]->GetFaceToElementMap(
                      i,m_forient[i],mapArray,signArray,
                      m_faces[i]->GetXmap(0)->GetEdgeNcoeffs(1),
                      m_faces[i]->GetXmap(0)->GetEdgeNcoeffs(0));
              }

              for(j = 0; j < m_coordim; j++)
              {
                  const Array<OneD, const NekDouble> &coeffs = 
                      (*m_faces[i])[j]->GetCoeffs();

                  for(k = 0; k < nFaceCoeffs; k++)
                  {
                      double v = signArray[k] * coeffs[k];
                      (m_xmap[j]->UpdateCoeffs())[mapArray[k]] = v;
                  }
              }
          }

          for(i = 0; i < m_coordim; ++i)
          {
              m_xmap[i]->BwdTrans(m_xmap[i]->GetCoeffs (),
                                  m_xmap[i]->UpdatePhys());
          }

          m_state = ePtsFilled;
      }

        /**
         * Generate the geometry factors for this element.
         */
        void Geometry3D::v_GenGeomFactors(
                const Array<OneD, const LibUtilities::BasisSharedPtr> &tbasis)
        {
            if (m_geomFactorsState != ePtsFilled)
            {
                GeomType      Gtype  = eRegular;
                GeomShapeType GSType = eQuadrilateral;

                v_FillGeom();

                // check to see if expansions are linear
                for(int i = 0; i < m_coordim; ++i)
                {
                    if (m_xmap[i]->GetBasisNumModes(0) != 2 ||
                        m_xmap[i]->GetBasisNumModes(1) != 2 ||
                        m_xmap[i]->GetBasisNumModes(2) != 2)
                    {
                        Gtype = eDeformed;
                    }
                }

                m_geomFactors = MemoryManager<GeomFactors3D>::AllocateSharedPtr(
                    Gtype, m_coordim, m_xmap, tbasis);

                m_geomFactorsState = ePtsFilled;
            }
        }

      /** 
       * @brief Given local collapsed coordinate Lcoord return the value of
       * physical coordinate in direction i.
       */
      NekDouble Geometry3D::v_GetCoord(
          const int i, const Array<OneD, const NekDouble> &Lcoord)
      {
          ASSERTL1(m_state == ePtsFilled, "Geometry is not in physical space.");
          return m_xmap[i]->PhysEvaluate(Lcoord);
      }


      //---------------------------------------
      // Helper functions
      //---------------------------------------

      /**
       * @brief Return the co-ordinate mapping for dimension i.
       */
      StdRegions::StdExpansion3DSharedPtr Geometry3D::GetXmap(const int i)
      {
          return v_GetXmap(i);
      }
      StdRegions::StdExpansion3DSharedPtr Geometry3D::v_GetXmap(const int i)
      {
          return m_xmap[i];
      }

      /**
       * @brief Return the dimension of this element.
       */
      int Geometry3D::v_GetShapeDim() const 
      {
          return 3;
      }

      /**
       * @brief Return the vertex ID of vertex i.
       */
      int Geometry3D::v_GetVid(const int i) const
      {
          ASSERTL2(i >= 0 && i <= m_verts.size() - 1,
                   "Vertex ID must be between 0 and "+
                   boost::lexical_cast<string>(m_verts.size() - 1));
          return m_verts[i]->GetVid();
      }

      /**
       * @brief Return edge i of this element.
       */
      const SegGeomSharedPtr Geometry3D::v_GetEdge(int i) const
      {
          ASSERTL2(i >= 0 && i <= m_edges.size() - 1,
                   "Edge ID must be between 0 and "+
                   boost::lexical_cast<string>(m_edges.size() - 1));
          return m_edges[i];
      }

      /**
       * @brief Return the orientation of edge i in this element.
       */
      inline StdRegions::Orientation Geometry3D::v_GetEorient(
          const int i) const
      {
          ASSERTL2(i >= 0 && i <= m_edges.size() - 1,
                   "Edge ID must be between 0 and "+
                   boost::lexical_cast<string>(m_edges.size() - 1));
          return m_eorient[i];
      }

      /**
       * @brief Return the ID of edge i in this element.
       */
      int Geometry3D::v_GetEid(int i) const
      {
          ASSERTL2(i >= 0 && i <= m_edges.size() - 1,
                   "Edge ID must be between 0 and "+
                   boost::lexical_cast<string>(m_edges.size() - 1));
          return m_edges[i]->GetEid();
      }

      /**
       * @brief Return face i in this element.
       */
      const Geometry2DSharedPtr Geometry3D::v_GetFace(int i) const
      {
          ASSERTL2((i >=0) && (i <= 5),"Edge id must be between 0 and 4");
          return m_faces[i];
      }

      /**
       * @brief Return the orientation of face i in this element.
       */
      StdRegions::Orientation Geometry3D::v_GetFaceOrient(const int i) const
      {
          ASSERTL2(i >= 0 && i <= m_faces.size() - 1,
                   "Face ID must be between 0 and "+
                   boost::lexical_cast<string>(m_faces.size() - 1));
          return m_forient[i];
      }

      /**
       * @brief Return the ID of face i in this element.
       */
      int Geometry3D::v_GetFid(int i) const
      {
          ASSERTL2(i >= 0 && i <= m_faces.size() - 1,
                   "Face ID must be between 0 and "+
                   boost::lexical_cast<string>(m_faces.size() - 1));
          return m_faces[i]->GetFid();
      }

      /**
       * @brief Return the ID of this element.
       */
      int Geometry3D::v_GetEid() const 
      {
          return m_eid;
      }

      /**
       * @brief Return the j-th basis of the i-th co-ordinate dimension.
       */
      const LibUtilities::BasisSharedPtr Geometry3D::v_GetBasis(
          const int i, const int j)
      {
          return m_xmap[i]->GetBasis(j);
      }

      /**
       * @brief Return a reference to the physical space of co-ordinate
       * dimension i.
       */
      Array<OneD,NekDouble> &Geometry3D::v_UpdatePhys(const int i)
      {
          return m_xmap[i]->UpdatePhys();
      }

      /**
       * @brief Return the local ID of a given edge.
       * 
       * The local ID of an edge is a number between 0 and the number of edges
       * in this element. If the edge is not found, this function returns -1.
       */
      int Geometry3D::v_WhichEdge(SegGeomSharedPtr edge)
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
       * @brief Return the local ID of a given face.
       * 
       * The local ID of a face is a number between 0 and the number of faces
       * in this element. If the face is not found, this function returns -1.
       */
      int Geometry3D::v_WhichFace(Geometry2DSharedPtr face)
      {
          int i = 0;

          Geometry2DVector::iterator f;
          for (i = 0, f = m_faces.begin(); f != m_faces.end(); ++f,++i)
          {
              if (*f == face)
              {
                  break;
              }
          }
          return i;
      }

      //---------------------------------------
      // Element connection functions
      //---------------------------------------

      void Geometry3D::v_AddElmtConnected(int gvo_id, int locid)
      { 
          CompToElmt ee(gvo_id,locid);
          m_elmtmap.push_back(ee);
      }

      int Geometry3D::v_NumElmtConnected() const
      {
          return int(m_elmtmap.size());
      }

      bool Geometry3D::v_IsElmtConnected(int gvo_id, int locid) const
      {
          std::list<CompToElmt>::const_iterator def;
          CompToElmt ee(gvo_id,locid);

          def = find(m_elmtmap.begin(),m_elmtmap.end(),ee);

          // Found the element connectivity object in the list
          return (def != m_elmtmap.end());
      }

      void Geometry3D::v_SetOwnData()
      {
          m_owndata = true; 
      }

  }; //end of namespace
}; //end of namespace

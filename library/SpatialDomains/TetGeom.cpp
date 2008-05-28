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

// #include "pchSpatialDomains.h"
// 
// #include <StdRegions/StdRegions.hpp>
// #include <SpatialDomains/SpatialDomains.hpp>
// 
// #include <SpatialDomains/TriGeom.h>
// #include <SpatialDomains/TetGeom.h>
// 
// #include <SpatialDomains/GeomFactors.h>
// #include <SpatialDomains/Geometry3D.h>
// #include <SpatialDomains/MeshComponents.h>
// #include <SpatialDomains/TriFaceComponent.h>

#include "pchSpatialDomains.h"

#include <SpatialDomains/TetGeom.h>



namespace Nektar
{
    namespace SpatialDomains
    {

         TetGeom::TetGeom(const VertexComponentSharedPtr verts[], const SegGeomSharedPtr edges[], const TriGeomSharedPtr faces[],
                          const StdRegions::EdgeOrientation eorient[], const StdRegions::FaceOrientation forient[])
         {
            m_GeomShapeType = eTet;
 
            /// Copy the vert shared pointers.
            m_verts.insert(m_verts.begin(), verts, verts+TetGeom::kNverts);

            /// Copy the edge shared pointers.
            m_edges.insert(m_edges.begin(), edges, edges+TetGeom::kNedges);

            /// Copy the face shared pointers
            m_tfaces.insert(m_tfaces.begin(), faces, faces+TetGeom::kNfaces);

            for (int i=0; i<kNedges; ++i)
            {
                m_eorient[i] = eorient[i];
            }

            //TODO: check FaceOrientation for Tetrahedron case
            for (int j=0; j<kNfaces; ++j)
            {
               m_forient[j] = forient[j];
            }

            m_coordim = verts[0]->GetCoordim();
            ASSERTL0(m_coordim > 2,"Cannot call function with dim == 2");
        }

        
        TetGeom::TetGeom (const TriGeomSharedPtr faces[],  const StdRegions::FaceOrientation forient[])
        {
            m_GeomShapeType = eTet;
        }

        TetGeom::TetGeom(const SegGeomSharedPtr edges[], const StdRegions::EdgeOrientation eorient[])          
        {
            m_GeomShapeType = eTet;

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

        void TetGeom::GenGeomFactors(void)
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
          //  m_geomfactors = MemoryManager<GeomFactors>::AllocateSharedPtr(Gtype, m_coordim, m_xmap);
        }

       void TetGeom::FillGeom()
       {
         // check to see if geometry structure is already filled
            if(m_state == ePtsFilled)
            {
                return;
            }

            int i,j; 
            int order0 = m_xmap[0]->GetBasisNumModes(0);
            int order1 = m_xmap[0]->GetBasisNumModes(1);
            Array<OneD, const NekDouble> coef;
            StdRegions::StdExpMap Map,MapV;

            // set side 1 
            m_xmap[0]->MapTo((*m_edges[0])[0]->GetNcoeffs(),
                 (*m_edges[0])[0]->GetBasisType(0),
                 0,m_eorient[0],Map);

            for(i = 0; i < m_coordim; ++i)
            {
                coef  = (*m_edges[0])[i]->GetCoeffs();
                for(j = 0; j < order0; ++j)
                {
                    m_xmap[i]->SetCoeff(Map[j],coef[j]);
                }
            }

            // set Vertices into side 2 
            (*m_edges[1])[0]->MapTo(m_eorient[1],MapV);

            for(i = 0; i < m_coordim; ++i)
            {
                (*m_edges[1])[i]->SetCoeff(MapV[0],(*m_verts[1])[i]);
                (*m_edges[1])[i]->SetCoeff(MapV[1],(*m_verts[2])[i]);
            }

            // set Vertices into side 3
            (*m_edges[2])[0]->MapTo(m_eorient[2],MapV);

            for(i = 0; i < m_coordim; ++i)
            {
                (*m_edges[2])[i]->SetCoeff(MapV[0],(*m_verts[0])[i]);
                (*m_edges[2])[i]->SetCoeff(MapV[1],(*m_verts[2])[i]);
            }

            // set side 2
            m_xmap[0]->MapTo((*m_edges[1])[0]->GetNcoeffs(),
                 (*m_edges[1])[0]->GetBasisType(0),
                 1,(m_eorient[1]),Map);
        
            for(i = 0; i < m_coordim; ++i)
            {
                coef  = (*m_edges[1])[i]->GetCoeffs();
                for(j = 0; j < order1; ++j)
                {
                    m_xmap[i]->SetCoeff(Map[j],coef[j]);
                }
            }

            // set side 3
            m_xmap[0]->MapTo((*m_edges[2])[0]->GetNcoeffs(),
                 (*m_edges[2])[0]->GetBasisType(0),
                 2,(m_eorient[2]),Map);

            for(i = 0; i < m_coordim; ++i)
            {
                coef  = (*m_edges[2])[i]->GetCoeffs();
                for(j = 0; j < order0; ++j)
                {
                    m_xmap[i]->SetCoeff(Map[j],coef[j]);
                }
            }

            for(i = 0; i < m_coordim; ++i)
            {
                m_xmap[i]->BwdTrans(m_xmap[i]->GetCoeffs(),
                                    m_xmap[i]->UpdatePhys());
            }

            m_state = ePtsFilled;

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
                NekDouble V = e30.dot(cp1020);
                NekDouble beta  = er0.dot(cp2030) / V;
                NekDouble gamma = er0.dot(cp3010) / V;
                NekDouble delta = er0.dot(cp1020) / V;

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

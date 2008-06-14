////////////////////////////////////////////////////////////////////////////////
//
//  File:  $Source: /usr/sci/projects/Nektar/cvs/Nektar++/library/SpatialDomains/HexGeom.cpp,v $
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

#include <SpatialDomains/HexGeom.h>

namespace Nektar
{
    namespace SpatialDomains
    {
        HexGeom::HexGeom(const QuadGeomSharedPtr faces[],  const StdRegions::FaceOrientation forient[])
        {
            m_GeomShapeType = eHexahedron;

            /// Copy the face shared pointers
            m_qfaces.insert(m_qfaces.begin(), faces, faces+HexGeom::kNfaces);

            for (int j=0; j<kNfaces; ++j)
            {
               m_forient[j] = forient[j];
            }

            m_coordim = faces[0]->GetEdge(0)->GetVertex(0)->GetCoordim();
            ASSERTL0(m_coordim > 2,"Cannot call function with dim == 2");
        }

        HexGeom::HexGeom(const VertexComponentSharedPtr verts[], const SegGeomSharedPtr edges[], const QuadGeomSharedPtr faces[],
                         const StdRegions::EdgeOrientation eorient[], const StdRegions::FaceOrientation forient[])
         {
            m_GeomShapeType = eHexahedron;
 
            /// Copy the vert shared pointers.
            m_verts.insert(m_verts.begin(), verts, verts+HexGeom::kNverts);

            /// Copy the edge shared pointers.
            m_edges.insert(m_edges.begin(), edges, edges+HexGeom::kNedges);

            /// Copy the face shared pointers
            m_qfaces.insert(m_qfaces.begin(), faces, faces+HexGeom::kNfaces);

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

        HexGeom::~HexGeom()
        {
        }

        void HexGeom::AddElmtConnected(int gvo_id, int locid)
        {
            CompToElmt ee(gvo_id,locid);
            m_elmtmap.push_back(ee);
        }

        int HexGeom::NumElmtConnected() const
        {
            return int(m_elmtmap.size());
        }

        bool HexGeom::IsElmtConnected(int gvo_id, int locid) const
        {
            std::list<CompToElmt>::const_iterator def;
            CompToElmt ee(gvo_id,locid);

            def = find(m_elmtmap.begin(),m_elmtmap.end(),ee);

            // Found the element connectivity object in the list
            return (def != m_elmtmap.end());
        }

        /** given local collapsed coordinate Lcoord return the value of
        physical coordinate in direction i **/


        NekDouble HexGeom::GetCoord(const int i, 
            const Array<OneD, const NekDouble> &Lcoord)
        {
            ASSERTL1(m_state == ePtsFilled,
                "Goemetry is not in physical space");

            return m_xmap[i]->PhysEvaluate(Lcoord);
        }

        // Set up GeoFac for this geometry using Coord quadrature distribution

        void HexGeom::GenGeomFactors(void)
        {
			// TODO: Insert code here.
		}


          /** \brief put all quadrature information into edge structure 
        and backward transform 

        Note verts, edges, and faces are listed according to anticlockwise
        convention but points in _coeffs have to be in array format from
        left to right.

        */
		void HexGeom::FillGeom()
		{
           // check to see if geometry structure is already filled
            if(m_state != ePtsFilled)
            {
                int i,j,k;

                int nFaceCoeffs = m_xmap[0]->GetFaceNcoeffs(0); //TODO: implement GetFaceNcoeffs()
                Array<OneD, unsigned int> mapArray (nFaceCoeffs);
                Array<OneD, int>          signArray(nFaceCoeffs);

                 for(i = 0; i < kNfaces; i++)
                {
                    m_qfaces[i]->FillGeom();
                    
                    //TODO: implement GetFaceToElementMap()
                    //  m_xmap[0]->GetFaceToElementMap(i,m_forient[i],mapArray,signArray); 
                    
                    nFaceCoeffs = m_xmap[0]->GetFaceNcoeffs(i);

                    for(j = 0 ; j < m_coordim; j++)
                    {
                        for(k = 0; k < nFaceCoeffs; k++)
                        {
                           //TODO : insert code here
                            //  (m_xmap[j]->UpdateCoeffs())[ mapArray[k] ] = signArray[k]*
                            //  ((*m_tfaces[i])[j]->GetCoeffs())[k];
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

		void HexGeom::GetLocCoords(const Array<OneD, const NekDouble> &coords, Array<OneD, NekDouble> &Lcoords)
		{
			// TODO: Insert code here.
		}

    }; //end of namespace
}; //end of namespace

//
// $Log: HexGeom.cpp,v $
// Revision 1.8  2008/06/12 21:22:43  delisi
// Added method stubs for GenGeomFactors, FillGeom, and GetLocCoords.
//
// Revision 1.7  2008/05/29 19:02:23  delisi
// Renamed eHex to eHexahedron.
//
// Revision 1.6  2008/05/28 21:52:27  jfrazier
// Added GeomShapeType initialization for the different shapes.
//
// Revision 1.5  2008/05/12 17:28:26  ehan
// Added virtual functions
//
// Revision 1.4  2008/04/06 06:00:37  bnelson
// Changed ConstArray to Array<const>
//
// Revision 1.3  2008/02/08 23:05:28  jfrazier
// More work on 3D components.
//
// Revision 1.2  2007/07/20 02:15:08  bnelson
// Replaced boost::shared_ptr with Nektar::ptr
//
// Revision 1.1  2006/05/04 18:59:00  kirby
// *** empty log message ***
//
// Revision 1.11  2006/04/09 02:08:35  jfrazier
// Added precompiled header.
//
// Revision 1.10  2006/03/12 07:42:02  sherwin
//
// Updated member names and StdRegions call. Still has not been compiled
//
// Revision 1.9  2006/02/19 01:37:33  jfrazier
// Initial attempt at bringing into conformance with the coding standard.  Still more work to be done.  Has not been compiled.
//
//

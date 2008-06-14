////////////////////////////////////////////////////////////////////////////////
//
//  File:  $Source: /usr/sci/projects/Nektar/cvs/Nektar++/library/SpatialDomains/PrismGeom.cpp,v $
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

#include <SpatialDomains/PrismGeom.h>

namespace Nektar
{
    namespace SpatialDomains
    {
        PrismGeom::PrismGeom(const TriGeomSharedPtr tfaces[], const QuadGeomSharedPtr qfaces[], const StdRegions::FaceOrientation forient[])
        {
            m_GeomShapeType = ePrism;

            /// Copy the quad face shared pointers
            m_qfaces.insert(m_qfaces.begin(), qfaces, qfaces+PrismGeom::kNfaces);

            /// Copy the triangle face shared pointers
            m_tfaces.insert(m_tfaces.begin(), tfaces, tfaces+PrismGeom::kNfaces);
            
            for (int j=0; j<kNfaces; ++j)
            {
               m_forient[j] = forient[j];
            }

            m_coordim = tfaces[0]->GetEdge(0)->GetVertex(0)->GetCoordim();
            ASSERTL0(m_coordim > 2,"Cannot call function with dim == 2");

        }

        PrismGeom::PrismGeom(const VertexComponentSharedPtr verts[], const SegGeomSharedPtr edges[], const TriGeomSharedPtr tfaces[],
                             const QuadGeomSharedPtr qfaces[],const StdRegions::EdgeOrientation eorient[],
                             const StdRegions::FaceOrientation forient[])
        {

             m_GeomShapeType = ePrism;
 
            /// Copy the vert shared pointers.
            m_verts.insert(m_verts.begin(), verts, verts+PrismGeom::kNverts);

            /// Copy the edge shared pointers.
            m_edges.insert(m_edges.begin(), edges, edges+PrismGeom::kNedges);

            /// Copy the quad face shared pointers
            m_qfaces.insert(m_qfaces.begin(), qfaces, qfaces+PrismGeom::kNfaces);

            /// Copy the triangle face shared pointers
            m_tfaces.insert(m_tfaces.begin(), tfaces, tfaces+PrismGeom::kNfaces);
            
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

        PrismGeom::~PrismGeom()
        {
        }

        void PrismGeom::AddElmtConnected(int gvo_id, int locid)
        {
            CompToElmt ee(gvo_id,locid);
            m_elmtmap.push_back(ee);
        }


        int PrismGeom::NumElmtConnected() const
        {
            return int(m_elmtmap.size());
        }


        bool PrismGeom::IsElmtConnected(int gvo_id, int locid) const
        {
            std::list<CompToElmt>::const_iterator def;
            CompToElmt ee(gvo_id,locid);

            def = find(m_elmtmap.begin(),m_elmtmap.end(),ee);

            // Found the element connectivity object in the list
            return (def != m_elmtmap.end());
        }

        /** given local collapsed coordinate Lcoord return the value of
        physical coordinate in direction i **/
        NekDouble PrismGeom::GetCoord(const int i, 
                                          const Array<OneD, const NekDouble> &Lcoord)
        {
            ASSERTL1(m_state == ePtsFilled,
                "Goemetry is not in physical space");

            return m_xmap[i]->PhysEvaluate(Lcoord);
        }

        // Set up GeoFac for this geometry using Coord quadrature distribution

        void PrismGeom::GenGeomFactors(void)
        {
			// TODO: Insert code here.
		}



          /** \brief put all quadrature information into edge structure 
        and backward transform 

        Note verts, edges, and faces are listed according to anticlockwise
        convention but points in _coeffs have to be in array format from
        left to right.

        */

        void PrismGeom::FillGeom()
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
                    m_tfaces[i]->FillGeom();
                    
                    //TODO: implement GetFaceToElementMap()
                    //  m_xmap[0]->GetFaceToElementMap(i,m_forient[i],mapArray,signArray); 
                    
                    nFaceCoeffs = m_xmap[0]->GetFaceNcoeffs(i);

                    for(j = 0 ; j < m_coordim; j++)
                    {
                        for(k = 0; k < nFaceCoeffs; k++)
                        {
                           //TODO : insert code here
//                             (m_xmap[j]->UpdateCoeffs())[ mapArray[k] ] = signArray[k]*
//                                        ((*m_tfaces[i])[j]->GetCoeffs())[k];
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

		void PrismGeom::GetLocCoords(const Array<OneD, const NekDouble> &coords, Array<OneD, NekDouble> &Lcoords)
		{
			// TODO: Insert code here.
		}

    }; //end of namespace
}; //end of namespace

//
// $Log: PrismGeom.cpp,v $
// Revision 1.6  2008/06/12 21:22:55  delisi
// Added method stubs for GenGeomFactors, FillGeom, and GetLocCoords.
//
// Revision 1.5  2008/06/11 16:10:12  delisi
// Added the 3D reader.
//
// Revision 1.4  2008/05/28 21:52:27  jfrazier
// Added GeomShapeType initialization for the different shapes.
//
// Revision 1.3  2008/04/06 06:00:38  bnelson
// Changed ConstArray to Array<const>
//
// Revision 1.2  2008/02/08 23:05:52  jfrazier
// More work on 3D components.
//
// Revision 1.1  2006/05/04 18:59:02  kirby
// *** empty log message ***
//
// Revision 1.10  2006/04/09 02:08:35  jfrazier
// Added precompiled header.
//
// Revision 1.9  2006/03/13 18:20:03  sherwin
//
// Fixed error in ResetGmat:
//
// Revision 1.8  2006/02/19 01:37:34  jfrazier
// Initial attempt at bringing into conformance with the coding standard.  Still more work to be done.  Has not been compiled.
//
//

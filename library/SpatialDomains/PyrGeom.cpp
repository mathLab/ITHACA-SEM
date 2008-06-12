////////////////////////////////////////////////////////////////////////////////
//
//  File:  $Source: /usr/sci/projects/Nektar/cvs/Nektar++/library/SpatialDomains/PyrGeom.cpp,v $
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

#include <SpatialDomains/PyrGeom.h>

namespace Nektar
{
    namespace SpatialDomains
    {

        PyrGeom::PyrGeom(const TriGeomSharedPtr tfaces[], const QuadGeomSharedPtr qfaces[], const StdRegions::FaceOrientation forient[])
        {
            m_GeomShapeType = ePyramid;
        }

        PyrGeom::~PyrGeom()
        {
        }
        void PyrGeom::AddElmtConnected(int gvo_id, int locid)
        {
            CompToElmt ee(gvo_id,locid);
            m_elmtmap.push_back(ee);
        }


        int PyrGeom::NumElmtConnected() const
        {
            return int(m_elmtmap.size());
        }


        bool PyrGeom::IsElmtConnected(int gvo_id, int locid) const
        {
            std::list<CompToElmt>::const_iterator def;
            CompToElmt ee(gvo_id,locid);

            def = find(m_elmtmap.begin(),m_elmtmap.end(),ee);

            // Found the element connectivity object in the list
            return (def != m_elmtmap.end());
        }

        /** given local collapsed coordinate Lcoord return the value of
        physical coordinate in direction i **/


        NekDouble PyrGeom::GetCoord(const int i, 
                                          const Array<OneD, const NekDouble> &Lcoord)
        {
            ASSERTL1(m_state == ePtsFilled,
                "Goemetry is not in physical space");

            return m_xmap[i]->PhysEvaluate(Lcoord);
        }

        // Set up GeoFac for this geometry using Coord quadrature distribution

        void PyrGeom::GenGeomFactors(void)
        {
			// TODO: Insert code here.
		}

		void PyrGeom::FillGeom()
		{
			// check to see if geometry structure is already filled
			if(m_state == ePtsFilled)
			{
				return;
			}

			// TODO: Insert code here.

			m_state = ePtsFilled;

		}

		void PyrGeom::GetLocCoords(const Array<OneD, const NekDouble> &coords, Array<OneD, NekDouble> &Lcoords)
		{
			// TODO: Insert code here.
		}

    }; //end of namespace
}; //end of namespace

//
// $Log: PyrGeom.cpp,v $
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
// Revision 1.9  2006/03/12 11:06:39  sherwin
//
// First complete copy of code standard code but still not compile tested
//
// Revision 1.8  2006/02/19 01:37:34  jfrazier
// Initial attempt at bringing into conformance with the coding standard.  Still more work to be done.  Has not been compiled.
//
//

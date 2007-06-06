////////////////////////////////////////////////////////////////////////////////
//
//  File:  QuadFaceComponent.cpp
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

#include <SpatialDomains/QuadFaceComponent.h>

namespace Nektar
{
    namespace SpatialDomains
    {
        QuadFaceComponent::QuadFaceComponent()
        {
        }

        QuadFaceComponent::QuadFaceComponent(const int coordim):
            Geometry2D(coordim),
            m_xmap(coordim)
        {
            const LibUtilities::BasisKey B(LibUtilities::eModified_A, 2,
                LibUtilities::PointsKey(3,LibUtilities::eGaussLobattoLegendre));


            for(int i = 0; i < m_coordim; ++i)
            {
                m_xmap[i] = MemoryManager<StdRegions::StdQuadExp>::AllocateSharedPtr(B,B);  
            }
        }

        QuadFaceComponent::~QuadFaceComponent()
        {
        }

        QuadFaceComponent::QuadFaceComponent(const QuadFaceComponent &T)
        {
            m_fid = T.m_fid;

            std::list<CompToElmt>::const_iterator def;
            for(def = T.m_elmtmap.begin(); def != T.m_elmtmap.end(); def++)
            {
                m_elmtmap.push_back(*def);    
            }
        }

        void QuadFaceComponent::AddElmtConnected(int gvo_id, int locid)
        {
            CompToElmt ee(gvo_id,locid);
            m_elmtmap.push_back(ee);
        }

        int QuadFaceComponent::NumElmtConnected() const
        {
            return int(m_elmtmap.size());
        }

        bool QuadFaceComponent::IsElmtConnected(int gvo_id, int locid) const
        {
            std::list<CompToElmt>::const_iterator def;
            CompToElmt ee(gvo_id,locid);

            def = find(m_elmtmap.begin(),m_elmtmap.end(),ee);

            // Found the element connectivity object in the list
            if(def != m_elmtmap.end())
            {
                return(true);
            }

            return(false);
        }

        double QuadFaceComponent::GetCoord(const int i, 
                                           const ConstArray<OneD,NekDouble> &Lcoord)
        {
            ASSERTL1(m_state == ePtsFilled,
                "Goemetry is not in physical space");

            return m_xmap[i]->PhysEvaluate(Lcoord);
        }
    }; //end of namespace
}; //end of namespace

//
// $Log: QuadFaceComponent.cpp,v $
// Revision 1.2  2007/05/28 21:48:42  sherwin
// Update for 2D functionality
//
// Revision 1.1  2006/05/04 18:59:02  kirby
// *** empty log message ***
//
// Revision 1.11  2006/04/09 02:08:35  jfrazier
// Added precompiled header.
//
// Revision 1.10  2006/04/05 18:08:56  sherwin
// Updates for SpatialDomains
//
// Revision 1.9  2006/03/25 00:58:29  jfrazier
// Many changes dealing with fundamental structure and reading/writing.
//
// Revision 1.8  2006/03/12 14:20:43  sherwin
//
// First compiling version of SpatialDomains and associated modifications
//
// Revision 1.7  2006/03/12 11:06:39  sherwin
//
// First complete copy of code standard code but still not compile tested
//
// Revision 1.6  2006/02/19 01:37:34  jfrazier
// Initial attempt at bringing into conformance with the coding standard.  Still more work to be done.  Has not been compiled.
//
//

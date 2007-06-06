////////////////////////////////////////////////////////////////////////////////
//
//  File: EdgeComponent.cpp
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

#include <SpatialDomains/EdgeComponent.h>

namespace Nektar
{
    namespace SpatialDomains
    {
        EdgeComponent::EdgeComponent()
        {
        }

        EdgeComponent::EdgeComponent(int id, const int coordim): 
        Geometry1D(coordim),
            m_xmap(coordim)
        {

            const LibUtilities::BasisKey B(LibUtilities::eModified_A, 2,
                LibUtilities::PointsKey(3,LibUtilities::eGaussLobattoLegendre));
            m_eid = id;

            for(int i = 0; i < m_coordim; ++i)
            {
                m_xmap[i] = MemoryManager<StdRegions::StdSegExp>::AllocateSharedPtr(B);
            }
        }

        EdgeComponent::EdgeComponent(int id, const int coordim,
            VertexComponentSharedPtr vertex[]): 
        Geometry1D(coordim),
            m_vertex(2) //always have two vertices per edge
        {
            m_eid = id;

            if (coordim > 0)
            {
                const LibUtilities::BasisKey B(LibUtilities::eModified_A, 2,
                    LibUtilities::PointsKey(3,LibUtilities::eGaussLobattoLegendre));

                m_xmap = Array<OneD, StdRegions::StdExpansion1DSharedPtr>(m_coordim);

                for(int i = 0; i < m_coordim; ++i)
                {
                    m_xmap[i] = MemoryManager<StdRegions::StdSegExp>::AllocateSharedPtr(B);
                }
            }

            m_vertex[0] = vertex[0];
            m_vertex[1] = vertex[1];
        }

        EdgeComponent::EdgeComponent(int id, const int coordim, 
            const int order, const int nquad):
        Geometry1D(coordim)
        {


            const LibUtilities::BasisKey B(LibUtilities::eModified_A, order,
                LibUtilities::PointsKey(nquad,LibUtilities::eGaussLobattoLegendre));

            m_xmap = Array<OneD, StdRegions::StdExpansion1DSharedPtr>(m_coordim);

            for(int i = 0; i < m_coordim; ++i)
            {
                m_xmap[i] = MemoryManager<StdRegions::StdSegExp>::AllocateSharedPtr(B);
            }

            m_state = eNotFilled;
        }


        EdgeComponent::~EdgeComponent()
        {
        }


        EdgeComponent::EdgeComponent(const EdgeComponent &T)
        {

            m_eid = T.m_eid;

            std::list<CompToElmt>::const_iterator def;
            for(def = T.m_elmtmap.begin(); def != T.m_elmtmap.end(); def++)
            {
                m_elmtmap.push_back(*def);    
            }
        }

        /** given local collapsed coordinate Lcoord return the value of
        physical coordinate in direction i **/

        NekDouble EdgeComponent::GetCoord(const int i, 
                                 const ConstArray<OneD, NekDouble> &Lcoord) 
        {

            ASSERTL1(m_state == ePtsFilled, "Goemetry is not in physical space");

            return m_xmap[i]->PhysEvaluate(Lcoord);
        }

        void EdgeComponent::AddElmtConnected(int gvo_id, int locid)
        {
            CompToElmt ee(gvo_id,locid);
            m_elmtmap.push_back(ee);
        }


        int EdgeComponent::NumElmtConnected() const
        {
            return int(m_elmtmap.size());
        }


        bool EdgeComponent::IsElmtConnected(int gvo_id, int locid) const
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

        /// \brief Get the orientation of edge1.
        ///
        /// If edge1 is connected to edge2 in the same direction as
        /// the points comprising edge1 then it is forward, otherwise
        /// it is backward.  
        ///
        /// For example, assume edge1 is comprised of points 1 and 2,
        /// and edge2 is comprised of points 2 and 3, then edge is
        /// forward.
        ///
        /// If edge1 is comprised of points 2 and 1 and edge2 is
        /// comprised of points 3 and 2, then edge1 is backward.

        StdRegions::EdgeOrientation EdgeComponent::GetEdgeOrientation(const EdgeComponent &edge1,  const EdgeComponent &edge2)
        {
            StdRegions::EdgeOrientation returnval = StdRegions::eForwards;

            /// Backward direction.  Vertex 0 is connected to edge 2.
            if ((*edge1.GetVertex(0) == *edge2.GetVertex(0)) || 
                (*edge1.GetVertex(0) == *edge2.GetVertex(1)))
            {
                returnval = StdRegions::eBackwards;
            }
            // Not forward either, then we have a problem.
            else if ((*edge1.GetVertex(1) != *edge2.GetVertex(0)) && 
                (*edge1.GetVertex(1) != *edge2.GetVertex(1)))
            {
                std::ostrstream errstrm;
                errstrm << "Connected edges do not share a vertex. Edges ";
                errstrm << edge1.GetEid() << ", " << edge2.GetEid();
                ASSERTL0(false, errstrm.str());
            }

            return returnval;
        }
    }; //end of namespace
}; //end of namespace

/** 
*    $Log: EdgeComponent.cpp,v $
*    Revision 1.18  2007/05/28 08:35:25  sherwin
*    Updated for localregions up to Project1D
*
*    Revision 1.17  2007/05/17 18:45:24  jfrazier
*    Minor changes to accommodate Array class.
*
*    Revision 1.16  2007/04/08 03:34:48  jfrazier
*    Updated to compile with SharedArray.  This has not been converted to SharedArray, just made to work with others that have been converted.
*
*    Revision 1.15  2007/04/04 21:49:24  sherwin
*    Update for SharedArray
*
*    Revision 1.14  2007/03/25 15:48:21  sherwin
*    UPdate LocalRegions to take new NekDouble and shared_array formats. Added new Demos
*
*    Revision 1.13  2007/03/14 21:24:08  sherwin
*    Update for working version of MultiRegions up to ExpList1D
*
*    Revision 1.12  2007/02/19 08:06:24  sherwin
*    Modified files to be consistent with new StdRegions prototypes and turned off 2D & 3D Calls.
*
*    Revision 1.11  2006/08/17 21:02:50  jfrazier
*    Eliminated bug where an array was indexed after it was deleted.
*
*    Revision 1.10  2006/08/16 23:34:42  jfrazier
*    *** empty log message ***
*
*    Revision 1.9  2006/07/02 17:16:17  sherwin
*
*    Modifications to make MultiRegions work for a connected domain in 2D (Tris)
*
*    Revision 1.8  2006/06/02 18:48:40  sherwin
*    Modifications to make ProjectLoc2D run bit there are bus errors for order > 3
*
*    Revision 1.7  2006/05/30 14:00:04  sherwin
*    Updates to make MultiRegions and its Demos work
*
*    Revision 1.6  2006/05/23 20:19:58  jfrazier
*    Added #pragma to show where problems currently exist.
*
*    Revision 1.5  2006/05/23 19:56:33  jfrazier
*    These build and run, but the expansion pieces are commented out
*    because they would not run.
*
*    Revision 1.4  2006/05/16 20:12:59  jfrazier
*    Minor fixes to correct bugs.
*
*    Revision 1.3  2006/05/09 13:37:01  jfrazier
*    Removed duplicate definition of shared vertex pointer.
*
*    Revision 1.2  2006/05/06 20:36:16  sherwin
*    Modifications to get LocalRegions/Project1D working
*
*    Revision 1.1  2006/05/04 18:58:59  kirby
*    *** empty log message ***
*
*    Revision 1.14  2006/04/11 23:18:11  jfrazier
*    Completed MeshGraph2D for tri's and quads.  Not thoroughly tested.
*
*    Revision 1.13  2006/04/09 02:08:34  jfrazier
*    Added precompiled header.
*
*    Revision 1.12  2006/04/04 23:12:37  jfrazier
*    More updates to readers.  Still work to do on MeshGraph2D to store tris and quads.
*
*    Revision 1.11  2006/03/25 00:58:28  jfrazier
*    Many changes dealing with fundamental structure and reading/writing.
*
*    Revision 1.10  2006/03/12 14:20:42  sherwin
*
*    First compiling version of SpatialDomains and associated modifications
*
*    Revision 1.9  2006/03/12 11:06:39  sherwin
*
*    First complete copy of code standard code but still not compile tested
*
*    Revision 1.8  2006/03/12 07:42:02  sherwin
*
*    Updated member names and StdRegions call. Still has not been compiled
*
**/

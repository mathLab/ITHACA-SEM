////////////////////////////////////////////////////////////////////////////////
//
//  File:  EdgeComponent.h
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
//  Description:  Specification of the Edge Component Class.
//
//
////////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_SPATIALDOMAINS_EDGECOMPONENT_H
#define NEKTAR_SPATIALDOMAINS_EDGECOMPONENT_H

#include <SpatialDomains/SpatialDomains.hpp>
#include <SpatialDomains/MeshGraph.h>
#include <SpatialDomains/MeshComponents.h>
#include <SpatialDomains/Geometry1D.h>

#include <StdRegions/StdSegExp.h>
#include <LibUtilities/Foundations/Basis.h>
namespace Nektar
{
    namespace SpatialDomains
    {
        class EdgeComponent: public Geometry1D
        {

        public:
            EdgeComponent();
            EdgeComponent(int id, const int coordim);
            EdgeComponent(int id, const int coordim, VertexComponentSharedPtr vertex[]);
            EdgeComponent(int id, const int coordim, const int order, const int nquad);

            ~EdgeComponent();
            EdgeComponent(const EdgeComponent &T);

            void AddElmtConnected(int gvoId, int locId);
            int NumElmtConnected() const;
            bool IsElmtConnected(int gvoId, int locId) const;

            inline int GetEid() const 
            {
                return m_eid;
            }

            inline const LibUtilities::BasisSharedPtr GetBasis(int i, int j) 
            {
                return m_xmap[i]->GetBasis(j);
            }

            inline const StdRegions::StdExpansion1DSharedPtr &GetXmap(int i)
            {
                return m_xmap[i];
            }

            inline NekDoubleSharedArray &UpdatePhys(int i){
                return m_xmap[i]->UpdatePhys();
            }

            inline VertexComponentSharedPtr GetVertex(int i) const
            {
                VertexComponentSharedPtr returnval;

                if (i >= 0 && i < int(m_vertex.size()))
                {
                    returnval = m_vertex[i];
                }

                return returnval;
            }


            StdRegions::StdExpansion1DSharedPtr operator[](const int i) const
            {
                if((i>=0)&& (i<m_coordim))
                {
                    return m_xmap[i];
                }

                ErrorUtil::Error(ErrorUtil::efatal,__FILE__,__LINE__,
                    "Invalid Index used in [] operator");
                return m_xmap[0]; //should never be reached
            }

            NekDouble GetCoord(const int i, const NekDoubleSharedArray &Lcoord);
                    
            /// \brief Get the orientation of edge1.
            ///
            /// Since both edges are passed, it does
            /// not need any information from the EdgeComponent instance.
            static StdRegions::EdgeOrientation GetEdgeOrientation(const EdgeComponent &edge1,
                const EdgeComponent &edge2);

        protected:
            int m_eid;
            std::list<CompToElmt> m_elmtmap;

	    SharedArray<StdRegions::StdExpansion1DSharedPtr> m_xmap;

            VertexVector m_vertex;
        private:
        };

        typedef boost::shared_ptr<EdgeComponent> EdgeComponentSharedPtr;
        typedef std::vector<EdgeComponentSharedPtr> EdgeComponentVector; 
    }//end of namespace
}//end of namespace

#endif //NEKTAR_SPATIALDOMAINS_EDGECOMPONENT_H

//
// $Log: EdgeComponent.h,v $
// Revision 1.7  2007/03/29 19:23:59  bnelson
// Replaced boost::shared_array with SharedArray
//
// Revision 1.6  2007/03/25 15:48:22  sherwin
// UPdate LocalRegions to take new NekDouble and shared_array formats. Added new Demos
//
// Revision 1.5  2007/03/02 12:01:58  sherwin
// Update for working version of LocalRegions/Project1D
//
// Revision 1.4  2007/02/19 08:06:25  sherwin
// Modified files to be consistent with new StdRegions prototypes and turned off 2D & 3D Calls.
//
// Revision 1.3  2006/05/09 13:37:01  jfrazier
// Removed duplicate definition of shared vertex pointer.
//
// Revision 1.2  2006/05/06 20:36:16  sherwin
// Modifications to get LocalRegions/Project1D working
//
// Revision 1.1  2006/05/04 18:58:59  kirby
// *** empty log message ***
//
// Revision 1.17  2006/05/02 21:21:11  sherwin
// Corrected libraries to compile new version of spatialdomains and demo Graph1D
//
// Revision 1.16  2006/04/11 23:18:11  jfrazier
// Completed MeshGraph2D for tri's and quads.  Not thoroughly tested.
//
// Revision 1.15  2006/04/09 02:08:34  jfrazier
// Added precompiled header.
//
// Revision 1.14  2006/04/04 23:12:37  jfrazier
// More updates to readers.  Still work to do on MeshGraph2D to store tris and quads.
//
// Revision 1.13  2006/03/25 00:58:28  jfrazier
// Many changes dealing with fundamental structure and reading/writing.
//
// Revision 1.12  2006/03/12 14:20:42  sherwin
//
// First compiling version of SpatialDomains and associated modifications
//
// Revision 1.11  2006/03/12 11:06:39  sherwin
//
// First complete copy of code standard code but still not compile tested
//
// Revision 1.10  2006/03/12 07:42:02  sherwin
//
// Updated member names and StdRegions call. Still has not been compiled
//
// Revision 1.9  2006/03/04 20:26:04  bnelson
// Added comments after #endif.
//
// Revision 1.8  2006/02/26 21:19:42  bnelson
// Fixed a variety of compiler errors caused by updates to the coding standard.
//
// Revision 1.7  2006/02/19 01:37:32  jfrazier
// Initial attempt at bringing into conformance with the coding standard.  Still more work to be done.  Has not been compiled.
//
//



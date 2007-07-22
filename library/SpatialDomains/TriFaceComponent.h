////////////////////////////////////////////////////////////////////////////////
//
//  File:  TriFaceComponent.h
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
#ifndef NEKTAR_SPATIALDOMAINS_TRIFACECOMPONENT_H
#define NEKTAR_SPATIALDOMAINS_TRIFACECOMPONENT_H

#include <SpatialDomains/SpatialDomains.hpp>
#include <SpatialDomains/MeshGraph.h>
#include <SpatialDomains/MeshComponents.h>
#include <SpatialDomains/Geometry2D.h>

#include <StdRegions/StdTriExp.h>

#include <list>

namespace Nektar
{
    namespace SpatialDomains
    {
        class TriFaceComponent: public Geometry2D
        {
        public:
            TriFaceComponent();
            TriFaceComponent(const int coordim);
            ~TriFaceComponent();
            TriFaceComponent(const TriFaceComponent &T);

            void AddElmtConnected(int gvo_id, int locid);
            int  NumElmtConnected() const;
            bool IsElmtConnected(int gvo_id, int locid) const;

            inline int GetFid() const 
        {
          return m_fid;
        }

            inline int GetCoordDim() const 
        {
          return m_coordim;
        }

            inline const LibUtilities::BasisSharedPtr GetBasis(const int i, const int j)
        {
                return m_xmap[i]->GetBasis(j);
            }
        
            inline Array<OneD,NekDouble> &UpdatePhys(const int i)
            {
                return m_xmap[i]->UpdatePhys();
            }

            NekDouble GetCoord(const int i, const ConstArray<OneD,NekDouble> &Lcoord);

        protected:
            int m_fid;
            bool m_ownverts;
            std::list<CompToElmt> m_elmtmap;

        Array<OneD, StdRegions::StdExpansion2DSharedPtr> m_xmap;
        private:
        };

        typedef boost::shared_ptr<TriFaceComponent> TriFaceComponentSharedPtr;
        typedef std::vector<TriFaceComponentSharedPtr> TriFaceComponentVector; 
    }; //end of namespace
}; //end of namespace

#endif //NEKTAR_SPATIALDOMAINS_TRIFACECOMPONENT_H

//
// $Log: TriFaceComponent.h,v $
// Revision 1.4  2007/07/20 02:15:09  bnelson
// Replaced boost::shared_ptr with Nektar::ptr
//
// Revision 1.3  2007/06/06 15:15:21  pvos
// Some minor updates for 2D routines
//
// Revision 1.2  2007/05/28 21:48:42  sherwin
// Update for 2D functionality
//
// Revision 1.1  2006/05/04 18:59:05  kirby
// *** empty log message ***
//
// Revision 1.8  2006/03/12 14:20:44  sherwin
//
// First compiling version of SpatialDomains and associated modifications
//
// Revision 1.7  2006/03/12 11:06:40  sherwin
//
// First complete copy of code standard code but still not compile tested
//
// Revision 1.6  2006/03/04 20:26:05  bnelson
// Added comments after #endif.
//
// Revision 1.5  2006/02/19 01:37:34  jfrazier
// Initial attempt at bringing into conformance with the coding standard.  Still more work to be done.  Has not been compiled.
//
//

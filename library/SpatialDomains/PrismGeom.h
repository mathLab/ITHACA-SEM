////////////////////////////////////////////////////////////////////////////////
//
//  File:  $Source: /usr/sci/projects/Nektar/cvs/Nektar++/library/SpatialDomains/PrismGeom.h,v $
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
#ifndef NEKTAR_SPATIALDOMAINS_PRISMGEOM_H
#define NEKTAR_SPATIALDOMAINS_PRISMGEOM_H

#include <StdRegions/StdPrismExp.h>
#include <StdRegions/StdRegions.hpp>
#include <SpatialDomains/SpatialDomains.hpp>

#include <SpatialDomains/GeomFactors.h>
#include <SpatialDomains/Geometry3D.h>
#include <SpatialDomains/MeshComponents.h>
#include <SpatialDomains/TriGeom.h>
#include <SpatialDomains/QuadGeom.h>

namespace Nektar
{
    namespace SpatialDomains
    {
        class PrismGeom;
        typedef boost::shared_ptr<PrismGeom> PrismGeomSharedPtr;
        typedef std::vector< PrismGeomSharedPtr > PrismGeomVector;
        typedef std::vector< PrismGeomSharedPtr >::iterator PrismGeomVectorIter;

        class PrismGeom: public LibUtilities::GraphVertexObject, public Geometry3D
        {
        public:
            PrismGeom(const TriGeomSharedPtr faces[], const StdRegions::FaceOrientation forient[]);
            ~PrismGeom();

            void AddElmtConnected(int gvo_id, int locid);
            int  NumElmtConnected() const;
            bool IsElmtConnected(int gvo_id, int locid) const;

            inline int GetEid() const 
            {
                return m_eid;
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

            static const int kNverts = 6;
            static const int kNedges = 9;
            static const int kNqfaces = 3;
            static const int kNtfaces = 2;
            static const int kNfaces = kNqfaces + kNtfaces;

        protected:

            VertexComponentVector           m_verts;
            EdgeComponentVector             m_edges;
            TriFaceComponentVector          m_tfaces;
            QuadFaceComponentVector         m_qfaces;
            StdRegions::EdgeOrientation     m_eorient [kNedges];
            StdRegions::FaceOrientation     m_forient[kNfaces];

            int m_eid;

            bool m_ownverts;
            std::list<CompToElmt> m_elmtmap;

            Array<OneD, StdRegions::StdExpansion3DSharedPtr> m_xmap;

        private:
            PrismGeom ();
        };
    }; //end of namespace
}; //end of namespace

#endif //NEKTAR_SPATIALDOMAINS_PRISMGEOM_H

//
// $Log: PrismGeom.h,v $
// Revision 1.2  2008/02/03 05:05:16  jfrazier
// Initial checkin of 3D components.
//
// Revision 1.1  2006/05/04 18:59:02  kirby
// *** empty log message ***
//
// Revision 1.15  2006/04/09 02:08:35  jfrazier
// Added precompiled header.
//
// Revision 1.14  2006/03/12 14:20:43  sherwin
//
// First compiling version of SpatialDomains and associated modifications
//
// Revision 1.13  2006/03/12 11:06:39  sherwin
//
// First complete copy of code standard code but still not compile tested
//
// Revision 1.12  2006/03/04 20:26:05  bnelson
// Added comments after #endif.
//
// Revision 1.11  2006/02/19 01:37:34  jfrazier
// Initial attempt at bringing into conformance with the coding standard.  Still more work to be done.  Has not been compiled.
//
//


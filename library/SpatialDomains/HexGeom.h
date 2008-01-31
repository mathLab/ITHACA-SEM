////////////////////////////////////////////////////////////////////////////////
//
//  File:  $Source: /usr/sci/projects/Nektar/cvs/Nektar++/library/SpatialDomains/HexGeom.h,v $
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
#ifndef NEKTAR_SPATIALDOMAINS_HEXGEOM_H
#define NEKTAR_SPATIALDOMAINS_HEXGEOM_H

#include <StdRegions/StdHexExp.h>
#include <StdRegions/StdRegions.hpp>
#include <SpatialDomains/SpatialDomains.hpp>

#include <SpatialDomains/MeshGraph.h>
#include <SpatialDomains/GeomFactors.h>
#include <SpatialDomains/Geometry3D.h>
#include <SpatialDomains/MeshComponents.h>
#include <SpatialDomains/EdgeComponent.h>
#include <SpatialDomains/QuadFaceComponent.h>


namespace Nektar
{
    namespace SpatialDomains
    {
        class HexGeom: public LibUtilities::GraphVertexObject, public Geometry3D
        {
        public:
            HexGeom ();
            HexGeom(double *verts);
            ~HexGeom();

        protected:
            static const int kNverts = 8;
            static const int kNedges = 12;
            static const int kNqfaces = 6;
            static const int kNtfaces = 0;
            static const int kNfaces = kNqfaces + kNtfaces;

            VertexComponent*    m_verts [kNverts];
            EdgeComponent*      m_edges [kNedges];
            QuadFaceComponent*  m_qfaces[kNqfaces];

            StdRegions::EdgeOrientation* m_eorient[kNedges];
            StdRegions::FaceOrientation* m_forient[kNfaces];


            StdRegions::StdHexExp** m_xmap;

            int m_eid;
        
        private:

        };

	// shorthand for boost pointer
        typedef boost::shared_ptr<HexGeom> HexGeomSharedPtr;
        typedef std::vector< HexGeomSharedPtr > HexGeomVector;
        typedef std::vector< HexGeomSharedPtr >::iterator HexGeomVectorIter;


    }; //end of namespace
}; //end of namespace

#endif //NEKTAR_SPATIALDOMAINS_HEXGEOM_H

//
// $Log: HexGeom.h,v $
// Revision 1.2  2007/07/20 02:15:08  bnelson
// Replaced boost::shared_ptr with Nektar::ptr
//
// Revision 1.1  2006/05/04 18:59:00  kirby
// *** empty log message ***
//
// Revision 1.16  2006/03/12 14:20:43  sherwin
//
// First compiling version of SpatialDomains and associated modifications
//
// Revision 1.15  2006/03/12 07:42:02  sherwin
//
// Updated member names and StdRegions call. Still has not been compiled
//
// Revision 1.14  2006/03/04 20:26:04  bnelson
// Added comments after #endif.
//
// Revision 1.13  2006/02/19 01:37:33  jfrazier
// Initial attempt at bringing into conformance with the coding standard.  Still more work to be done.  Has not been compiled.
//
//

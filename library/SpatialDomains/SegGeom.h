////////////////////////////////////////////////////////////////////////////////
//
//  File:  $Source: /usr/sci/projects/Nektar/cvs/Nektar++/libs/SpatialDomains/SegGeom.h,v $
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
#ifndef NEKTAR_SPATIALDOMAINS_SEGGEOM_H
#define NEKTAR_SPATIALDOMAINS_SEGGEOM_H

#include <StdRegions/StdRegions.hpp>
#include <SpatialDomains/SpatialDomains.hpp>

#include <SpatialDomains/MeshGraph.h>
#include <SpatialDomains/GeoFac.h>
#include <SpatialDomains/Geometry1D.h>
#include <SpatialDomains/MeshComponents.h>
#include <SpatialDomains/EdgeComponent.h>

namespace Nektar
{
    namespace SpatialDomains
    {
        class SegGeom: public EdgeComponent
        {
        public:
            SegGeom();
            SegGeom(int id, VertexComponentSharedPtr vert1, VertexComponentSharedPtr  vert2);

            SegGeom(SegGeom &in);

            ~SegGeom();

            void    FillGeom ();

            StdRegions::ShapeType DetShapeType() 
            {
                return StdRegions::eSegment;
            }

            GeoFac *GenXGeoFac(void);

            void GetLocCoords(double *Lcoords, const double *coords);
            inline void SetOwnData()
            {
                m_owndata = true; 
            }

            void WriteToFile(std::ofstream &outfile, const int dumpVar);

        protected:
            static const int kNverts = 2;
            SpatialDomains::VertexComponentSharedPtr m_verts[kNverts];

        private:
            bool m_owndata;   ///< Boolean indicating whether object owns the data
            virtual StdRegions::ShapeType V_DetShapeType() 
            {
                return DetShapeType();
            }
        };

	// shorthand for boost pointer
	typedef boost::shared_ptr<SegGeom> SegGeomSharedPtr;
        typedef std::vector< SegGeomSharedPtr > SegGeomVector;
        typedef std::vector< SegGeomSharedPtr >::iterator SegGeomVectorIter;

    }; //end of namespace
}; //end of namespace

#endif //NEKTAR_SPATIALDOMAINS_SEGGEOM_H

//
// $Log: SegGeom.h,v $
// Revision 1.3  2006/05/29 17:05:17  sherwin
// Updated to use shared_ptr around Geom types - added typedef
//
// Revision 1.2  2006/05/09 13:37:01  jfrazier
// Removed duplicate definition of shared vertex pointer.
//
// Revision 1.1  2006/05/04 18:59:04  kirby
// *** empty log message ***
//
// Revision 1.25  2006/04/04 23:12:38  jfrazier
// More updates to readers.  Still work to do on MeshGraph2D to store tris and quads.
//
// Revision 1.24  2006/03/25 00:58:29  jfrazier
// Many changes dealing with fundamental structure and reading/writing.
//
// Revision 1.23  2006/03/12 14:20:43  sherwin
//
// First compiling version of SpatialDomains and associated modifications
//
// Revision 1.22  2006/03/12 11:06:40  sherwin
//
// First complete copy of code standard code but still not compile tested
//
// Revision 1.21  2006/03/04 20:26:05  bnelson
// Added comments after #endif.
//
// Revision 1.20  2006/02/19 01:37:34  jfrazier
// Initial attempt at bringing into conformance with the coding standard.  Still more work to be done.  Has not been compiled.
//
//

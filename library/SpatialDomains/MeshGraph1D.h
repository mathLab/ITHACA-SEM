///////////////////////////////////////////////////////////////////////////////
//
//  File:  $Source: /usr/sci/projects/Nektar/cvs/Nektar++/libs/SpatialDomains/MeshGraph1D.h,v $
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
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_SPATIALDOMAINS_MESHGRAPH1D_H
#define NEKTAR_SPATIALDOMAINS_MESHGRAPH1D_H

#include <SpatialDomains/MeshGraph.h>
#include <SpatialDomains/SegGeom.h>

#include <list>

namespace Nektar
{
    namespace SpatialDomains
    {
        typedef boost::shared_ptr<SegGeom> SharedSegGeomPtr;
        typedef std::vector< SharedSegGeomPtr > SegGeomVector;

        class MeshGraph1D: public MeshGraph
        {
        public:

            MeshGraph1D();
            virtual ~MeshGraph1D();

            void Read(std::string &infilename);
            void Read(TiXmlDocument &doc);

            void Write(std::string &outfilename);

            inline const int GetCoordim(void){
                return GetSpaceDimension();
            }

            inline const SegGeomVector &GetSeggeoms(void) const
            {
                return m_seggeoms;
            }

            inline bool GetGeofac_defined(void)
            {
                return m_geofac_defined;
            }

            void GenXGeoFac();

        protected:
            bool   m_geofac_defined;

        private:
            SegGeomVector m_seggeoms;
        };
    }; //end of namespace
}; //end of namespace

#endif //NEKTAR_SPATIALDOMAINS_MESHGRAPH1D_H

//
// $Log: MeshGraph1D.h,v $
// Revision 1.15  2006/05/02 21:21:11  sherwin
// Corrected libraries to compile new version of spatialdomains and demo Graph1D
//
// Revision 1.14  2006/04/04 23:12:37  jfrazier
// More updates to readers.  Still work to do on MeshGraph2D to store tris and quads.
//
// Revision 1.13  2006/03/25 00:58:29  jfrazier
// Many changes dealing with fundamental structure and reading/writing.
//
// Revision 1.12  2006/03/12 14:20:43  sherwin
//
// First compiling version of SpatialDomains and associated modifications
//
// Revision 1.11  2006/03/12 07:42:03  sherwin
//
// Updated member names and StdRegions call. Still has not been compiled
//
// Revision 1.10  2006/03/04 20:26:04  bnelson
// Added comments after #endif.
//
// Revision 1.9  2006/02/19 01:37:33  jfrazier
// Initial attempt at bringing into conformance with the coding standard.  Still more work to be done.  Has not been compiled.
//
//

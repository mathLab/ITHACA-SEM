///////////////////////////////////////////////////////////////////////////////
//
//  File:  $Source: /usr/sci/projects/Nektar/cvs/Nektar++/library/SpatialDomains/MeshGraph1D.h,v $
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

#include <SpatialDomains/SpatialDomainsDeclspec.h>
namespace Nektar
{
    namespace SpatialDomains
    {
        class MeshGraph1D: public MeshGraph
        {
        public:

            SPATIAL_DOMAINS_EXPORT MeshGraph1D();
            SPATIAL_DOMAINS_EXPORT MeshGraph1D(const LibUtilities::SessionReaderSharedPtr &pSession);
            SPATIAL_DOMAINS_EXPORT virtual ~MeshGraph1D();

            SPATIAL_DOMAINS_EXPORT void ReadGeometry(const std::string &infilename);
            SPATIAL_DOMAINS_EXPORT void ReadGeometry(TiXmlDocument &doc);
            SPATIAL_DOMAINS_EXPORT void ReadElements(TiXmlDocument &doc);
            SPATIAL_DOMAINS_EXPORT void ReadComposites(TiXmlDocument &doc);
            SPATIAL_DOMAINS_EXPORT void ResolveGeomRef(const std::string
                    &prevToken, const std::string &token, Composite& composite);

            inline int GetCoordim(void)
            {
                return GetSpaceDimension();
            }

            inline const SegGeomMap &GetSeggeoms(void) const
            {
                return m_segGeoms;
            }

            inline int GetVidFromElmt(const int vert, const int elmt) const
            {
                ASSERTL2((elmt >=0)&&(elmt < m_segGeoms.size()),
                    "eid is out of range");

                //return m_segGeoms[elmt]->GetVid(vert);
                return m_segGeoms.find(elmt)->second->GetVid(vert);
            }

        protected:

        private:
        };

        typedef boost::shared_ptr<MeshGraph1D> MeshGraph1DSharedPtr;

    }; //end of namespace
}; //end of namespace

#endif //NEKTAR_SPATIALDOMAINS_MESHGRAPH1D_H

//
// $Log: MeshGraph1D.h,v $
// Revision 1.10  2008/06/09 21:33:04  jfrazier
// Moved segment vector to base MeshGraph class since it is used by all derived types.
//
// Revision 1.9  2008/05/29 19:07:39  delisi
// Removed the Write(..) methods, so it is only in the base MeshGraph class. Also, added a line to set the global ID of the geometry object for every element read in.
//
// Revision 1.8  2007/09/20 22:25:06  jfrazier
// Added expansion information to meshgraph class.
//
// Revision 1.7  2007/06/07 23:55:24  jfrazier
// Intermediate revisions to add parsing for boundary conditions file.
//
// Revision 1.6  2007/05/28 16:15:01  sherwin
// Updated files in MultiRegions to make 1D demos work
//
// Revision 1.5  2007/03/14 21:24:08  sherwin
// Update for working version of MultiRegions up to ExpList1D
//
// Revision 1.4  2006/07/02 17:16:17  sherwin
//
// Modifications to make MultiRegions work for a connected domain in 2D (Tris)
//
// Revision 1.3  2006/06/01 14:15:30  sherwin
// Added typdef of boost wrappers and made GeoFac a boost shared pointer.
//
// Revision 1.2  2006/05/30 14:00:04  sherwin
// Updates to make MultiRegions and its Demos work
//
// Revision 1.1  2006/05/04 18:59:01  kirby
// *** empty log message ***
//
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

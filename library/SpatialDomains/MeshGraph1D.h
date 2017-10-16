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
            SPATIAL_DOMAINS_EXPORT MeshGraph1D(const LibUtilities::SessionReaderSharedPtr &pSession, const DomainRangeShPtr &rng = NullDomainRangeShPtr);
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

        typedef std::shared_ptr<MeshGraph1D> MeshGraph1DSharedPtr;

    }; //end of namespace
}; //end of namespace

#endif //NEKTAR_SPATIALDOMAINS_MESHGRAPH1D_H

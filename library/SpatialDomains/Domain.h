////////////////////////////////////////////////////////////////////////////////
//
//  File: Domain.h
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
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
////////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_SPATIALDOMAINS_DOMAIN_H
#define NEKTAR_SPATIALDOMAINS_DOMAIN_H

#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <SpatialDomains/MeshGraph.h>

class TiXmlDocument;

namespace Nektar
{
    namespace SpatialDomains
    {
        enum BoundaryType
        {
            eUnknown=0,
            eDirichlet,
            eNeumann,
            eRobin,
            eCauchy,

            eDummy,
            eBoundaryTypeLastElement = eDummy-1
        };

        // Corresponds to the entries above.  These are the tags within the domain definition
        // corresponding to type of BC above.
        const char BoundaryTypeNameMap[] =
        {
            'U',// Just a placeholder to get the correct index for the rest.
            'D',
            'N',
            'R',
            'C'
        };

        struct BoundaryEntry
        {
            BoundaryType m_BoundaryType;
            std::vector< Composite > m_BoundaryComposites;
        };

        typedef std::shared_ptr< BoundaryEntry > BoundarySharedPtr;
        typedef std::vector< BoundarySharedPtr > BoundaryVector;
        typedef std::vector< Composite > CompositeVector;

        class Domain
        {
        public:
            // Must have a MeshGraph from which the composites
            // and associated items can be obtained.
            Domain(MeshGraph *meshGraph);
            virtual ~Domain();

            void Read(std::string &infilename);
            void Read(TiXmlDocument &doc);
            void Write(std::string &outfilename);

            CompositeVector GetDomain(void) const { return m_Domain; };
            BoundaryVector  GetBoundaries(void) const { return m_Boundaries; };

            inline void SetFileName(const std::string &inString)
            {
                m_FileName = inString;
            };


        protected:
            //std::string     m_filename;
            MeshGraph      *m_meshGraph;
            CompositeVector m_domain;
            BoundaryVector  m_boundaries;

        private:
            Domain(){ NEKERROR(ErrorUtil::efatal, "Must provide a meshgraph to create a Domain."); };   // Don't call this.
        };
    };
}

#endif //NEKTAR_SPATIALDOMAINS_DOMAIN_H

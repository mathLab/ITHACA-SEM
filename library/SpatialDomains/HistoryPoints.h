////////////////////////////////////////////////////////////////////////////////
//
//  File:  History.h
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
//  Description: Stores history points.
//
//
////////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_SPATIALDOMAINS_HISTORY_H
#define NEKTAR_SPATIALDOMAINS_HISTORY_H

#include <string>
#include <iostream>

#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
#include <SpatialDomains/MeshGraph.h>
#include <SpatialDomains/SpatialDomainsDeclspec.h>
class TiXmlElement;
class TiXmlDocument;

class MeshGraph;

namespace Nektar
{
    namespace SpatialDomains
    {
        class History
        {
            public:
                SPATIAL_DOMAINS_EXPORT History(const MeshGraph *meshGraph);
                SPATIAL_DOMAINS_EXPORT void Read(std::string &infilename);
                SPATIAL_DOMAINS_EXPORT void Read(TiXmlDocument &doc);

                SPATIAL_DOMAINS_EXPORT int GetNumHistoryPoints() const;
                SPATIAL_DOMAINS_EXPORT VertexComponentSharedPtr GetHistoryPoint(int idx) const;

            protected:
                void ReadHistoryPoints(TiXmlElement *history);

            private:
                VertexComponentVector       m_historyPoints;

                /// The mesh graph to use for referencing geometry info.
                const MeshGraph *m_meshGraph;
        };

        typedef boost::shared_ptr<History> HistorySharedPtr;
    }
}
#endif

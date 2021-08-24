////////////////////////////////////////////////////////////////////////////////
//
//  File:  MeshComponents.h
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
//
////////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_SPATIALDOMAINS_MESHCOMPONENTS_H
#define NEKTAR_SPATIALDOMAINS_MESHCOMPONENTS_H

#include <SpatialDomains/Geometry0D.h>
#include <LibUtilities/LinearAlgebra/NekPoint.hpp>
#include <SpatialDomains/SpatialDomainsDeclspec.h>
#include <set>
#include <LibUtilities/Foundations/Graph.h>

namespace Nektar
{
      namespace SpatialDomains
      {
        // ------------------------------------------------------------------------
        /// Structure holding graphvertexobject id and local element facet id
        class CompToElmt
        {
        public:

            CompToElmt(int id, int locid):
              m_id(id),
                  m_locId(locid)
              {
                  m_id = id;
                  m_locId = locid;
              }

              ~CompToElmt()
              {
                  m_id = -1;
                  m_locId = -1;
              }

              inline int GetId()
              {
                  return m_id;
              }

              SPATIAL_DOMAINS_EXPORT friend bool operator  == (const CompToElmt &x, const CompToElmt &y);
              SPATIAL_DOMAINS_EXPORT friend bool operator  != (const CompToElmt &x, const CompToElmt &y);

        protected:
            int m_id;
            int m_locId;

        private:

        };


        // -----------------------------------------------------------------------
        // WireFrame

        class WireframeEdgeComponent: public LibUtilities::GraphEdgeObject
        {
        public:
            WireframeEdgeComponent(int gvoid1, int gvoid2)
            {
                m_gvoid1 = gvoid1;
                m_gvoid2 = gvoid2;
            }

            ~WireframeEdgeComponent()
            {
            }

            void GetConnectivity(int &gvoid1, int &gvoid2) const
            {
                gvoid1 = m_gvoid1;
                gvoid2 = m_gvoid2;
            }

        protected:
        private:
        };

    } //end of namespace
} //end of namespace

#endif //NEKTAR_SPATIALDOMAINS_MESHCOMPONENTS_H


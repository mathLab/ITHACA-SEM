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
#ifndef NEKTAR_SPATIALDOMAINS_MESHCOMPONENTS_H
#define NEKTAR_SPATIALDOMAINS_MESHCOMPONENTS_H

#include <SpatialDomains/Geometry.h>
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

        // --------------------------------------------------------------------
        /// Vertex Component
        class VertexComponent: public Geometry, public NekPoint <NekDouble>
        {
        public:
                SPATIAL_DOMAINS_EXPORT VertexComponent(const int coordim, const int vid,
                    NekDouble x, NekDouble y, NekDouble z);
                SPATIAL_DOMAINS_EXPORT VertexComponent(){}
                SPATIAL_DOMAINS_EXPORT ~VertexComponent();
                SPATIAL_DOMAINS_EXPORT VertexComponent(const VertexComponent &T);

                SPATIAL_DOMAINS_EXPORT void AddElmtConnected(int gvo_id, int locid);
                SPATIAL_DOMAINS_EXPORT int  NumElmtConnected() const;
                SPATIAL_DOMAINS_EXPORT bool IsElmtConnected(int gvo_id, int locid) const;
                SPATIAL_DOMAINS_EXPORT void GetCoords(NekDouble &x, NekDouble &y, NekDouble &z);
                SPATIAL_DOMAINS_EXPORT void GetCoords(Array<OneD,NekDouble> &coords);
                SPATIAL_DOMAINS_EXPORT void UpdatePosition(NekDouble x, NekDouble y, NekDouble z);


                inline int GetCoordim() const
                {
                    return m_coordim;
                }

                inline int GetVid() const
                {
                    return m_vid;
                }

                inline void SetVid(const int vid)
                {
                    m_vid = vid;
                }

                // Math routines
                SPATIAL_DOMAINS_EXPORT void   Mult (VertexComponent& a, VertexComponent& b);
                SPATIAL_DOMAINS_EXPORT void   Add  (VertexComponent& a, VertexComponent& b);
                SPATIAL_DOMAINS_EXPORT void   Sub  (VertexComponent& a, VertexComponent& b);
                SPATIAL_DOMAINS_EXPORT NekDouble dist  (VertexComponent& a);
                SPATIAL_DOMAINS_EXPORT NekDouble dot   (VertexComponent& a);

                SPATIAL_DOMAINS_EXPORT friend bool operator == (const VertexComponent &x, const VertexComponent &y);
                SPATIAL_DOMAINS_EXPORT friend bool operator == (const VertexComponent &x, const VertexComponent *y);
                SPATIAL_DOMAINS_EXPORT friend bool operator == (const VertexComponent *x, const VertexComponent &y);
                SPATIAL_DOMAINS_EXPORT friend bool operator != (const VertexComponent &x, const VertexComponent &y);
                SPATIAL_DOMAINS_EXPORT friend bool operator != (const VertexComponent &x, const VertexComponent *y);
                SPATIAL_DOMAINS_EXPORT friend bool operator != (const VertexComponent *x, const VertexComponent &y);

            protected:
                int m_vid;
                int m_coordim;
                std::list<CompToElmt> m_elmtMap;
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

        typedef boost::shared_ptr< VertexComponent >  VertexComponentSharedPtr;
        typedef std::set< VertexComponentSharedPtr >  VertexComponentSet;
        typedef std::vector< VertexComponentSharedPtr >  VertexComponentVector;
        typedef std::vector< VertexComponentSharedPtr >::iterator  VertexComponentVectorIter;

    }; //end of namespace
}; //end of namespace

#endif //NEKTAR_SPATIALDOMAINS_MESHCOMPONENTS_H


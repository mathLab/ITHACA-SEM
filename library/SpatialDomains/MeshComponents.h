////////////////////////////////////////////////////////////////////////////////
//
//  File:  $Source: /usr/sci/projects/Nektar/cvs/Nektar++/libs/SpatialDomains/MeshComponents.h,v $
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

#include <SpatialDomains/SpatialDomains.hpp>
#include <LibUtilities/NekPoint.hpp>
#include <LibUtilities/Vmath.hpp>
#include <LibUtilities/Blas.hpp>
#include <LibUtilities/NekMemoryManager.hpp>
#include <set>

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
                  m_locid(locid)
              {
                  m_id = id;
                  m_locid = locid;
              }

              ~CompToElmt()
              {
                  m_id = -1; 
                  m_locid = -1;
              }

              inline int GetId()
              {
                  return m_id;
              }

              friend bool operator  == (const CompToElmt &x, const CompToElmt &y);
              friend bool operator  != (const CompToElmt &x, const CompToElmt &y);

        protected:
            int m_id;
            int m_locid;

        private:

        };

        // -----------------------------------------------------------------------
        /// Vertex Component
        class VertexComponent: public LibUtilities::NekPoint<double, 3>
        {
        public:
            //Temp debug constructor
            VertexComponent()
            {
            }

            VertexComponent::VertexComponent(const int coordim, const int vid, 
                double x, double y, double z);
            ~VertexComponent();
            VertexComponent(const VertexComponent &T);

            void AddElmtConnected(int gvo_id, int locid);
            int  NumElmtConnected() const;
            bool IsElmtConnected(int gvo_id, int locid) const;
            void UpdatePosition(double x, double y, double z);

            inline int  GetCoordim()
            {
                return(m_coordim);
            }

            inline int  GetVid()
            {
                return m_vid ;
            }

            inline void SetVid(const int vid)
            {
                m_vid = vid;
            }

            // Math routines
            void   Mult (VertexComponent& a, VertexComponent& b);
            void   Add  (VertexComponent& a, VertexComponent& b);
            void   Sub  (VertexComponent& a, VertexComponent& b);
            double dot  (VertexComponent& a);

            friend bool operator == (const VertexComponent &x, const VertexComponent &y);
            friend bool operator == (const VertexComponent &x, const VertexComponent *y);
            friend bool operator == (const VertexComponent *x, const VertexComponent &y);
            friend bool operator != (const VertexComponent &x, const VertexComponent &y);
            friend bool operator != (const VertexComponent &x, const VertexComponent *y);
            friend bool operator != (const VertexComponent *x, const VertexComponent &y);

        protected:
            int m_vid;
            int m_coordim;
            std::list<CompToElmt> m_elmtmap;

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

        typedef boost::shared_ptr< VertexComponent >  VertexComponentSharedPtr;
        typedef std::set< VertexComponentSharedPtr >  VertexComponentSet;
        typedef std::vector< VertexComponentSharedPtr >  VertexComponentVector;
        typedef std::vector< VertexComponentSharedPtr >::iterator  VertexComponentVectorIter;

    }; //end of namespace
}; //end of namespace

#endif //NEKTAR_SPATIALDOMAINS_MESHCOMPONENTS_H

//
// $Log: MeshComponents.h,v $
// Revision 1.31  2006/05/02 21:21:11  sherwin
// Corrected libraries to compile new version of spatialdomains and demo Graph1D
//
// Revision 1.30  2006/04/11 23:18:11  jfrazier
// Completed MeshGraph2D for tri's and quads.  Not thoroughly tested.
//
// Revision 1.29  2006/03/25 00:58:29  jfrazier
// Many changes dealing with fundamental structure and reading/writing.
//
// Revision 1.28  2006/03/12 14:20:43  sherwin
//
// First compiling version of SpatialDomains and associated modifications
//
// Revision 1.27  2006/03/12 07:42:02  sherwin
//
// Updated member names and StdRegions call. Still has not been compiled
//
// Revision 1.26  2006/03/04 20:26:04  bnelson
// Added comments after #endif.
//
// Revision 1.25  2006/02/26 21:19:43  bnelson
// Fixed a variety of compiler errors caused by updates to the coding standard.
//
// Revision 1.24  2006/02/19 01:37:33  jfrazier
// Initial attempt at bringing into conformance with the coding standard.  Still more work to be done.  Has not been compiled.
//
//


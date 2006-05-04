////////////////////////////////////////////////////////////////////////////////
//
//  File:  $Source: /usr/sci/projects/Nektar/cvs/Nektar++/libs/SpatialDomains/MeshComponents.cpp,v $
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
#include "pchSpatialDomains.h"

#include <SpatialDomains/MeshComponents.h>

namespace Nektar
{
    namespace SpatialDomains
    {
        VertexComponent::VertexComponent(const int coordim, const int vid, 
            double x, double y, double z)
        {
            m_coordim = coordim;
            m_vid     = vid;

            (*this)(0) = x;
            (*this)(1) = y;
            (*this)(2) = z;
        }

        VertexComponent::~VertexComponent(){
        }

        // copy constructor
        VertexComponent::VertexComponent(const VertexComponent &T)
            :LibUtilities::NekPoint<double,3>(T)
        {
            m_vid = T.m_vid;
            m_coordim = T.m_coordim;

            std::list<CompToElmt>::const_iterator def;
            for(def = T.m_elmtmap.begin(); def != T.m_elmtmap.end(); def++)
            {
                m_elmtmap.push_back(*def);    
            }
        }

        void VertexComponent::AddElmtConnected(int gvo_id, int locid)
        {
            CompToElmt ee(gvo_id,locid);
            m_elmtmap.push_back(ee);
        }

        int VertexComponent::NumElmtConnected() const
        {
            return int(m_elmtmap.size());
        }

        bool VertexComponent::IsElmtConnected(int gvo_id, int locid) const
        {

            std::list<CompToElmt>::const_iterator def;
            CompToElmt ee(gvo_id,locid);

            def = find(m_elmtmap.begin(),m_elmtmap.end(),ee);

            // Found the element connectivity object in the list
            if(def != m_elmtmap.end())
            {
                return(true);
            }
            return(false);
        }

        void VertexComponent::UpdatePosition(double x, double y, double z)
        {
            (*this)(0) = x;
            (*this)(1) = x;
            (*this)(2) = x;
        }  

        // _this = a x b 
        void VertexComponent::Mult(VertexComponent& a,VertexComponent& b)
        {
            (*this)(0) = a[1]*b[2] - a[2]*b[1];
            (*this)(1) = a[2]*b[0] - a[0]*b[2];
            (*this)(2) = a[0]*b[1] - a[1]*b[0];
            m_coordim = 3; 
        }

        // _output = this.a 
        double VertexComponent::dot(VertexComponent& a)
        {
            return (x()*a.x() + y()*a.y() + z()*a.z());
        }

        /// Determine equivalence by the ids.  No matter what the position,
        /// if the ids are the same, then they are equivalent, and vice versa.
        bool operator  == (const VertexComponent &x, const VertexComponent &y)
        {
            return (x.m_vid == y.m_vid);
        }

        bool operator  == (const VertexComponent &x, const VertexComponent *y)
        {
            return (x.m_vid == y->m_vid);
        }

        bool operator  == (const VertexComponent *x, const VertexComponent &y)
        {
            return (x->m_vid == y.m_vid);
        }

        bool operator != (const VertexComponent &x, const VertexComponent &y)
        {
            return (x.m_vid != y.m_vid);
        }

        bool operator  != (const VertexComponent &x, const VertexComponent *y)
        {
            return (x.m_vid != y->m_vid);
        }

        bool operator  != (const VertexComponent *x, const VertexComponent &y)
        {
            return (x->m_vid != y.m_vid);
        }

        bool operator  == (const CompToElmt &x, const CompToElmt &y)
        {
            return (x.m_id == y.m_id) || (x.m_locid == y.m_locid);
        }

        bool operator  != (const CompToElmt &x, const CompToElmt &y)
        {
            return (x.m_id != y.m_id);
        }
    }; //end of namespace
}; //end of namespace

//
// $Log: MeshComponents.cpp,v $
// Revision 1.19  2006/04/09 02:08:35  jfrazier
// Added precompiled header.
//
// Revision 1.18  2006/03/25 00:58:28  jfrazier
// Many changes dealing with fundamental structure and reading/writing.
//
// Revision 1.17  2006/03/12 14:20:43  sherwin
//
// First compiling version of SpatialDomains and associated modifications
//
// Revision 1.16  2006/02/26 21:19:42  bnelson
// Fixed a variety of compiler errors caused by updates to the coding standard.
//
// Revision 1.15  2006/02/19 01:37:33  jfrazier
// Initial attempt at bringing into conformance with the coding standard.  Still more work to be done.  Has not been compiled.
//
//

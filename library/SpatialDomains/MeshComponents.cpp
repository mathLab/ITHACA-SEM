////////////////////////////////////////////////////////////////////////////////
//
//  File:  $Source: /usr/sci/projects/Nektar/cvs/Nektar++/library/SpatialDomains/MeshComponents.cpp,v $
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
        VertexComponent::VertexComponent(const VertexComponent &T): NekPoint<NekDouble>(T)
        {
            m_vid = T.m_vid;
            m_coordim = T.m_coordim;

            std::list<CompToElmt>::const_iterator def;
            for(def = T.m_elmtMap.begin(); def != T.m_elmtMap.end(); def++)
            {
                m_elmtMap.push_back(*def);
            }
        }

        void VertexComponent::AddElmtConnected(int gvo_id, int locid)
        {
            CompToElmt ee(gvo_id,locid);
            m_elmtMap.push_back(ee);
        }

        int VertexComponent::NumElmtConnected() const
        {
            return int(m_elmtMap.size());
        }

        bool VertexComponent::IsElmtConnected(int gvo_id, int locid) const
        {

            std::list<CompToElmt>::const_iterator def;
            CompToElmt ee(gvo_id,locid);

            def = find(m_elmtMap.begin(),m_elmtMap.end(),ee);

            // Found the element connectivity object in the list
            if(def != m_elmtMap.end())
            {
                return(true);
            }
            return(false);
        }

        void VertexComponent::GetCoords(double &x, double &y, double &z)
        {
            switch(m_coordim)
            {
            case 3:
                z = (*this)(2);
            case 2:
                y = (*this)(1);
            case 1:
                x = (*this)(0);
                break;
            }
        }

        void VertexComponent::GetCoords(Array<OneD,NekDouble> &coords)
        {
            switch(m_coordim)
            {
            case 3:
                coords[2] = (*this)(2);
            case 2:
                coords[1] = (*this)(1);
            case 1:
                coords[0] = (*this)(0);
                break;
            }
        }


        void VertexComponent::UpdatePosition(double x, double y, double z)
        {
            (*this)(0) = x;
            (*this)(1) = y;
            (*this)(2) = z;
        }

        // _this = a + b
        void VertexComponent::Add(VertexComponent& a,VertexComponent& b)
        {
            (*this)(0) = a[0] + b[0];
            (*this)(1) = a[1] + b[1];
            (*this)(2) = a[2] + b[2];
            m_coordim = std::max(a.GetCoordim(),b.GetCoordim());
        }

        // _this = a + b
        void VertexComponent::Sub(VertexComponent& a,VertexComponent& b)
        {
            (*this)(0) = a[0] - b[0];
            (*this)(1) = a[1] - b[1];
            (*this)(2) = a[2] - b[2];
            m_coordim = std::max(a.GetCoordim(),b.GetCoordim());
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
        NekDouble VertexComponent::dist(VertexComponent& a)
        {
            return sqrt((x()-a.x())*(x()-a.x()) + (y()-a.y())*(y()-a.y()) + (z()-a.z())*(z()-a.z()));
        }

        // _output = this.a
        NekDouble VertexComponent::dot(VertexComponent& a)
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
            return (x.m_id == y.m_id) || (x.m_locId == y.m_locId);
        }

        bool operator  != (const CompToElmt &x, const CompToElmt &y)
        {
            return (x.m_id != y.m_id);
        }
    }; //end of namespace
}; //end of namespace

//
// $Log: MeshComponents.cpp,v $
// Revision 1.6  2008/03/06 05:12:20  ehan
// Changed NekVector<NekDouble, 3> to NekVector<NekDouble, ThreeD> to match the changes of enum to struct typename dim.
//
// Revision 1.5  2007/08/02 12:03:27  sherwin
// Fixed order of coordinate assignment in VertexComponent::GetCoords
//
// Revision 1.4  2006/10/15 06:18:57  sherwin
// Moved NekPoint out of namespace LibUtilities
//
// Revision 1.3  2006/05/09 13:37:01  jfrazier
// Removed duplicate definition of shared vertex pointer.
//
// Revision 1.2  2006/05/07 11:26:38  sherwin
// Modifications to get the demo LocalRegions::Project2D to compile
//
// Revision 1.1  2006/05/04 18:59:01  kirby
// *** empty log message ***
//
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

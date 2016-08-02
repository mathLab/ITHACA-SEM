////////////////////////////////////////////////////////////////////////////////
//
//  File:  PointGeom.cpp
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
//  Description: Point geometry information
//
//
////////////////////////////////////////////////////////////////////////////////
//#include "pchSpatialDomains.h"

#include <SpatialDomains/PointGeom.h>
#include <SpatialDomains/SegGeom.h>


#include <StdRegions/StdRegions.hpp>


#include <fstream>

namespace Nektar
{
    namespace SpatialDomains
    {
        PointGeom::PointGeom()
                : NekPoint<NekDouble>(0.0, 0.0, 0.0)
        {
            m_shapeType = LibUtilities::ePoint;
            m_coordim = 0;
            m_vid = 0;
        }

        PointGeom::PointGeom(const int coordim, const int vid,
            NekDouble x, NekDouble y, NekDouble z)
                : NekPoint<NekDouble>(x,y,z)
        {
            m_shapeType = LibUtilities::ePoint;
            m_coordim = coordim;
            m_vid     = vid;
            m_globalID = vid;

            (*this)(0) = x;
            (*this)(1) = y;
            (*this)(2) = z;
        }


        // copy constructor
        PointGeom::PointGeom(const PointGeom &T): NekPoint<NekDouble>(T)
        {
            m_shapeType = T.m_shapeType;
            m_vid = T.m_vid;
            m_coordim = T.m_coordim;
            m_globalID = T.m_globalID;

            std::list<CompToElmt>::const_iterator def;
            for(def = T.m_elmtMap.begin(); def != T.m_elmtMap.end(); def++)
            {
                m_elmtMap.push_back(*def);
            }
        }


        PointGeom::~PointGeom()
        {
        }


        void PointGeom::AddElmtConnected(int gvo_id, int locid)
        {
            CompToElmt ee(gvo_id,locid);
            m_elmtMap.push_back(ee);
        }

        int PointGeom::NumElmtConnected() const
        {
            return int(m_elmtMap.size());
        }

        bool PointGeom::IsElmtConnected(int gvo_id, int locid) const
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

        void PointGeom::GetCoords(NekDouble &x, NekDouble &y, NekDouble &z)
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

        void PointGeom::GetCoords(Array<OneD,NekDouble> &coords)
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


        void PointGeom::UpdatePosition(NekDouble x, NekDouble y, NekDouble z)
        {
            (*this)(0) = x;
            (*this)(1) = y;
            (*this)(2) = z;
        }

        // _this = a + b
        void PointGeom::Add(PointGeom& a,PointGeom& b)
        {
            (*this)(0) = a[0] + b[0];
            (*this)(1) = a[1] + b[1];
            (*this)(2) = a[2] + b[2];
            m_coordim = std::max(a.GetCoordim(),b.GetCoordim());
        }

        // _this = a + b
        void PointGeom::Sub(PointGeom& a,PointGeom& b)
        {
            (*this)(0) = a[0] - b[0];
            (*this)(1) = a[1] - b[1];
            (*this)(2) = a[2] - b[2];
            m_coordim = std::max(a.GetCoordim(),b.GetCoordim());
        }

        // _this = a x b
        void PointGeom::Mult(PointGeom& a,PointGeom& b)
        {
            (*this)(0) = a[1]*b[2] - a[2]*b[1];
            (*this)(1) = a[2]*b[0] - a[0]*b[2];
            (*this)(2) = a[0]*b[1] - a[1]*b[0];
            m_coordim = 3;
        }

        // _output = this.a
        NekDouble PointGeom::dist(PointGeom& a)
        {
            return sqrt((x()-a.x())*(x()-a.x()) + (y()-a.y())*(y()-a.y()) + (z()-a.z())*(z()-a.z()));
        }

        // _output = this.a
        NekDouble PointGeom::dot(PointGeom& a)
        {
            return (x()*a.x() + y()*a.y() + z()*a.z());
        }

        /// Determine equivalence by the ids.  No matter what the position,
        /// if the ids are the same, then they are equivalent, and vice versa.
        bool operator  == (const PointGeom &x, const PointGeom &y)
        {
            return (x.m_vid == y.m_vid);
        }

        bool operator  == (const PointGeom &x, const PointGeom *y)
        {
            return (x.m_vid == y->m_vid);
        }

        bool operator  == (const PointGeom *x, const PointGeom &y)
        {
            return (x->m_vid == y.m_vid);
        }

        bool operator != (const PointGeom &x, const PointGeom &y)
        {
            return (x.m_vid != y.m_vid);
        }

        bool operator  != (const PointGeom &x, const PointGeom *y)
        {
            return (x.m_vid != y->m_vid);
        }

        bool operator  != (const PointGeom *x, const PointGeom &y)
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

        int PointGeom::v_GetVid(int id) const
        {
            return m_vid;
        }

        PointGeomSharedPtr PointGeom::v_GetVertex(int i) const
        {
            ASSERTL0(i == 0, "Index other than 0 is meaningless.");
            // shared_this_ptr() returns const PointGeom, which cannot be
            // returned.
            return PointGeomSharedPtr(new PointGeom(*this));
        }

        /// \brief Get the orientation of point1; to be used later 
        /// for normal convention
        ///
        /// If edge1 is connected to edge2 in the same direction as
        /// the points comprising edge1 then it is forward, otherwise
        /// it is backward.
        ///
        /// For example, assume edge1 is comprised of points 1 and 2,
        /// and edge2 is comprised of points 2 and 3, then edge1 is
        /// forward.
        ///
        /// If edge1 is comprised of points 2 and 1 and edge2 is
        /// comprised of points 3 and 2, then edge1 is backward.
        ///
        /// Since both edges are passed, it does
        /// not need any information from the EdgeComponent instance.

        StdRegions::Orientation PointGeom::GetPointOrientation(const SegGeom &edge1,  const SegGeom &edge2)
        {
            StdRegions::Orientation returnval = StdRegions::eFwd;

            /// Backward direction.  Vertex 0 is connected to edge 2.
            if ((*edge1.GetVertex(0) == *edge2.GetVertex(0)) ||
                (*edge1.GetVertex(0) == *edge2.GetVertex(1)))
            {
                returnval = StdRegions::eBwd;
            }

            // Not forward either, then we have a problem.
            else if ((*edge1.GetVertex(1) != *edge2.GetVertex(0)) &&
                (*edge1.GetVertex(1) != *edge2.GetVertex(1)))
            {
                std::ostringstream errstrm;
                errstrm << "Connected edges do not share a vertex. Edges ";
                errstrm << edge1.GetEid() << ", " << edge2.GetEid();
                ASSERTL0(false, errstrm.str());
            }

            return returnval;
        }

        void PointGeom::v_GenGeomFactors()
        {

        }

        NekDouble PointGeom::v_GetCoord(const int i, const Array<OneD,const NekDouble> &Lcoord)
        {
            return GetCoord(i,Lcoord);
        }

        NekDouble PointGeom::v_GetLocCoords(const Array<OneD,const NekDouble> &coords, Array<OneD,NekDouble> &Lcoords)
        {
            return GetLocCoords(coords,Lcoords);
        }

    }; //end of namespace
}; //end of namespace

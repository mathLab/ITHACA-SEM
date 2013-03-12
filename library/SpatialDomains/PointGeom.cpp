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
        {
            m_shapeType = LibUtilities::ePoint;
        }


        PointGeom::~PointGeom()
        {
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

        VertexComponentSharedPtr PointGeom::v_GetVertex(const int i) const
        {
            VertexComponentSharedPtr returnval;

            if (i >= 0 && i < kNverts)
            {
                returnval = m_verts[i];
            }

            return returnval;
        }


        int PointGeom::v_GetVid(int i) const
        {
            ASSERTL2((i ==0),"Verted id must be 0");
            return m_verts[i]->GetVid();
        }


        NekDouble PointGeom::v_GetCoord(const int i, const Array<OneD,const NekDouble> &Lcoord)
        {
            return GetCoord(i,Lcoord);
        }

        void PointGeom::v_GetLocCoords(const Array<OneD,const NekDouble> &coords, Array<OneD,NekDouble> &Lcoords)
        {
            GetLocCoords(coords,Lcoords);
        }

        int PointGeom::v_GetNumVerts() const
        {
            return kNverts;
        }



    }; //end of namespace
}; //end of namespace

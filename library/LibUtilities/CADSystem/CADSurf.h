////////////////////////////////////////////////////////////////////////////////
//
//  File: CADSurf.h
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
//  Description: CAD object surface.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_CADSYSTEM_CADSURF_H
#define NEKTAR_LIB_UTILITIES_CADSYSTEM_CADSURF_H

#include <boost/shared_ptr.hpp>

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/LibUtilitiesDeclspec.h>
#include <LibUtilities/Memory/NekMemoryManager.hpp>

#include <LibUtilities/CADSystem/OpenCascade.h>

namespace Nektar {
namespace LibUtilities {

class CADSurf
{
    public:
        friend class MemoryManager<CADSurf>;

        /**
         * @brief Default constructor.
         */
        CADSurf(int i, TopoDS_Shape in,
                std::vector<std::vector<std::pair<int,int> > > ein);

        /**
         * @brief Get the IDs of the edges which bound the surface.
         *
         * The edges are organsised into two vectors, which are grouped into the
         * continuous loops of the bounding edges, then the edges, which are a
         * pair of integers. The first item is the edge ID and the second is an
         * integer that indicates whether this edge is orientated forwards or
         * backwards on this surface to form the loop.
         */
        std::vector<std::vector<std::pair<int,int> > > GetEdges()
        {
            return m_edges;
        }

        /**
         * @brief Get the limits of the parametric space for the surface.
         *
         * @return Array of 4 entries with parametric umin,umax,vmin,vmax.
         */
        Array<OneD, NekDouble> GetBounds()
        {
            Array<OneD,NekDouble> b(4);

            b[0] = m_occSurface.FirstUParameter();
            b[1] = m_occSurface.LastUParameter();
            b[2] = m_occSurface.FirstVParameter();
            b[3] = m_occSurface.LastVParameter();

            return b;
        }

        Array<OneD, NekDouble> N    (Array<OneD, NekDouble> uv);
        Array<OneD, NekDouble> D1   (Array<OneD, NekDouble> uv);
        Array<OneD, NekDouble> D2   (Array<OneD, NekDouble> uv);
        Array<OneD, NekDouble> P    (Array<OneD, NekDouble> uv);
        Array<OneD, NekDouble> locuv(Array<OneD, NekDouble> p);

    private:
        /// ID of surface.
        int m_ID;
        /// OpenCascade object for surface.
        BRepAdaptor_Surface m_occSurface;
        /// Alternate OpenCascade object for surface. Used by reverse lookup.
        Handle(Geom_Surface) m_s;
        /// List of bounding edges in loops with orientation.
        std::vector<std::vector<std::pair<int,int> > > m_edges;
};

typedef boost::shared_ptr<CADSurf> CADSurfSharedPtr;

}
}

#endif

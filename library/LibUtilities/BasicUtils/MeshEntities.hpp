///////////////////////////////////////////////////////////////////////////////
//
//  File: MeshEntities.hpp
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
//  Description: Mesh entities that are used in MeshPartition and
//  other I/O routiens in SpatialDomains
//
////////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_LIBUTILITIES_BASICUTILS_MESHENTITIES_HPP
#define NEKTAR_LIBUTILITIES_BASICUTILS_MESHENTITIES_HPP

#include <LibUtilities/Foundations/PointsType.h>
namespace Nektar
{
    namespace LibUtilities
    {
        
        struct MeshVertex
        {
            int id;
            NekDouble x;
            NekDouble y;
            NekDouble z;
        };
        
        struct MeshEdge
        {
            NekInt id;
            NekInt v0;
            NekInt v1;
        };

        struct MeshTri
        {
            NekInt id;
            NekInt e[3];
        };

        struct MeshQuad
        {
            NekInt id;
            NekInt e[4];
        };

        struct MeshTet
        {
            NekInt id;
            NekInt f[4];
        };

        struct MeshPyr
        {
            NekInt id;
            NekInt f[5];
        };

        struct MeshPrism
        {
            NekInt id;
            NekInt f[5];
        };

        struct MeshHex
        {
            NekInt id;
            NekInt f[6];
        };
        
        struct MeshCurvedInfo
        {
            NekInt        id;
            NekInt        entityid;
            NekInt        npoints;
            PointsType    ptype;     
            NekInt        ptid;     // id of point data map
            NekInt        ptoffset; // pffset of data entry for this curve
        };
        
        struct MeshCurvedPts
        {
            NekInt id;
            std::vector<NekInt>     index;
            std::vector<MeshVertex> pts; 
        };
            

    }
}
#endif

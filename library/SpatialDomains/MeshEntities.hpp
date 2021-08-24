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
    namespace SpatialDomains
    {
        // Note: the following structs are defined using 64 bit ints so
        // that the structs are memory aligned in both 64 bit and 32
        // bit machines. All Structs (with the exception of
        // MeshCurvedPts) are therefore aligned to 8 byte blocks.
        //
        // Note the MeshCurvedPts are exported as a list of NekInt64
        // and MeshVetexs and so the struct does not comply with the
        // above.


        struct MeshVertex
        {
            NekInt64 id;
            NekDouble x;
            NekDouble y;
            NekDouble z;
        };

        struct MeshEdge
        {
            NekInt64 id;
            NekInt64 v0;
            NekInt64 v1;
        };

        struct MeshTri
        {
            NekInt64 id;
            NekInt64 e[3];
        };

        struct MeshQuad
        {
            NekInt64 id;
            NekInt64 e[4];
        };

        struct MeshTet
        {
            NekInt64 id;
            NekInt64 f[4];
        };

        struct MeshPyr
        {
            NekInt64 id;
            NekInt64 f[5];
        };

        struct MeshPrism
        {
            NekInt64 id;
            NekInt64 f[5];
        };

        struct MeshHex
        {
            NekInt64 id;
            NekInt64 f[6];
        };

        struct MeshCurvedInfo
        {
            NekInt64        id;        /// Id of this curved information
            NekInt64        entityid;  /// The entity id corresponding to the global edge/curve
            NekInt64        npoints;   /// The number of points in this curved entity. 
            NekInt64        ptid;      /// the id of point data map (currently always 0 since we are using just one set). 
            NekInt64        ptoffset;  /// point offset of data entry for this curve

            // An int (instead of a PointsType) defining the point
            // type from a PointsKey enum list. Since that we are
            // using a memory aligned structure which is suitable for
            // 32 and 64 bit machines.
            NekInt64        ptype; 
        };

        struct MeshCurvedPts
        {
            NekInt64 id;     /// id of this Point set

            /// Mapping to access the pts value. Given a 'ptoffset'
            /// value the npoints subsquent values provide the
            /// indexing on how to obtain the MeshVertex structure
            /// definiting the actually x,y,z values of each point in
            /// the curved entity. i.e. a list of edge values are found from
            //// pts[index[ptoffset +i] ]  0 <= i < npoints;
            std::vector<NekInt64>   index;  /// mapping to access pts value. 
            
            /// A list of MeshVertex entities containing the x,y,z
            /// values of unique points used in the curved entitites. 
            std::vector<MeshVertex> pts;
        };

        struct MeshEntity
        {
            int id;
            int origId;
            std::vector<unsigned int> list;
            bool ghost;
        };
    }
}
#endif

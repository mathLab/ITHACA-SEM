////////////////////////////////////////////////////////////////////////////////
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
namespace Nektar
{
    namespace LibUtilities
    {
        
        struct MeshEntity
        {
            int id;
            char type;
            std::vector<unsigned int> list;
        };
        
        struct MeshVertex
        {
            int id;
            NekDouble x;
            NekDouble y;
            NekDouble z;
        };
        
        struct MeshEdge
        {
            int id;
            int v0;
            int v1;
        };
        
        struct MeshFace
        {
            int id;
            std::vector<int> edgeList;
        };

        struct MeshElement
        {
            int id;
            char type;
            std::vector<int> list;
        };
        
        struct MeshCurved
        {
            int id;
            std::string entitytype;
            int entityid;
            std::string type;
            int npoints;
            std::string data;
        };
        
        struct MeshComposite
        {
            int id;
            char type;
            std::vector<int> list;
        };
    }
}
#endif

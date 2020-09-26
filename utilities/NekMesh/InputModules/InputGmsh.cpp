////////////////////////////////////////////////////////////////////////////////
//
//  File: InputGmsh.cpp
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
//  Description: GMSH converter.
//
////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <iostream>
using namespace std;

#include <tuple>

#include <NekMeshUtils/MeshElements/Element.h>
#include <NekMeshUtils/MeshElements/Point.h>
#include <NekMeshUtils/MeshElements/Line.h>
#include <NekMeshUtils/MeshElements/Triangle.h>
#include <NekMeshUtils/MeshElements/Quadrilateral.h>
#include <NekMeshUtils/MeshElements/Tetrahedron.h>
#include <NekMeshUtils/MeshElements/Pyramid.h>
#include <NekMeshUtils/MeshElements/Prism.h>
#include <NekMeshUtils/MeshElements/Hexahedron.h>

#include "InputGmsh.h"

using namespace Nektar::NekMeshUtils;

namespace Nektar
{
namespace Utilities
{

ModuleKey InputGmsh::className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eInputModule, "msh"), InputGmsh::create, "Reads Gmsh msh file.");

std::map<unsigned int, ElmtConfig> InputGmsh::elmMap = InputGmsh::GenElmMap();

/**
 * @brief Reorder Gmsh nodes so that they appear in a tensor product format
 * suitable for the interior of a Nektar++ quadrilateral.
 *
 * For an example, consider a second order quadrilateral. This routine will
 * produce a permutation map that applies the following permutation:
 *
 * ```
 *   3---6---2         6---7---8
 *   |       |         |       |
 *   7   8   5   --->  3   4   5
 *   |       |         |       |
 *   0---4---1         0---1---2
 *
 *     Gmsh          tensor-product
 * ```
 *
 * We assume that Gmsh uses a recursive ordering system, so that interior nodes
 * are reordered by calling this function recursively.
 *
 * @param nodes   The integer IDs of the nodes to be reordered, in Gmsh format.
 * @param n       The number of nodes in one coordinate direction. If this is
 *                zero we assume no reordering needs to be done and return the
 *                identity permutation.
 *
 * @return The nodes vector in tensor-product ordering.
 */
std::vector<int> quadTensorNodeOrdering(const std::vector<int> &nodes, int n)
{
    std::vector<int> nodeList;

    // Triangle
    nodeList.resize(nodes.size());

    // Vertices and edges
    nodeList[0] = nodes[0];
    if (n > 1)
    {
        nodeList[n - 1]     = nodes[1];
        nodeList[n * n - 1] = nodes[2];
        nodeList[n * (n - 1)] = nodes[3];
    }
    for (int i = 1; i < n - 1; ++i)
    {
        nodeList[i] = nodes[4 + i - 1];
    }
    for (int i = 1; i < n - 1; ++i)
    {
        nodeList[n * n - 1 - i] = nodes[4 + 2 * (n - 2) + i - 1];
    }

    // Interior (recursion)
    if (n > 2)
    {
        // Reorder interior nodes
        std::vector<int> interior((n - 2) * (n - 2));
        std::copy(
            nodes.begin() + 4 + 4 * (n - 2), nodes.end(), interior.begin());
        interior = quadTensorNodeOrdering(interior, n - 2);

        // Copy into full node list
        for (int j = 1; j < n - 1; ++j)
        {
            nodeList[j * n] = nodes[4 + 3 * (n - 2) + n - 2 - j];
            for (int i = 1; i < n - 1; ++i)
            {
                nodeList[j * n + i] = interior[(j - 1) * (n - 2) + (i - 1)];
            }
            nodeList[(j + 1) * n - 1] = nodes[4 + (n - 2) + j - 1];
        }
    }

    return nodeList;
}

/**
 * @brief Reorder Gmsh nodes so that they appear in a "tensor product" format
 * suitable for the interior of a Nektar++ triangle.
 *
 * For an example, consider a third-order triangle. This routine will produce a
 * permutation map that applies the following permutation:
 *
 * ```
 *   2                 9
 *   | \               | \
 *   7   6             7   7
 *   |     \           |     \
 *   8   9   5         4   5   6
 *   |         \       |         \
 *   0---3---4---1     0---1---2---3
 *
 *       Gmsh          tensor-product
 * ```
 *
 * We assume that Gmsh uses a recursive ordering system, so that interior nodes
 * are reordered by calling this function recursively.
 *
 * @param nodes   The integer IDs of the nodes to be reordered, in Gmsh format.
 * @param n       The number of nodes in one coordinate direction. If this is
 *                zero we assume no reordering needs to be done and return the
 *                identity permutation.
 *
 * @return The nodes vector in tensor-product ordering.
 */
std::vector<int> triTensorNodeOrdering(const std::vector<int> &nodes, int n)
{
    std::vector<int> nodeList;
    int cnt2;

    nodeList.resize(nodes.size());

    // Vertices
    nodeList[0] = nodes[0];
    if (n > 1)
    {
        nodeList[n - 1] = nodes[1];
        nodeList[n * (n + 1) / 2 - 1] = nodes[2];
    }

    // Edges
    int cnt = n;
    for (int i = 1; i < n - 1; ++i)
    {
        nodeList[i]               = nodes[3 + i - 1];
        nodeList[cnt]             = nodes[3 + 3 * (n - 2) - i];
        nodeList[cnt + n - i - 1] = nodes[3 + (n - 2) + i - 1];
        cnt += n - i;
    }

    // Interior (recursion)
    if (n > 3)
    {
        // Reorder interior nodes
        std::vector<int> interior((n - 3) * (n - 2) / 2);
        std::copy(
            nodes.begin() + 3 + 3 * (n - 2), nodes.end(), interior.begin());
        interior = triTensorNodeOrdering(interior, n - 3);

        // Copy into full node list
        cnt  = n;
        cnt2 = 0;
        for (int j = 1; j < n - 2; ++j)
        {
            for (int i = 0; i < n - j - 2; ++i)
            {
                nodeList[cnt + i + 1] = interior[cnt2 + i];
            }
            cnt += n - j;
            cnt2 += n - 2 - j;
        }
    }

    return nodeList;
}

typedef std::tuple<int, int, int> Mode;
struct cmpop
{
    bool operator()(Mode const &a, Mode const &b) const
    {
        if (std::get<0>(a) < std::get<0>(b))
        {
            return true;
        }
        if (std::get<0>(a) > std::get<0>(b))
        {
            return false;
        }
        if (std::get<1>(a) < std::get<1>(b))
        {
            return true;
        }
        if (std::get<1>(a) > std::get<1>(b))
        {
            return false;
        }
        if (std::get<2>(a) < std::get<2>(b))
        {
            return true;
        }

        return false;
    }
};

/**
 * @brief Reorder Gmsh nodes so that they appear in a "tensor product" format
 * suitable for the interior of a Nektar++ tetrahedron.
 *
 * For an example, consider a second order tetrahedron. This routine will
 * produce a permutation map that applies the following permutation:
 *
 * ```
 *               2                           5
 *             ,/|`\                       ,/|`\
 *           ,/  |  `\                   ,/  |  `\
 *         ,6    '.   `5               ,3    '.   `4
 *       ,/       8     `\           ,/       8     `\
 *     ,/         |       `\       ,/         |       `\
 *    0--------4--'.--------1     0--------1--'.--------2
 *     `\.         |      ,/       `\.         |      ,/
 *        `\.      |    ,9            `\.      |    ,7
 *           `7.   '. ,/                 `6.   '. ,/
 *              `\. |/                      `\. |/
 *                 `3                          `9
 *
 *             Gmsh                     tensor-product
 * ```
 *
 * We assume that Gmsh uses a recursive ordering system, so that interior nodes
 * are reordered by calling this function recursively.
 *
 * @param nodes   The integer IDs of the nodes to be reordered, in Gmsh format.
 * @param n       The number of nodes in one coordinate direction. If this is
 *                zero we assume no reordering needs to be done and return the
 *                identity permutation.
 *
 * @return The nodes vector in tensor-product ordering.
 */
std::vector<int> tetTensorNodeOrdering(const std::vector<int> &nodes, int n)
{
    std::vector<int> nodeList;
    int nTri = n*(n+1)/2;
    int nTet = n*(n+1)*(n+2)/6;

    nodeList.resize(nodes.size());
    nodeList[0] = nodes[0];

    if (n == 1)
    {
        return nodeList;
    }

    // Vertices
    nodeList[n - 1]    = nodes[1];
    nodeList[nTri - 1] = nodes[2];
    nodeList[nTet - 1] = nodes[3];

    if (n == 2)
    {
        return nodeList;
    }

    // Set up a map that takes (a,b,c) -> m to help us figure out where things
    // are inside the tetrahedron.
    std::map<Mode, int, cmpop> tmp;

    for (int k = 0, cnt = 0; k < n; ++k)
    {
        for (int j = 0; j < n - k; ++j)
        {
            for (int i = 0; i < n - k - j; ++i)
            {
                tmp[Mode(i,j,k)] = cnt++;
            }
        }
    }

    // Edges
    for (int i = 1; i < n-1; ++i)
    {
        int eI = i-1;
        nodeList[tmp[Mode(i,0,0)]]     = nodes[4 + eI];
        nodeList[tmp[Mode(n-1-i,i,0)]] = nodes[4 + (n-2) + eI];
        nodeList[tmp[Mode(0,n-1-i,0)]] = nodes[4 + 2*(n-2) + eI];
        nodeList[tmp[Mode(0,0,n-1-i)]] = nodes[4 + 3*(n-2) + eI];
        nodeList[tmp[Mode(0,i,n-1-i)]] = nodes[4 + 4*(n-2) + eI];
        nodeList[tmp[Mode(i,0,n-1-i)]] = nodes[4 + 5*(n-2) + eI];
    }

    if (n == 3)
    {
        return nodeList;
    }

    // For faces, we use the triTensorNodeOrdering routine to make our lives
    // slightly easier.
    int nFacePts = (n-3)*(n-2)/2;

    // Grab face points and reorder into a tensor-product type format
    vector<vector<int> > tmpNodes(4);
    int offset = 4 + 6*(n-2);

    for (int i = 0; i < 4; ++i)
    {
        tmpNodes[i].resize(nFacePts);
        for (int j = 0; j < nFacePts; ++j)
        {
            tmpNodes[i][j] = nodes[offset++];
        }
        tmpNodes[i] = triTensorNodeOrdering(tmpNodes[i], n-3);
    }

    if (n > 4)
    {
        // Now align faces
        vector<int> triVertId(3), toAlign(3);
        triVertId[0] = 0;
        triVertId[1] = 1;
        triVertId[2] = 2;

        // Faces 0,2: triangle verts {0,2,1} --> {0,1,2}
        HOTriangle<int> hoTri(triVertId, tmpNodes[0]);
        toAlign[0] = 0;
        toAlign[1] = 2;
        toAlign[2] = 1;

        hoTri.Align(toAlign);
        tmpNodes[0] = hoTri.surfVerts;

        hoTri.surfVerts = tmpNodes[2];
        hoTri.Align(toAlign);
        tmpNodes[2] = hoTri.surfVerts;

        // Face 3: triangle verts {1,2,0} --> {0,1,2}
        toAlign[0] = 1;
        toAlign[1] = 2;
        toAlign[2] = 0;

        hoTri.surfVerts = tmpNodes[3];
        hoTri.Align(toAlign);
        tmpNodes[3] = hoTri.surfVerts;
    }

    // Now apply faces. Note that faces 3 and 2 are swapped between Gmsh and
    // Nektar++ order.
    for (int j = 1, cnt = 0; j < n-2; ++j)
    {
        for (int i = 1; i < n-j-1; ++i, ++cnt)
        {
            nodeList[tmp[Mode(i,j,0)]]       = tmpNodes[0][cnt];
            nodeList[tmp[Mode(i,0,j)]]       = tmpNodes[1][cnt];
            nodeList[tmp[Mode(n-1-i-j,i,j)]] = tmpNodes[3][cnt];
            nodeList[tmp[Mode(0,i,j)]]       = tmpNodes[2][cnt];
        }
    }

    if (n == 4)
    {
        return nodeList;
    }

    // Finally, recurse on interior volume
    vector<int> intNodes, tmpInt;
    for (int i = offset; i < nTet; ++i)
    {
        intNodes.push_back(nodes[i]);
    }
    tmpInt = tetTensorNodeOrdering(intNodes, n-4);

    for (int k = 1, cnt = 0; k < n - 2; ++k)
    {
        for (int j = 1; j < n - k - 1; ++j)
        {
            for (int i = 1; i < n - k - j - 1; ++i)
            {
                nodeList[tmp[Mode(i,j,k)]] = tmpInt[cnt++];
            }
        }
    }

    return nodeList;
}

/**
 * @brief Reorder Gmsh nodes so that they appear in a "tensor product" format
 * suitable for the interior of a Nektar++ prism. This routine is specifically
 * designed for *interior* Gmsh nodes only.
 *
 * Prisms are a bit of a special case, in that interior nodes are heterogeneous
 * in the number of points they use. As an example, for a second-order interior
 * of a fourth-order prism, this routine calculates the mapping taking us
 *
 * ```
 *               ,6                        ,8
 *            .-' | \                   .-' | \
 *         ,-3    |   \              ,-5    |   \
 *      ,-'  | \  |     \         ,-'  | \  |     \
 *     2     |   \7-------8      2     |   \6-------7
 *     | \   | ,-' \   ,-'       | \   | ,-' \   ,-'
 *     |   \,5------,0'          |   \,3------,4'
 *     |,-'  \   ,-'             |,-'  \   ,-'
 *     1-------4'                0-------1'
 *
 *             Gmsh                 tensor-product
 * ```
 *
 * @todo Make this routine work for higher orders. The Gmsh ordering seems
 *       a little different than other elements, therefore the functionality is
 *       (for now) hard-coded to orders less than 5.
 *
 * @param nodes   The integer IDs of the nodes to be reordered, in Gmsh format.
 * @param n       The number of nodes in one coordinate direction. If this is
 *                zero we assume no reordering needs to be done and return the
 *                identity permutation.
 *
 * @return The nodes vector in tensor-product ordering.
 */
std::vector<int> prismTensorNodeOrdering(const std::vector<int> &nodes, int n)
{
    std::vector<int> nodeList;

    if (n == 0)
    {
        return nodeList;
    }

    nodeList.resize(nodes.size());

    if (n == 2)
    {
        nodeList[0] = nodes[1];
        nodeList[1] = nodes[0];
        return nodeList;
    }

    // For some reason, this ordering is different. Whereas Gmsh usually orders
    // vertices first, followed by edges, this ordering goes VVE-VVE-VVE for the
    // three edges that contain edge-interior information, but also vertex
    // orientation is not the same as the original Gmsh prism. The below is done
    // by looking at the ordering by hand and is hence why we don't support
    // order > 4 right now.
    if (n == 3)
    {
        nodeList[0] = nodes[1];
        nodeList[1] = nodes[4];
        nodeList[2] = nodes[2];
        nodeList[3] = nodes[5];
        nodeList[4] = nodes[0];
        nodeList[5] = nodes[3];
        nodeList[6] = nodes[7];
        nodeList[7] = nodes[8];
        nodeList[8] = nodes[6];
    }

    ASSERTL0(n < 4, "Prism Gmsh input and output is incomplete for orders "
                    "larger than 4");

    return nodeList;
}

/**
 * @brief Reorder Gmsh nodes so that they appear in a "tensor product" format
 * suitable for the interior of a Nektar++ hexahedron.
 *
 * For an example, consider a second order hexahedron. This routine will produce
 * a permutation map that applies the following permutation:
 *
 * ```
 *    3----13----2         6----7-----8
 *    |\         |\        |\         |\
 *    |15    24  | 14      |15    16  | 18
 *    9  \ 20    11 \      3  \  4    5  \
 *    |   7----19+---6     |   24---25+---26
 *    |22 |  26  | 23|     |12 |  13  | 14|
 *    0---+-8----1   |     0---+-1----2   |
 *     \ 17    25 \  18     \ 21    22 \  23
 *     10 |  21    12|       9 |  10    11|
 *       \|         \|        \|         \|
 *        4----16----5         18---19----20
 *
 *          Gmsh            tensor-product
 * ```
 *
 * We assume that Gmsh uses a recursive ordering system, so that interior nodes
 * are reordered by calling this function recursively.
 *
 * @param nodes   The integer IDs of the nodes to be reordered, in Gmsh format.
 * @param n       The number of nodes in one coordinate direction. If this is
 *                zero we assume no reordering needs to be done and return the
 *                identity permutation.
 *
 * @return The nodes vector in tensor-product ordering.
 */
std::vector<int> hexTensorNodeOrdering(const std::vector<int> &nodes, int n)
{
    int i, j, k;
    std::vector<int> nodeList;

    nodeList.resize(nodes.size());
    nodeList[0] = nodes[0];

    if (n == 1)
    {
        return nodeList;
    }

    // Vertices: same order as Nektar++
    nodeList[n - 1]               = nodes[1];
    nodeList[n*n -1]              = nodes[2];
    nodeList[n*(n-1)]             = nodes[3];
    nodeList[n*n*(n-1)]           = nodes[4];
    nodeList[n - 1 + n*n*(n-1)]   = nodes[5];
    nodeList[n*n -1 + n*n*(n-1)]  = nodes[6];
    nodeList[n*(n-1) + n*n*(n-1)] = nodes[7];

    if (n == 2)
    {
        return nodeList;
    }

    int hexEdges[12][2] = {
        { 0, 1 }, { n-1, n }, { n*n-1, -1 }, { n*(n-1), -n },
        { 0, n*n }, { n-1, n*n }, { n*n - 1, n*n }, { n*(n-1), n*n },
        { n*n*(n-1), 1 }, { n*n*(n-1) + n-1, n }, { n*n*n-1, -1 },
        { n*n*(n-1) + n*(n-1), -n }
    };
    int hexFaces[6][3] = {
        { 0, 1, n }, { 0, 1, n*n }, { n-1, n, n*n },
        { n*(n-1), 1, n*n }, { 0, n, n*n }, { n*n*(n-1), 1, n }
    };
    int gmshToNekEdge[12] = {0, -3, 4, 1, 5, 2, 6, 7, 8, -11, 9, 10};

    // Edges
    int offset = 8;
    for (int i = 0; i < 12; ++i)
    {
        int e = abs(gmshToNekEdge[i]);

        if (gmshToNekEdge[i] >= 0)
        {
            for (int j = 1; j < n-1; ++j)
            {
                nodeList[hexEdges[e][0] + j*hexEdges[e][1]] = nodes[offset++];
            }
        }
        else
        {
            for (int j = 1; j < n-1; ++j)
            {
                nodeList[hexEdges[e][0] + (n-j-1)*hexEdges[e][1]] = nodes[offset++];
            }
        }
    }

    // Faces
    int gmsh2NekFace[6] = {0, 1, 4, 2, 3, 5};

    // Map which defines orientation between Gmsh and Nektar++ faces.
    StdRegions::Orientation faceOrient[6] = {
        StdRegions::eDir1FwdDir2_Dir2FwdDir1,
        StdRegions::eDir1FwdDir1_Dir2FwdDir2,
        StdRegions::eDir1FwdDir2_Dir2FwdDir1,
        StdRegions::eDir1FwdDir1_Dir2FwdDir2,
        StdRegions::eDir1BwdDir1_Dir2FwdDir2,
        StdRegions::eDir1FwdDir1_Dir2FwdDir2};

    for (i = 0; i < 6; ++i)
    {
        int n2 = (n-2)*(n-2);
        int face = gmsh2NekFace[i];
        offset   = 8 + 12 * (n-2) + i * n2;

        // Create a list of interior face nodes for this face only.
        vector<int> faceNodes(n2);
        for (j = 0; j < n2; ++j)
        {
            faceNodes[j] = nodes[offset + j];
        }

        // Now get the reordering of this face, which puts Gmsh
        // recursive ordering into Nektar++ row-by-row order.
        faceNodes = quadTensorNodeOrdering(faceNodes, n-2);
        vector<int> tmp(n2);

        // Finally reorient the face according to the geometry
        // differences.
        if (faceOrient[i] == StdRegions::eDir1FwdDir1_Dir2FwdDir2)
        {
            // Orientation is the same, just copy.
            tmp = faceNodes;
        }
        else if (faceOrient[i] == StdRegions::eDir1FwdDir2_Dir2FwdDir1)
        {
            // Tranposed faces
            for (j = 0; j < n-2; ++j)
            {
                for (k = 0; k < n-2; ++k)
                {
                    tmp[j * (n-2) + k] = faceNodes[k * (n-2) + j];
                }
            }
        }
        else if (faceOrient[i] == StdRegions::eDir1BwdDir1_Dir2FwdDir2)
        {
            for (j = 0; j < n-2; ++j)
            {
                for (k = 0; k < n-2; ++k)
                {
                    tmp[j * (n-2) + k] = faceNodes[j * (n-2) + (n - k - 3)];
                }
            }
        }

        // Now put this into the right place in the output array
        for (k = 1; k < n-1; ++k)
        {
            for (j = 1; j < n-1; ++j)
            {
                nodeList[hexFaces[face][0] + j*hexFaces[face][1] + k*hexFaces[face][2]]
                    = faceNodes[(k-1)*(n-2) + j-1];
            }
        }
    }

    // Finally, recurse on interior volume
    vector<int> intNodes, tmpInt;
    for (int i = 8 + 12 * (n-2) + 6 * (n-2) * (n-2); i < n*n*n; ++i)
    {
        intNodes.push_back(nodes[i]);
    }

    if (intNodes.size())
    {
        tmpInt = hexTensorNodeOrdering(intNodes, n-2);
        for (int k = 1, cnt = 0; k < n - 1; ++k)
        {
            for (int j = 1; j < n - 1; ++j)
            {
                for (int i = 1; i < n - 1; ++i)
                {
                    nodeList[i + j * n + k * n * n] = tmpInt[cnt++];
                }
            }
        }
    }

    return nodeList;
}


/**
 * @brief Set up InputGmsh object.
 *
 */
InputGmsh::InputGmsh(MeshSharedPtr m)
    : InputModule(m), m_version(0.0), m_prevId(-1), m_maxTagId(-1)
{
}

InputGmsh::~InputGmsh()
{
}

/**
 * @brief Representation of Gmsh entity so that we can extract physical tag IDs.
 */
struct GmshEntity
{
    double minX, minY, minZ;
    double maxX, maxY, maxZ;
    std::vector<int> physicalTags;
};

/**
 * Gmsh file contains a list of nodes and their coordinates, along with a list
 * of elements and those nodes which define them. We read in and store the list
 * of nodes in #m_node and store the list of elements in #m_element. Each new
 * element is supplied with a list of entries from m_node which defines the
 * element. Finally some mesh statistics are printed.
 *
 * @param   pFilename           Filename of Gmsh file to read.
 */
void InputGmsh::Process()
{
    // Open the file stream.
    OpenStream();

    m_mesh->m_expDim   = 0;
    m_mesh->m_spaceDim = 0;
    string line;
    int fileType      = 0;
    int nVBlocks      = 0;
    int nVertices     = 0;
    int nEBlocks      = 0;
    int nElements     = 0;
    int tag           = 0;
    int elm_type      = 0;
    int tmp           = 0;

    if (m_mesh->m_verbose)
    {
        cout << "InputGmsh: Start reading file..." << endl;
    }

    std::vector<std::map<int, GmshEntity>> entityMap(4);

    while (!m_mshFile.eof())
    {
        getline(m_mshFile, line);
        stringstream s(line);
        string word;
        s >> word;

        // Process file format.
        if (word == "$MeshFormat")
        {
            getline(m_mshFile, line);
            stringstream s(line);
            s >> m_version;
            s >> fileType;
            ASSERTL0(fileType == 0, "Cannot read binary Gmsh files.")
            ASSERTL0(m_version <= 4.1, ".msh file format versions greater than 4.1 are not currently supported.")
        }
        // Process entities (v4+)
        else if (word == "$Entities")
        {
            size_t nEntity[4];
            getline(m_mshFile, line);
            stringstream s(line);

            s >> nEntity[0] >> nEntity[1] >> nEntity[2] >> nEntity[3];

            for (int i = 0; i < 4; ++i)
            {
                for (int j = 0; j < nEntity[i]; ++j)
                {
                    getline(m_mshFile, line);
                    stringstream si(line);
                    GmshEntity ent;

                    int entityId;
                    si >> entityId >> ent.minX >> ent.minY >> ent.minZ;

                    if (i > 0)
                    {
                        si >> ent.maxX >> ent.maxY >> ent.maxZ;
                    }

                    int nPhysTags, tmp;
                    si >> nPhysTags;

                    for (int k = 0; k < nPhysTags; ++k)
                    {
                        si >> tmp;
                        ent.physicalTags.push_back(tmp);
                    }

                    entityMap[i][entityId] = ent;
                }
            }
        }
        // Process nodes.
        else if (word == "$Nodes")
        {
            getline(m_mshFile, line);
            stringstream s(line);

            if (m_version >= 4.0)
            {
                s >> nVBlocks;

                for (int i = 0; i < nVBlocks; ++i)
                {
                    getline(m_mshFile, line);
                    stringstream si(line);
                    si >> tmp >> tmp >> tmp >> nVertices;

                    if (m_version == 4.0)
                    {
                        for (int i = 0; i < nVertices; ++i)
                        {
                            ReadNextNode();
                        }
                    }
                    else if (m_version == 4.1)
                    {
                        ReadNextNodeBlock(nVertices);
                    }
                }
            }
            else
            {
                s >> nVertices;

                for (int i = 0; i < nVertices; ++i)
                {
                    ReadNextNode();
                }
            }
        }
        // Process elements
        else if (word == "$Elements")
        {
            getline(m_mshFile, line);
            stringstream s(line);

            if (m_version >= 4.0)
            {
                s >> nEBlocks;

                for (int i = 0; i < nEBlocks; ++i)
                {
                    getline(m_mshFile, line);
                    stringstream si(line);

                    int tagDim;
                    if (m_version == 4.0)
                    {
                        si >> tag >> tagDim >> elm_type >> nElements;
                    }
                    else
                    {
                        si >> tagDim >> tag >> elm_type >> nElements;
                    }
                    // Query tag in map & don't bother constructing non-physical
                    // surfaces.
                    std::vector<int> physIds = entityMap[tagDim][tag].physicalTags;

                    if (physIds.size() == 0)
                    {
                        for (int j = 0; j < nElements; ++j)
                        {
                            getline(m_mshFile, line);
                        }
                    }
                    else
                    {
                        tag = physIds[0];
                        for (int j = 0; j < nElements; ++j)
                        {
                            ReadNextElement(tag, elm_type);
                        }
                    }

                }
            }
            else
            {
                s >> nElements;
                for (int i = 0; i < nElements; ++i)
                {
                    ReadNextElement();
                }
            }
        }
    }
    m_mshFile.reset();

    // Go through element and remap tags if necessary.
    map<int, map<LibUtilities::ShapeType, int> > compMap;

    for (int i = 0; i < m_mesh->m_element[m_mesh->m_expDim].size(); ++i)
    {
        ElementSharedPtr el          = m_mesh->m_element[m_mesh->m_expDim][i];
        LibUtilities::ShapeType type = el->GetConf().m_e;

        vector<int> tags = el->GetTagList();
        int tag          = tags[0];

        auto cIt = compMap.find(tag);

        if (cIt == compMap.end())
        {
            compMap[tag][type] = tag;
            continue;
        }

        // Reset tag for this element.
        auto sIt = cIt->second.find(type);
        if (sIt == cIt->second.end())
        {
            m_maxTagId++;
            cIt->second[type] = m_maxTagId;
            tags[0] = m_maxTagId;
            el->SetTagList(tags);
        }
        else if (sIt->second != tag)
        {
            tags[0] = sIt->second;
            el->SetTagList(tags);
        }
    }

    bool printInfo = false;
    for (auto &cIt : compMap)
    {
        if (cIt.second.size() > 1)
        {
            printInfo = true;
            break;
        }
    }

    if (printInfo)
    {
        cout << "Multiple elements in composite detected; remapped:" << endl;
        for (auto &cIt : compMap)
        {
            if (cIt.second.size() > 1)
            {
                auto sIt = cIt.second.begin();
                cout << "- Tag " << cIt.first << " => " << sIt->second << " ("
                     << LibUtilities::ShapeTypeMap[sIt->first] << ")";
                sIt++;

                for (; sIt != cIt.second.end(); ++sIt)
                {
                    cout << ", " << sIt->second << " ("
                         << LibUtilities::ShapeTypeMap[sIt->first] << ")";
                }

                cout << endl;
            }
        }
    }

    // Process rest of mesh.
    ProcessVertices();
    ProcessEdges();
    ProcessFaces();
    ProcessElements();
    ProcessComposites();
}

/**
 * Read in next node block for v4 format
 */
void InputGmsh::ReadNextNodeBlock(int nVertices)
{
    string line;
    double x = 0, y = 0, z = 0;
    vector<int> id(nVertices);

    for (int i = 0; i < nVertices; ++i)
    {
        getline(m_mshFile, line);
        stringstream st(line);
        st >> id[i];
    }

    for (int i = 0; i < nVertices; ++i)
    {
        getline(m_mshFile, line);
        stringstream st(line);
        st >> x >> y >> z;

        SaveNode(id[i], x, y, z);
    }
}

/**
 * Read in next node
 */
void InputGmsh::ReadNextNode()
{
    string line;
    getline(m_mshFile, line);
    stringstream st(line);
    double x = 0, y = 0, z = 0;
    int id = 0;
    st >> id >> x >> y >> z;

    SaveNode(id, x, y, z);
}

/**
 * Save node into mesh
 */
void InputGmsh::SaveNode(int id, NekDouble x, NekDouble y, NekDouble z)
{
    if ((x * x) > 0.000001 && m_mesh->m_spaceDim < 1)
    {
        m_mesh->m_spaceDim = 1;
    }
    if ((y * y) > 0.000001 && m_mesh->m_spaceDim < 2)
    {
        m_mesh->m_spaceDim = 2;
    }
    if ((z * z) > 0.000001 && m_mesh->m_spaceDim < 3)
    {
        m_mesh->m_spaceDim = 3;
    }

    id -= 1; // counter starts at 0

    if (!m_idMap.size())
    {
        if (id - m_prevId == 1)
        {
            m_prevId = id;
        }
        else
        {
            // Build m_idMap so far
            for (int i = 0; i < m_mesh->m_node.size(); ++i)
            {
                m_idMap[i] = i;
            }

            m_idMap[id] = m_mesh->m_node.size();
        }
    }
    else
    {
        m_idMap[id] = m_mesh->m_node.size();
    }

    m_mesh->m_node.push_back(std::shared_ptr<Node>(new Node(id, x, y, z)));
}

/**
 * Read in next element
 */
void InputGmsh::ReadNextElement(int tag, int elm_type)
{
    string line;
    getline(m_mshFile, line);
    stringstream st(line);
    int id = 0, num_tag = 0, num_nodes = 0;

    st >> id;
    id -= 1; // counter starts at 0

    if (m_version < 4.0)
    {
        st >> elm_type >> num_tag;
    }

    auto it = elmMap.find(elm_type);
    if (it == elmMap.end())
    {
        cerr << "Error: element type " << elm_type << " not supported" << endl;
        abort();
    }

    // Read element tags (version 2 only)
    vector<int> tags;
    for (int j = 0; j < num_tag; ++j)
    {
        int tag = 0;
        st >> tag;
        tags.push_back(tag);
    }
    if (m_version >= 4.0)
    {
        tags.push_back(tag);
    }
    tags.resize(1);

    m_maxTagId = max(m_maxTagId, tags[0]);

    // Read element node list
    vector<NodeSharedPtr> nodeList;
    num_nodes = GetNnodes(elm_type);
    for (int k = 0; k < num_nodes; ++k)
    {
        int node = 0;
        st >> node;
        node -= 1; // counter starts at 0
        if (!m_idMap.size())
        {
            nodeList.push_back(m_mesh->m_node[node]);
        }
        else
        {
            nodeList.push_back(m_mesh->m_node[m_idMap[node]]);
        }
    }

    // Look up reordering.
    auto oIt = m_orderingMap.find(elm_type);

    // If it's not created, then create it.
    if (oIt == m_orderingMap.end())
    {
        oIt =
            m_orderingMap.insert(make_pair(elm_type, CreateReordering(elm_type)))
                .first;
    }

    // Apply reordering map where necessary.
    if (oIt->second.size() > 0)
    {
        vector<int> &mapping      = oIt->second;
        vector<NodeSharedPtr> tmp = nodeList;
        for (int i = 0; i < mapping.size(); ++i)
        {
            nodeList[i] = tmp[mapping[i]];
        }
    }

    // Create element
    ElementSharedPtr E = GetElementFactory().CreateInstance(
        it->second.m_e, it->second, nodeList, tags);

    // Determine mesh expansion dimension
    if (E->GetDim() > m_mesh->m_expDim)
    {
        m_mesh->m_expDim = E->GetDim();
    }
    m_mesh->m_element[E->GetDim()].push_back(E);
}

/**
 * For a given msh ID, return the corresponding number of nodes.
 */
int InputGmsh::GetNnodes(unsigned int InputGmshEntity)
{
    int nNodes;

    auto it = elmMap.find(InputGmshEntity);

    if (it == elmMap.end())
    {
        cerr << "Unknown element type " << InputGmshEntity << endl;
        abort();
    }

    switch (it->second.m_e)
    {
        case LibUtilities::ePoint:
            nNodes = Point::GetNumNodes(it->second);
            break;
        case LibUtilities::eSegment:
            nNodes = Line::GetNumNodes(it->second);
            break;
        case LibUtilities::eTriangle:
            nNodes = Triangle::GetNumNodes(it->second);
            break;
        case LibUtilities::eQuadrilateral:
            nNodes = Quadrilateral::GetNumNodes(it->second);
            break;
        case LibUtilities::eTetrahedron:
            nNodes = Tetrahedron::GetNumNodes(it->second);
            it->second.m_faceCurveType = LibUtilities::eNodalTriEvenlySpaced;
            break;
        case LibUtilities::ePyramid:
            nNodes = Pyramid::GetNumNodes(it->second);
            break;
        case LibUtilities::ePrism:
            nNodes = Prism::GetNumNodes(it->second);
            break;
        case LibUtilities::eHexahedron:
            nNodes = Hexahedron::GetNumNodes(it->second);
            break;
        default:
            cerr << "Unknown element type!" << endl;
            abort();
            break;
    }

    return nNodes;
}

/**
 * @brief Create a reordering map for a given element.
 *
 * Since Gmsh and Nektar++ have different vertex, edge and face
 * orientations, we need to reorder the nodes in a Gmsh MSH file so that
 * they work with the Nektar++ orderings, since this is what is used in
 * the elements defined in the converter.
 */
vector<int> InputGmsh::CreateReordering(unsigned int InputGmshEntity)
{
    auto it = elmMap.find(InputGmshEntity);

    if (it == elmMap.end())
    {
        cerr << "Unknown element type " << InputGmshEntity << endl;
        abort();
    }

    // For specific elements, call the appropriate function to perform
    // the renumbering.
    switch (it->second.m_e)
    {
        case LibUtilities::eTriangle:
            return TriReordering(it->second);
            break;
        case LibUtilities::eQuadrilateral:
            return QuadReordering(it->second);
            break;
        case LibUtilities::eTetrahedron:
            return TetReordering(it->second);
            break;
        case LibUtilities::ePrism:
            return PrismReordering(it->second);
            break;
        case LibUtilities::eHexahedron:
            return HexReordering(it->second);
            break;
        case LibUtilities::eSegment:
            return LineReordering(it->second);
        default:
            break;
    }

    // Default: no reordering.
    vector<int> returnVal;
    return returnVal;
}

vector<int> InputGmsh::LineReordering(ElmtConfig conf)
{
    const int order = conf.m_order;

    vector<int> mapping(order+1);
    for(int i = 0; i < order+1; i++)
    {
        mapping[i] = i;
    }

    return mapping;
}

/**
 * @brief Create a reordering for triangles.
 */
vector<int> InputGmsh::TriReordering(ElmtConfig conf)
{
    const int order = conf.m_order;
    const int n     = order - 1;

    // Copy vertices.
    vector<int> mapping(3);
    for (int i = 0; i < 3; ++i)
    {
        mapping[i] = i;
    }

    if (order == 1)
    {
        return mapping;
    }

    // Curvilinear edges.
    mapping.resize(3 + 3 * n);

    for (int i = 3; i < 3 + 3 * n; ++i)
    {
        mapping[i] = i;
    }

    if (!conf.m_faceNodes)
    {
        return mapping;
    }

    // Interior nodes.
    vector<int> interior(n * (n - 1) / 2);
    for (int i = 0; i < interior.size(); ++i)
    {
        interior[i] = i + 3 + 3 * n;
    }

    if (interior.size() > 0)
    {
        interior = triTensorNodeOrdering(interior, n - 1);
    }

    mapping.insert(mapping.end(), interior.begin(), interior.end());
    return mapping;
}

/**
 * @brief Create a reordering for quadrilaterals.
 */
vector<int> InputGmsh::QuadReordering(ElmtConfig conf)
{
    const int order = conf.m_order;
    const int n     = order - 1;

    // Copy vertices.
    vector<int> mapping(4);
    for (int i = 0; i < 4; ++i)
    {
        mapping[i] = i;
    }

    if (order == 1)
    {
        return mapping;
    }

    // Curvilinear edges.
    mapping.resize(4 + 4 * n);

    for (int i = 4; i < 4 + 4 * n; ++i)
    {
        mapping[i] = i;
    }

    if (!conf.m_faceNodes)
    {
        return mapping;
    }

    // Interior nodes.
    vector<int> interior(n * n);
    for (int i = 0; i < interior.size(); ++i)
    {
        interior[i] = i + 4 + 4 * n;
    }

    if (interior.size() > 0)
    {
        interior = quadTensorNodeOrdering(interior, n);
    }
    mapping.insert(mapping.end(), interior.begin(), interior.end());
    return mapping;
}

/**
 * @brief Create a reordering for tetrahedra.
 */
vector<int> InputGmsh::TetReordering(ElmtConfig conf)
{
    const int order = conf.m_order;
    const int n     = order - 1;
    const int n2    = n * (n - 1) / 2;

    int i, j;
    vector<int> mapping(4);

    // Copy vertices.
    for (i = 0; i < 4; ++i)
    {
        mapping[i] = i;
    }

    if (order == 1)
    {
        return mapping;
    }

    // Curvilinear edges.
    mapping.resize(4 + 6 * n);

    // Curvilinear edges.
    static int gmshToNekEdge[6] = {0, 1, 2, 3, 5, 4};
    static int gmshToNekRev[6]  = {0, 0, 1, 1, 1, 1};

    // Reorder edges.
    int offset, cnt = 4;
    for (i = 0; i < 6; ++i)
    {
        offset = 4 + n * gmshToNekEdge[i];

        if (gmshToNekRev[i])
        {
            for (int j = 0; j < n; ++j)
            {
                mapping[offset + n - j - 1] = cnt++;
            }
        }
        else
        {
            for (int j = 0; j < n; ++j)
            {
                mapping[offset + j] = cnt++;
            }
        }
    }

    if (conf.m_faceNodes == false || n2 == 0)
    {
        return mapping;
    }

    // Curvilinear faces.
    mapping.resize(4 + 6 * n + 4 * n2);

    static int gmshToNekFace[4] = {0, 1, 3, 2};

    vector<int> triVertId(3);
    triVertId[0] = 0;
    triVertId[1] = 1;
    triVertId[2] = 2;

    // Loop over Gmsh faces
    for (i = 0; i < 4; ++i)
    {
        int face    = gmshToNekFace[i];
        int offset2 = 4 + 6 * n + i * n2;
        offset      = 4 + 6 * n + face * n2;

        // Create a list of interior face nodes for this face only.
        vector<int> faceNodes(n2);
        vector<int> toAlign(3);
        for (j = 0; j < n2; ++j)
        {
            faceNodes[j] = offset2 + j;
        }

        // Now get the reordering of this face, which puts Gmsh
        // recursive ordering into Nektar++ row-by-row order.
        vector<int> tmp = triTensorNodeOrdering(faceNodes, n - 1);
        HOTriangle<int> hoTri(triVertId, tmp);

        // Apply reorientation
        if (i == 0 || i == 2)
        {
            // Triangle verts {0,2,1} --> {0,1,2}
            toAlign[0] = 0;
            toAlign[1] = 2;
            toAlign[2] = 1;
            hoTri.Align(toAlign);
        }
        else if (i == 3)
        {
            // Triangle verts {1,2,0} --> {0,1,2}
            toAlign[0] = 1;
            toAlign[1] = 2;
            toAlign[2] = 0;
            hoTri.Align(toAlign);
        }

        // Fill in mapping.
        for (j = 0; j < n2; ++j)
        {
            mapping[offset + j] = hoTri.surfVerts[j];
        }
    }

    if (conf.m_volumeNodes == false)
    {
        return mapping;
    }

    const int nInt = (order - 3) * (order - 2) * (order - 1) / 6;
    if (nInt <= 0)
    {
        return mapping;
    }

    if (nInt == 1)
    {
        mapping.push_back(mapping.size());
        return mapping;
    }

    int ntot = (order+1)*(order+2)*(order+3)/6;
    vector<int> interior;

    for (int i = 4 + 6 * n + 4 * n2; i < ntot; ++i)
    {
        interior.push_back(i);
    }

    if (interior.size() > 0)
    {
        interior = tetTensorNodeOrdering(interior, order-3);
    }

    mapping.insert(mapping.end(), interior.begin(), interior.end());

    return mapping;
}

/**
 * @brief Create a reordering for prisms.
 *
 * Note that whilst Gmsh MSH files have the capability to support
 * high-order prisms, presently Gmsh does not seem to be capable of
 * generating higher than second-order prismatic meshes, so most of the
 * following is untested.
 */
vector<int> InputGmsh::PrismReordering(ElmtConfig conf)
{
    const int order = conf.m_order;
    const int n     = order - 1;

    int i;
    vector<int> mapping(6);

    // To get from Gmsh -> Nektar++ prism, coordinates axes are
    // different; need to mirror in the triangular faces, and then
    // reorder vertices to make ordering anticlockwise on base quad.
    static int nekToGmshVerts[6] = {3, 4, 1, 0, 5, 2};
    // Inverse (gmsh vert -> Nektar++ vert): 3->0, 4->1, 1->2, 0->3, 5->4, 2->5

    for (i = 0; i < 6; ++i)
    {
        mapping[i] = nekToGmshVerts[i];
    }

    if (order == 1)
    {
        return mapping;
    }

    // Curvilinear edges.
    mapping.resize(6 + 9 * n);

    // Nektar++ edges:
    //   0 =    1 =    2 =    3 =    4 =    5 =    6 =    7 =    8 =
    //   {0,1}, {1,2}, {3,2}, {0,3}, {0,4}, {1,4}, {2,5}, {3,5}, {4,5}
    // Gmsh edges (from Geo/MPrism.h of Gmsh source):
    //   {0,1}, {0,2}, {0,3}, {1,2}, {1,4}, {2,5}, {3,4}, {3,5}, {4,5}
    // Apply inverse of gmshToNekVerts map:
    //   {3,2}, {3,5}, {3,0}, {2,5}, {2,1}, {5,4}, {0,1}, {0,4}, {1,4}
    // Nektar++ mapping (negative indicates reverse orientation):
    //   = 2    = 7    = -3   = 6    = -1   = -8   = 0    = 4    = 5
    static int gmshToNekEdge[9] = {2, 7, 3, 6, 1, 8, 0, 4, 5};
    static int gmshToNekRev[9]  = {0, 0, 1, 0, 1, 1, 0, 0, 0};

    // Reorder edges.
    int offset, cnt = 6;
    for (i = 0; i < 9; ++i)
    {
        offset = 6 + n * gmshToNekEdge[i];

        if (gmshToNekRev[i])
        {
            for (int j = 0; j < n; ++j)
            {
                mapping[offset + n - j - 1] = cnt++;
            }
        }
        else
        {
            for (int j = 0; j < n; ++j)
            {
                mapping[offset + j] = cnt++;
            }
        }
    }

    if (conf.m_faceNodes == false)
    {
        return mapping;
    }

    int nTriInt  = n * (n - 1) / 2;
    int nQuadInt = n * n;

    // Nektar++ faces:
    //   0 = {0,1,2,3}, 1 = {0,1,4}, 2 = {1,2,5,4}, 3 = {3,2,5}, 4 = {0,3,5,4}
    // Gmsh faces (from Geo/MPrism.h of Gmsh source):
    //   {0,2,1}, {3,4,5}, {0,1,4,3}, {0,3,5,2}, {1,2,5,4}
    // Apply inverse of gmshToNekVerts map:
    //   {3,5,2}, {0,1,4}, {3,2,1,0}, {3,0,4,5}, {2,5,4,1}
    //   = 3      = 1      = 0        = 4        = 2
    // This gives gmsh -> Nektar++ faces:
    static int gmshToNekFace[5] = {3, 1, 0, 4, 2};

    // Face offsets
    vector<int> offsets(5), offsets2(5);
    offset = 6 + 9 * n;

    // Offsets in the gmsh order: face ordering is TTQQQ
    offsets[0] = offset;
    offsets[1] = offset + nTriInt;
    offsets[2] = offset + 2 * nTriInt;
    offsets[3] = offset + 2 * nTriInt + nQuadInt;
    offsets[4] = offset + 2 * nTriInt + 2 * nQuadInt;

    // Offsets in the Nektar++ order: face ordering is QTQTQ
    offsets2[0] = offset;
    offsets2[1] = offset + nQuadInt;
    offsets2[2] = offset + nQuadInt + nTriInt;
    offsets2[3] = offset + 2 * nQuadInt + nTriInt;
    offsets2[4] = offset + 2 * nQuadInt + 2 * nTriInt;

    mapping.resize(6 + 9 * n + 3 * nQuadInt + 2 * nTriInt);

    offset = 6 + 9 * n;

    vector<int> triVertId(3);
    triVertId[0] = 0;
    triVertId[1] = 1;
    triVertId[2] = 2;

    vector<int> quadVertId(4);
    quadVertId[0] = 0;
    quadVertId[1] = 1;
    quadVertId[2] = 2;
    quadVertId[3] = 3;

    for (i = 0; i < 5; ++i)
    {
        int face    = gmshToNekFace[i];
        int offset2 = offsets[i];
        offset      = offsets2[face];

        bool tri = i < 2;
        int nFacePts = tri ? nTriInt : nQuadInt;

        if (nFacePts == 0)
        {
            continue;
        }

        // Create a list of interior face nodes for this face only.
        vector<int> faceNodes(nFacePts);
        vector<int> toAlign(tri ? 3 : 4);
        for (int j = 0; j < nFacePts; ++j)
        {
            faceNodes[j] = offset2 + j;
        }

        if (tri)
        {
            // Now get the reordering of this face, which puts Gmsh
            // recursive ordering into Nektar++ row-by-row order.
            vector<int> tmp = triTensorNodeOrdering(faceNodes, n - 1);
            HOTriangle<int> hoTri(triVertId, tmp);

            // Apply reorientation
            if (i == 0)
            {
                // Triangle verts {0,2,1} --> {0,1,2}
                toAlign[0] = 0;
                toAlign[1] = 2;
                toAlign[2] = 1;
                hoTri.Align(toAlign);
            }

            // Fill in mapping.
            for (int j = 0; j < nTriInt; ++j)
            {
                mapping[offset + j] = hoTri.surfVerts[j];
            }
        }
        else
        {
            vector<int> tmp = quadTensorNodeOrdering(faceNodes, n);
            HOQuadrilateral<int> hoQuad(quadVertId, tmp);

            // Apply reorientation
            if (i == 2)
            {
                toAlign[0] = 3;
                toAlign[1] = 2;
                toAlign[2] = 1;
                toAlign[3] = 0;
            }
            else if (i == 3)
            {
                toAlign[0] = 1;
                toAlign[1] = 0;
                toAlign[2] = 3;
                toAlign[3] = 2;
            }
            else if (i == 4)
            {
                toAlign[0] = 3;
                toAlign[1] = 0;
                toAlign[2] = 1;
                toAlign[3] = 2;
            }

            hoQuad.Align(toAlign);

            // Fill in mapping.
            for (int j = 0; j < nQuadInt; ++j)
            {
                mapping[offset + j] = hoQuad.surfVerts[j];
            }
        }
    }

    if (conf.m_volumeNodes == false)
    {
        return mapping;
    }

    // Interior points
    offset = offsets[4] + nQuadInt;
    vector<int> intPoints, tmp;

    for (int i = offset; i < (order+1) * (order+1) * (order+2) / 2; ++i)
    {
        intPoints.push_back(i);
    }

    // Reorder interior points
    tmp = prismTensorNodeOrdering(intPoints, order - 1);
    mapping.insert(mapping.end(), tmp.begin(), tmp.end());

    return mapping;
}

/**
 * @brief Create a reordering for hexahedra.
 */
vector<int> InputGmsh::HexReordering(ElmtConfig conf)
{
    const int order = conf.m_order;
    const int n     = order - 1;
    const int n2    = n * n;
    int i, j, k;

    vector<int> mapping;

    // Map taking Gmsh edges to Nektar++ edges.
    static int gmshToNekEdge[12] = {0, -3, 4, 1, 5, 2, 6, 7, 8, -11, 9, 10};

    // Push back vertices.
    mapping.resize(8);
    for (i = 0; i < 8; ++i)
    {
        mapping[i] = i;
    }

    if (order == 1)
    {
        return mapping;
    }

    // Curvilinear edges
    mapping.resize(8 + 12 * n);

    // Reorder edges.
    int cnt = 8, offset;
    for (i = 0; i < 12; ++i)
    {
        int edge = gmshToNekEdge[i];
        offset   = 8 + n * abs(edge);

        if (edge < 0)
        {
            for (int j = 0; j < n; ++j)
            {
                mapping[offset + n - j - 1] = cnt++;
            }
        }
        else
        {
            for (int j = 0; j < n; ++j)
            {
                mapping[offset + j] = cnt++;
            }
        }
    }

    if (conf.m_faceNodes == false || n2 == 0)
    {
        return mapping;
    }

    // Curvilinear face nodes.
    mapping.resize(8 + 12 * n + 6 * n2);

    // Map which takes Gmsh -> Nektar++ faces in the local element.
    static int gmsh2NekFace[6] = {0, 1, 4, 2, 3, 5};

    // Map which defines orientation between Gmsh and Nektar++ faces.
    StdRegions::Orientation faceOrient[6] = {
        StdRegions::eDir1FwdDir2_Dir2FwdDir1,
        StdRegions::eDir1FwdDir1_Dir2FwdDir2,
        StdRegions::eDir1FwdDir2_Dir2FwdDir1,
        StdRegions::eDir1FwdDir1_Dir2FwdDir2,
        StdRegions::eDir1BwdDir1_Dir2FwdDir2,
        StdRegions::eDir1FwdDir1_Dir2FwdDir2};

    for (i = 0; i < 6; ++i)
    {
        int face    = gmsh2NekFace[i];
        int offset2 = 8 + 12 * n + i * n2;
        offset      = 8 + 12 * n + face * n2;

        // Create a list of interior face nodes for this face only.
        vector<int> faceNodes(n2);
        for (j = 0; j < n2; ++j)
        {
            faceNodes[j] = offset2 + j;
        }

        // Now get the reordering of this face, which puts Gmsh
        // recursive ordering into Nektar++ row-by-row order.
        vector<int> tmp = quadTensorNodeOrdering(faceNodes, n);

        // Finally reorient the face according to the geometry
        // differences.
        if (faceOrient[i] == StdRegions::eDir1FwdDir1_Dir2FwdDir2)
        {
            // Orientation is the same, just copy.
            for (j = 0; j < n2; ++j)
            {
                mapping[offset + j] = tmp[j];
            }
        }
        else if (faceOrient[i] == StdRegions::eDir1FwdDir2_Dir2FwdDir1)
        {
            // Tranposed faces
            for (j = 0; j < n; ++j)
            {
                for (k = 0; k < n; ++k)
                {
                    mapping[offset + j * n + k] = tmp[k * n + j];
                }
            }
        }
        else if (faceOrient[i] == StdRegions::eDir1BwdDir1_Dir2FwdDir2)
        {
            for (j = 0; j < n; ++j)
            {
                for (k = 0; k < n; ++k)
                {
                    mapping[offset + j * n + k] = tmp[j * n + (n - k - 1)];
                }
            }
        }
    }

    if (conf.m_volumeNodes == false)
    {
        return mapping;
    }

    const int totPoints = (order + 1) * (order + 1) * (order + 1);
    vector<int> interior;
    for (i = 8 + 12 * n + 6 * n2; i < totPoints; ++i)
    {
        interior.push_back(i);
    }

    interior = hexTensorNodeOrdering(interior, order - 1);
    mapping.insert(mapping.end(), interior.begin(), interior.end());

    return mapping;
}

/*
 * @brief Populate the element map #elmMap.
 *
 * This function primarily populates the element mapping #elmMap,
 * which takes a msh ID used by Gmsh and translates to element type,
 * element order and whether the element is incomplete (i.e. whether
 * it contains solely boundary nodes, or just face nodes). Note that
 * some of these elements, such as prisms of order >= 3, cannot yet be
 * generated by Gmsh.
 */
std::map<unsigned int, ElmtConfig> InputGmsh::GenElmMap()
{
    using namespace LibUtilities;
    std::map<unsigned int, ElmtConfig> tmp;

    //                    Elmt type,   order,  face, volume
    tmp[  1] = ElmtConfig(eSegment,        1, false, false);
    tmp[  2] = ElmtConfig(eTriangle,       1, false, false);
    tmp[  3] = ElmtConfig(eQuadrilateral,  1, false, false);
    tmp[  4] = ElmtConfig(eTetrahedron,    1, false, false);
    tmp[  5] = ElmtConfig(eHexahedron,     1, false, false);
    tmp[  6] = ElmtConfig(ePrism,          1, false, false);
    tmp[  7] = ElmtConfig(ePyramid,        1, false, false);
    tmp[  8] = ElmtConfig(eSegment,        2,  true, false);
    tmp[  9] = ElmtConfig(eTriangle,       2, false, false);
    tmp[ 10] = ElmtConfig(eQuadrilateral,  2,  true, false);
    tmp[ 11] = ElmtConfig(eTetrahedron,    2, false, false);
    tmp[ 12] = ElmtConfig(eHexahedron,     2,  true,  true);
    tmp[ 13] = ElmtConfig(ePrism,          2,  true, false);
    tmp[ 14] = ElmtConfig(ePyramid,        2,  true, false);
    tmp[ 15] = ElmtConfig(ePoint,          1, false, false);
    tmp[ 16] = ElmtConfig(eQuadrilateral,  2, false, false);
    tmp[ 17] = ElmtConfig(eHexahedron,     2, false, false);
    tmp[ 18] = ElmtConfig(ePrism,          2, false, false);
    tmp[ 19] = ElmtConfig(ePyramid,        2, false, false);
    tmp[ 20] = ElmtConfig(eTriangle,       3, false, false);
    tmp[ 21] = ElmtConfig(eTriangle,       3,  true, false);
    tmp[ 22] = ElmtConfig(eTriangle,       4, false, false);
    tmp[ 23] = ElmtConfig(eTriangle,       4,  true, false);
    tmp[ 24] = ElmtConfig(eTriangle,       5, false, false);
    tmp[ 25] = ElmtConfig(eTriangle,       5,  true, false);
    tmp[ 26] = ElmtConfig(eSegment,        3,  true, false);
    tmp[ 27] = ElmtConfig(eSegment,        4,  true, false);
    tmp[ 28] = ElmtConfig(eSegment,        5,  true, false);
    tmp[ 29] = ElmtConfig(eTetrahedron,    3,  true, false);
    tmp[ 30] = ElmtConfig(eTetrahedron,    4,  true,  true);
    tmp[ 31] = ElmtConfig(eTetrahedron,    5,  true,  true);
    tmp[ 32] = ElmtConfig(eTetrahedron,    4,  true, false);
    tmp[ 33] = ElmtConfig(eTetrahedron,    5,  true, false);
    tmp[ 36] = ElmtConfig(eQuadrilateral,  3,  true, false);
    tmp[ 37] = ElmtConfig(eQuadrilateral,  4,  true, false);
    tmp[ 38] = ElmtConfig(eQuadrilateral,  5,  true, false);
    tmp[ 39] = ElmtConfig(eQuadrilateral,  3, false, false);
    tmp[ 40] = ElmtConfig(eQuadrilateral,  4, false, false);
    tmp[ 41] = ElmtConfig(eQuadrilateral,  5, false, false);
    tmp[ 42] = ElmtConfig(eTriangle,       6,  true, false);
    tmp[ 43] = ElmtConfig(eTriangle,       7,  true, false);
    tmp[ 44] = ElmtConfig(eTriangle,       8,  true, false);
    tmp[ 45] = ElmtConfig(eTriangle,       9,  true, false);
    tmp[ 46] = ElmtConfig(eTriangle,      10,  true, false);
    tmp[ 47] = ElmtConfig(eQuadrilateral,  6,  true, false);
    tmp[ 48] = ElmtConfig(eQuadrilateral,  7,  true, false);
    tmp[ 49] = ElmtConfig(eQuadrilateral,  8,  true, false);
    tmp[ 50] = ElmtConfig(eQuadrilateral,  9,  true, false);
    tmp[ 51] = ElmtConfig(eQuadrilateral, 10,  true, false);
    tmp[ 52] = ElmtConfig(eTriangle,       6, false, false);
    tmp[ 53] = ElmtConfig(eTriangle,       7, false, false);
    tmp[ 54] = ElmtConfig(eTriangle,       8, false, false);
    tmp[ 55] = ElmtConfig(eTriangle,       9, false, false);
    tmp[ 56] = ElmtConfig(eTriangle,      10, false, false);
    tmp[ 57] = ElmtConfig(eQuadrilateral,  6, false, false);
    tmp[ 58] = ElmtConfig(eQuadrilateral,  7, false, false);
    tmp[ 59] = ElmtConfig(eQuadrilateral,  8, false, false);
    tmp[ 60] = ElmtConfig(eQuadrilateral,  9, false, false);
    tmp[ 61] = ElmtConfig(eQuadrilateral, 10, false, false);
    tmp[ 62] = ElmtConfig(eSegment,        6,  true, false);
    tmp[ 63] = ElmtConfig(eSegment,        7,  true, false);
    tmp[ 64] = ElmtConfig(eSegment,        8,  true, false);
    tmp[ 65] = ElmtConfig(eSegment,        9,  true, false);
    tmp[ 66] = ElmtConfig(eSegment,       10,  true, false);
    tmp[ 71] = ElmtConfig(eTetrahedron,    6,  true,  true);
    tmp[ 72] = ElmtConfig(eTetrahedron,    7,  true,  true);
    tmp[ 73] = ElmtConfig(eTetrahedron,    8,  true,  true);
    tmp[ 74] = ElmtConfig(eTetrahedron,    9,  true,  true);
    tmp[ 75] = ElmtConfig(eTetrahedron,   10,  true,  true);
    tmp[ 79] = ElmtConfig(eTetrahedron,    6,  true, false);
    tmp[ 80] = ElmtConfig(eTetrahedron,    7,  true, false);
    tmp[ 81] = ElmtConfig(eTetrahedron,    8,  true, false);
    tmp[ 82] = ElmtConfig(eTetrahedron,    9,  true, false);
    tmp[ 83] = ElmtConfig(eTetrahedron,   10,  true, false);
    tmp[ 90] = ElmtConfig(ePrism,          3,  true,  true);
    tmp[ 91] = ElmtConfig(ePrism,          4,  true,  true);
    tmp[ 92] = ElmtConfig(eHexahedron,     3,  true,  true);
    tmp[ 93] = ElmtConfig(eHexahedron,     4,  true,  true);
    tmp[ 94] = ElmtConfig(eHexahedron,     5,  true,  true);
    tmp[ 95] = ElmtConfig(eHexahedron,     6,  true,  true);
    tmp[ 96] = ElmtConfig(eHexahedron,     7,  true,  true);
    tmp[ 97] = ElmtConfig(eHexahedron,     8,  true,  true);
    tmp[ 98] = ElmtConfig(eHexahedron,     9,  true,  true);
    tmp[ 99] = ElmtConfig(eHexahedron,     3,  true, false);
    tmp[100] = ElmtConfig(eHexahedron,     4,  true, false);
    tmp[101] = ElmtConfig(eHexahedron,     5,  true, false);
    tmp[102] = ElmtConfig(eHexahedron,     6,  true, false);
    tmp[103] = ElmtConfig(eHexahedron,     7,  true, false);
    tmp[104] = ElmtConfig(eHexahedron,     8,  true, false);
    tmp[105] = ElmtConfig(eHexahedron,     9,  true, false);
    tmp[106] = ElmtConfig(ePrism,          5,  true,  true);
    tmp[107] = ElmtConfig(ePrism,          6,  true,  true);
    tmp[108] = ElmtConfig(ePrism,          7,  true,  true);
    tmp[109] = ElmtConfig(ePrism,          8,  true,  true);
    tmp[110] = ElmtConfig(ePrism,          9,  true,  true);
    tmp[111] = ElmtConfig(ePrism,          3,  true, false);
    tmp[112] = ElmtConfig(ePrism,          4,  true, false);
    tmp[113] = ElmtConfig(ePrism,          5,  true, false);
    tmp[114] = ElmtConfig(ePrism,          6,  true, false);
    tmp[115] = ElmtConfig(ePrism,          7,  true, false);
    tmp[116] = ElmtConfig(ePrism,          8,  true, false);
    tmp[117] = ElmtConfig(ePrism,          9,  true, false);
    tmp[118] = ElmtConfig(ePyramid,        3,  true,  true);
    tmp[119] = ElmtConfig(ePyramid,        4,  true,  true);
    tmp[120] = ElmtConfig(ePyramid,        5,  true,  true);
    tmp[121] = ElmtConfig(ePyramid,        6,  true,  true);
    tmp[122] = ElmtConfig(ePyramid,        7,  true,  true);
    tmp[123] = ElmtConfig(ePyramid,        8,  true,  true);
    tmp[124] = ElmtConfig(ePyramid,        9,  true,  true);
    tmp[125] = ElmtConfig(ePyramid,        3,  true, false);
    tmp[126] = ElmtConfig(ePyramid,        4,  true, false);
    tmp[127] = ElmtConfig(ePyramid,        5,  true, false);
    tmp[128] = ElmtConfig(ePyramid,        6,  true, false);
    tmp[129] = ElmtConfig(ePyramid,        7,  true, false);
    tmp[130] = ElmtConfig(ePyramid,        7,  true, false);
    tmp[131] = ElmtConfig(ePyramid,        8,  true, false);

    return tmp;
}

}
}

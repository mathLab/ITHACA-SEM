////////////////////////////////////////////////////////////////////////////////
//
//  File: HOSAlignment.h
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
//  Description: HO aligment routines.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKMESHUTILS_MESHELEMENTS_ALIGNMENT
#define NEKMESHUTILS_MESHELEMENTS_ALIGNMENT

namespace Nektar
{
namespace NekMeshUtils
{
/**
 * @brief A lightweight struct for dealing with high-order triangle
 * alignment.
 *
 * The logic underlying these routines is taken from the original Nektar
 * code.
 */
template <typename T> struct HOTriangle
{
    HOTriangle(std::vector<int> pVertId, std::vector<T> pSurfVerts)
        : vertId(pVertId), surfVerts(pSurfVerts)
    {
    }
    HOTriangle(std::vector<int> pVertId) : vertId(pVertId)
    {
    }

    /// The triangle vertex IDs
    std::vector<int> vertId;

    /// The triangle surface vertices -- templated so that this can
    /// either be nodes or IDs.
    std::vector<T> surfVerts;

    /**
     * @brief Rotates the triangle of data points inside #surfVerts
     * counter-clockwise nrot times.
     *
     * @param nrot Number of times to rotate triangle.
     */
    void Rotate(int nrot)
    {
        int n, i, j, cnt;
        int np = ((int)sqrt(8.0 * surfVerts.size() + 1.0) - 1) / 2;
        std::vector<T> tmp(np * np);

        for (n = 0; n < nrot; ++n)
        {
            for (cnt = i = 0; i < np; ++i)
            {
                for (j = 0; j < np - i; ++j, cnt++)
                {
                    tmp[i * np + j] = surfVerts[cnt];
                }
            }
            for (cnt = i = 0; i < np; ++i)
            {
                for (j = 0; j < np - i; ++j, cnt++)
                {
                    surfVerts[cnt] = tmp[(np - 1 - i - j) * np + i];
                }
            }
        }
    }

    /**
     * @brief Reflect data points inside #surfVerts.
     *
     * This applies a mapping essentially doing the following
     * reordering:
     *
     * 9          9
     * 7 8    ->  8 7
     * 4 5 6      6 5 4
     * 0 1 2 3    3 2 1 0
     */
    void Reflect()
    {
        int i, j, cnt;
        int np = ((int)sqrt(8.0 * surfVerts.size() + 1.0) - 1) / 2;
        std::vector<T> tmp(np * np);

        for (cnt = i = 0; i < np; ++i)
        {
            for (j = 0; j < np - i; ++j, cnt++)
            {
                tmp[i * np + np - i - 1 - j] = surfVerts[cnt];
            }
        }

        for (cnt = i = 0; i < np; ++i)
        {
            for (j = 0; j < np - i; ++j, cnt++)
            {
                surfVerts[cnt] = tmp[i * np + j];
            }
        }
    }

    /**
     * @brief Align this surface to a given vertex ID.
     */
    void Align(std::vector<int> vertId)
    {
        if (vertId[0] == this->vertId[0])
        {
            if (vertId[1] == this->vertId[1] || vertId[1] == this->vertId[2])
            {
                if (vertId[1] == this->vertId[2])
                {
                    Rotate(1);
                    Reflect();
                }
            }
        }
        else if (vertId[0] == this->vertId[1])
        {
            if (vertId[1] == this->vertId[0] || vertId[1] == this->vertId[2])
            {
                if (vertId[1] == this->vertId[0])
                {
                    Reflect();
                }
                else
                {
                    Rotate(2);
                }
            }
        }
        else if (vertId[0] == this->vertId[2])
        {
            if (vertId[1] == this->vertId[0] || vertId[1] == this->vertId[1])
            {
                if (vertId[1] == this->vertId[1])
                {
                    Rotate(2);
                    Reflect();
                }
                else
                {
                    Rotate(1);
                }
            }
        }
    }
};

typedef HOTriangle<NodeSharedPtr> HOSurf;
typedef std::shared_ptr<HOSurf> HOSurfSharedPtr;

/**
 * Hash class for high-order surfaces.
 */
struct HOSurfHash : std::unary_function<HOSurfSharedPtr, std::size_t>
{
    /**
     * Calculate hash of a given high-order surface p by taking
     * successive hashes of the vertex IDs.
     */
    std::size_t operator()(HOSurfSharedPtr const &p) const
    {
        std::vector<int> ids = p->vertId;
        std::sort(ids.begin(), ids.end());
        return hash_range(ids.begin(), ids.end());
    }
};

NEKMESHUTILS_EXPORT bool operator==(HOSurfSharedPtr const &p1,
                                    HOSurfSharedPtr const &p2);

typedef std::unordered_set<HOSurfSharedPtr, HOSurfHash> HOSurfSet;

/**
 * @brief A lightweight struct for dealing with high-order quadrilateral
 * alignment.
 */
template <typename T> struct HOQuadrilateral
{
    HOQuadrilateral(std::vector<int> pVertId, std::vector<T> pSurfVerts)
        : vertId(pVertId), surfVerts(pSurfVerts)
    {
    }

    HOQuadrilateral(std::vector<int> pVertId) : vertId(pVertId)
    {
    }

    /// The quadrilateral vertex IDs
    std::vector<int> vertId;

    /// The quadrilateral surface vertices -- templated so that this can either
    /// be nodes or IDs.
    std::vector<T> surfVerts;

    void ReverseX()
    {
        int np = (int)(sqrt((NekDouble)surfVerts.size()) + 0.5);
        for (int i = 0; i < np; ++i)
        {
            for (int j = 0; j < np/2; ++j)
            {
                swap(surfVerts[i*np + j], surfVerts[i*np + np-j-1]);
            }
        }
    }

    void ReverseY()
    {
        int np = (int)(sqrt((NekDouble)surfVerts.size()) + 0.5);
        // Reverse y direction
        for (int j = 0; j < np; ++j)
        {
            for (int i = 0; i < np/2; ++i)
            {
                swap(surfVerts[i*np + j], surfVerts[(np-i-1)*np + j]);
            }
        }
    }

    void Transpose()
    {
        int np = (int)(sqrt((NekDouble)surfVerts.size()) + 0.5);
        std::vector<T> tmp(surfVerts.size());

        for (int i = 0; i < np; ++i)
        {
            for (int j = 0; j < np; ++j)
            {
                tmp[i*np+j] = surfVerts[j*np+i];
            }
        }

        surfVerts = tmp;
    }

    /**
     * @brief Align this surface to a given vertex ID.
     */
    void Align(std::vector<int> vertId)
    {

        int vmap[4] = {-1, -1, -1, -1};

        // Determine which vertices map to vertId
        for (int i = 0; i < 4; ++i)
        {
            for (int j = 0; j < 4; ++j)
            {
                if (this->vertId[j] == vertId[i])
                {
                    vmap[i] = j;
                    break;
                }
            }

            ASSERTL1(vmap[i] != -1,
                     "Could not determine mapping between vertex IDs");
        }

        StdRegions::Orientation orient = StdRegions::eNoOrientation;

        if (vmap[1] == (vmap[0]+1) % 4)
        {
            switch (vmap[0])
            {
                case 0:
                    orient = StdRegions::eDir1FwdDir1_Dir2FwdDir2;
                    break;
                case 1:
                    orient = StdRegions::eDir1BwdDir2_Dir2FwdDir1;
                    break;
                case 2:
                    orient = StdRegions::eDir1BwdDir1_Dir2BwdDir2;
                    break;
                case 3:
                    orient = StdRegions::eDir1FwdDir2_Dir2BwdDir1;
                    break;
            }
        }
        else
        {
            switch (vmap[0])
            {
                case 0:
                    orient = StdRegions::eDir1FwdDir2_Dir2FwdDir1;
                    break;
                case 1:
                    orient = StdRegions::eDir1BwdDir1_Dir2FwdDir2;
                    break;
                case 2:
                    orient = StdRegions::eDir1BwdDir2_Dir2BwdDir1;
                    break;
                case 3:
                    orient = StdRegions::eDir1FwdDir1_Dir2BwdDir2;
                    break;
            }
        }

        if (orient == StdRegions::eDir1FwdDir2_Dir2FwdDir1 ||
            orient == StdRegions::eDir1BwdDir2_Dir2FwdDir1 ||
            orient == StdRegions::eDir1FwdDir2_Dir2BwdDir1 ||
            orient == StdRegions::eDir1BwdDir2_Dir2BwdDir1)
        {
            Transpose();
        }

        if (orient == StdRegions::eDir1BwdDir1_Dir2FwdDir2 ||
            orient == StdRegions::eDir1BwdDir1_Dir2BwdDir2 ||
            orient == StdRegions::eDir1FwdDir2_Dir2BwdDir1 ||
            orient == StdRegions::eDir1BwdDir2_Dir2BwdDir1)
        {
            ReverseX();
        }

        if (orient == StdRegions::eDir1FwdDir1_Dir2BwdDir2 ||
            orient == StdRegions::eDir1BwdDir1_Dir2BwdDir2 ||
            orient == StdRegions::eDir1BwdDir2_Dir2FwdDir1 ||
            orient == StdRegions::eDir1BwdDir2_Dir2BwdDir1)
        {
            ReverseY();
        }
    }
};

}
}
#endif

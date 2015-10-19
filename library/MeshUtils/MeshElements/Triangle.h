////////////////////////////////////////////////////////////////////////////////
//
//  File: Mesh.h
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
//  Description: Mesh manipulation objects.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef MESHUTILS_MESHELEMENTS_TRIANGLE
#define MESHUTILS_MESHELEMENTS_TRIANGLE

namespace Nektar
{
namespace Utilities
{
    /**
     * @brief A lightweight struct for dealing with high-order triangle
     * alignment.
     *
     * The logic underlying these routines is taken from the original Nektar
     * code.
     */
    template<typename T>
    struct HOTriangle
    {
        HOTriangle(vector<int> pVertId,
                   vector<T>   pSurfVerts) :
            vertId(pVertId), surfVerts(pSurfVerts) {}
        HOTriangle(vector<int> pVertId) : vertId(pVertId) {}

        /// The triangle vertex IDs
        vector<int> vertId;

        /// The triangle surface vertices -- templated so that this can
        /// either be nodes or IDs.
        vector<T> surfVerts;

        /**
         * @brief Rotates the triangle of data points inside #surfVerts
         * counter-clockwise nrot times.
         *
         * @param nrot Number of times to rotate triangle.
         */
        void Rotate(int nrot)
        {
            int n, i, j, cnt;
            int np = ((int)sqrt(8.0*surfVerts.size()+1.0)-1)/2;
            vector<T> tmp(np*np);

            for (n = 0; n < nrot; ++n)
            {
                for (cnt = i = 0; i < np; ++i)
                {
                    for (j = 0; j < np-i; ++j, cnt++)
                    {
                        tmp[i*np+j] = surfVerts[cnt];
                    }
                }
                for (cnt = i = 0; i < np; ++i)
                {
                    for (j = 0; j < np-i; ++j,cnt++)
                    {
                        surfVerts[cnt] = tmp[(np-1-i-j)*np+i];
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
            int np = ((int)sqrt(8.0*surfVerts.size()+1.0)-1)/2;
            vector<T> tmp(np*np);

            for (cnt = i = 0; i < np; ++i)
            {
                for (j = 0; j < np-i; ++j,cnt++)
                {
                    tmp[i*np+np-i-1-j] = surfVerts[cnt];
                }
            }

            for (cnt = i = 0; i < np; ++i)
            {
                for(j = 0; j < np-i; ++j,cnt++)
                {
                    surfVerts[cnt] = tmp[i*np+j];
                }
            }
        }

        /**
         * @brief Align this surface to a given vertex ID.
         */
        void Align(vector<int> vertId)
        {
            if (vertId[0] == this->vertId[0])
            {
                if (vertId[1] == this->vertId[1] ||
                    vertId[1] == this->vertId[2])
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
                if (vertId[1] == this->vertId[0] ||
                    vertId[1] == this->vertId[2])
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
                if (vertId[1] == this->vertId[0] ||
                    vertId[1] == this->vertId[1])
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

    /**
     * @brief A 2-dimensional three-sided element.
     */
    class Triangle : public Element {
    public:
        /// Creates an instance of this class
        static ElementSharedPtr create(
            ElmtConfig                 pConf,
            std::vector<NodeSharedPtr> pNodeList,
            std::vector<int>           pTagList)
        {
            ElementSharedPtr e = boost::shared_ptr<Element>(
                new Triangle(pConf, pNodeList, pTagList));
            vector<EdgeSharedPtr> m_edges = e->GetEdgeList();
            for (int i = 0; i < m_edges.size(); ++i)
            {
                m_edges[i]->m_elLink.push_back(pair<ElementSharedPtr, int>(e,i));
            }
            return e;
        }
        /// Element type
        static LibUtilities::ShapeType m_type;

        Triangle(ElmtConfig                 pConf,
                 std::vector<NodeSharedPtr> pNodeList,
                 std::vector<int>           pTagList);
        Triangle(const Triangle& pSrc);
        virtual ~Triangle() {}

        virtual SpatialDomains::GeometrySharedPtr GetGeom(int coordDim);
        virtual void Complete(int order);

        static unsigned int GetNumNodes(ElmtConfig pConf);
    };


}
}

#endif

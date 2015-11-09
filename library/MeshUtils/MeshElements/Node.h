////////////////////////////////////////////////////////////////////////////////
//
//  File: Node.h
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

#ifndef MESHUTILS_MESHELEMENTS_NODE
#define MESHUTILS_MESHELEMENTS_NODE

namespace Nektar
{
namespace MeshUtils
{
    // Forwards declaration for Element class.
    class Element;
    /// Shared pointer to an element.
    typedef boost::shared_ptr<Element> ElementSharedPtr;

    class Node;
    typedef boost::shared_ptr<Node> NodeSharedPtr;

    /**
     * @brief Represents a point in the domain.
     *
     * Such points may either be element vertices, or simply control
     * points on high-order edges/faces, although this information is not
     * contained within this class.
     */
    class Node {
    public:
        /// Create a new node at a specified coordinate.
        Node(int pId, NekDouble pX, NekDouble pY, NekDouble pZ)
            : m_id(pId), m_x(pX), m_y(pY), m_z(pZ), m_geom() {}
        /// Copy an existing node.
        //Node(const Node& pSrc)
        //    : m_id(pSrc.m_id), m_x(pSrc.m_x), m_y(pSrc.m_y),
        //      m_z(pSrc.m_z), m_geom() {}
        Node() : m_id(0), m_x(0.0), m_y(0.0), m_z(0.0), m_geom() {}
        ~Node() {}

        /// Reset the local id;
        void SetID(int pId)
        {
            m_id = pId;
        }

        /// Get the local id;
        int GetID(void)
        {
            return m_id;
        }

        /// Define node ordering based on ID.
        bool operator<(const Node& pSrc)
        {
            return (m_id < pSrc.m_id);
        }
        /// Define node equality based on coordinate.
        bool operator==(const Node& pSrc)
        {
            return m_x == pSrc.m_x && m_y == pSrc.m_y && m_z == pSrc.m_z;
        }

        Node operator+(const Node &pSrc) const
        {
            return Node(m_id, m_x+pSrc.m_x, m_y+pSrc.m_y, m_z+pSrc.m_z);
        }

        Node operator-(const Node &pSrc) const
        {
            return Node(m_id, m_x-pSrc.m_x, m_y-pSrc.m_y, m_z-pSrc.m_z);
        }

        Node operator*(const Node &pSrc) const
        {
            return Node(m_id, m_x*pSrc.m_x, m_y*pSrc.m_y, m_z*pSrc.m_z);
        }

        Node operator*(const NekDouble &alpha) const
        {
            return Node(m_id, alpha*m_x, alpha*m_y, alpha*m_z);
        }

        Node operator/(const NekDouble &alpha) const
        {
            return Node(m_id, m_x/alpha, m_y/alpha, m_z/alpha);
        }

        void operator+=(const Node &pSrc)
        {
            m_x += pSrc.m_x;
            m_y += pSrc.m_y;
            m_z += pSrc.m_z;
        }

        void operator*=(const NekDouble &alpha)
        {
            m_x *= alpha;
            m_y *= alpha;
            m_z *= alpha;
        }

        void operator/=(const NekDouble &alpha)
        {
            m_x /= alpha;
            m_y /= alpha;
            m_z /= alpha;
        }

        NekDouble abs2() const
        {
            return m_x*m_x+m_y*m_y+m_z*m_z;
        }

        NekDouble dot(const Node &pSrc) const
        {
            return m_x*pSrc.m_x + m_y*pSrc.m_y + m_z*pSrc.m_z;
        }


        Node curl(const Node &pSrc) const
        {
            return Node(m_id, m_y*pSrc.m_z - m_z*pSrc.m_y,
                        m_z*pSrc.m_x-m_x*pSrc.m_z, m_x*pSrc.m_y-m_y*pSrc.m_x);
        }

        void SetCADCurve(int i, NekDouble t)
        {
            CADCurve[i] = t;
        }

        void SetCADSurf(int i, Array<OneD, NekDouble> uv)
        {
            CADSurf[i] = uv;
        }

        Array<OneD, NekDouble> GetLoc()
        {
            Array<OneD, NekDouble> out(3);
            out[0] = m_x; out[1] = m_y; out[2] = m_z;
            return out;
        }

        NekDouble GetCADCurve(int i)
        {
            std::map<int, NekDouble>::iterator search =
                            CADCurve.find(i);
            ASSERTL0(search->first == i, "node not on this curve");

            return search->second;
        }

        Array<OneD, NekDouble> GetCADSurf(int i)
        {
            //I dont know why I ahev to do this to get it to work
            //this really needs bound checking
            std::map<int, Array<OneD, NekDouble> >::iterator search =
                        CADSurf.find(i);
            ASSERTL0(search->first == i,"surface not found");

            return search->second;
        }

        std::vector<int> GetListCADCurve()
        {
            std::vector<int> list;
            std::map<int, NekDouble >::iterator c;
            for(c = CADCurve.begin(); c != CADCurve.end(); c++)
            {
                list.push_back(c->first);
            }
            return list;
        }

        std::vector<int> GetListCADSurf()
        {
            std::vector<int> list;
            std::map<int, Array<OneD, NekDouble> >::iterator s;
            for(s = CADSurf.begin(); s != CADSurf.end(); s++)
            {
                list.push_back(s->first);
            }
            return list;
        }

        NekDouble Distance(NodeSharedPtr &p)
        {
            return sqrt((m_x-p->m_x)*(m_x-p->m_x)+(m_y-p->m_y)*(m_y-p->m_y)
                        +(m_z-p->m_z)*(m_z-p->m_z));
        }

        NekDouble Angle(NodeSharedPtr &a, NodeSharedPtr &b)
        {
            Array<OneD,NekDouble> va(3),vb(3),cn(3);
            va[0] = a->m_x - m_x;
            va[1] = a->m_y - m_y;
            va[2] = a->m_z - m_z;
            vb[0] = b->m_x - m_x;
            vb[1] = b->m_y - m_y;
            vb[2] = b->m_z - m_z;

            NekDouble lva = sqrt(va[0]*va[0] + va[1]*va[1] + va[2]*va[2]);
            NekDouble lvb = sqrt(vb[0]*vb[0] + vb[1]*vb[1] + vb[2]*vb[2]);

            NekDouble aw = 1.0/(lva*lvb);

            NekDouble cosw = (va[0]*vb[0] + va[1]*vb[1] + va[2]*vb[2])*aw;

            cn[0] = vb[1]*va[2] - vb[2]*va[1];
            cn[1] = vb[2]*va[0] - vb[0]*va[2];
            cn[2] = vb[0]*va[1] - vb[1]*va[0];

            NekDouble lcn = sqrt(cn[0]*cn[0] + cn[1]*cn[1] + cn[2]*cn[2]);

            NekDouble sinw = aw*lcn;

            NekDouble an = atan2(sinw,cosw);

            if(an < 0)
                an += 6.2831853071796;

            return an;
        }

        void Move(Array<OneD, NekDouble> l, int s, Array<OneD, NekDouble> uv)
        {
            m_x = l[0]; m_y = l[1]; m_z = l[2];
            CADSurf[s] = uv;
        }

        /// Generate a %SpatialDomains::PointGeom for this node.
        SpatialDomains::PointGeomSharedPtr GetGeom(int coordDim)
        {
            SpatialDomains::PointGeomSharedPtr ret =
                MemoryManager<SpatialDomains::PointGeom>
                    ::AllocateSharedPtr(coordDim,m_id,m_x,m_y,m_z);

            return ret;
        }

        /// ID of node.
        int m_id;
        /// X-coordinate.
        NekDouble m_x;
        /// Y-coordinate.
        NekDouble m_y;
        /// Z-coordinate.
        NekDouble m_z;
        ///list of cadcurves the node lies on
        std::map<int, NekDouble> CADCurve;
        ///list of cadsurfs the node lies on
        std::map<int, Array<OneD, NekDouble> > CADSurf;
        ///edge link count;
        int m_elCount;
        ///this could probably do with a better method but its easy for now
        vector<NodeSharedPtr> m_connectingNodes;

    private:
        SpatialDomains::PointGeomSharedPtr m_geom;
    };
    /// Shared pointer to a Node.

    bool operator==(NodeSharedPtr const &p1, NodeSharedPtr const &p2);
    bool operator< (NodeSharedPtr const &p1, NodeSharedPtr const &p2);
    bool operator!=(NodeSharedPtr const &p1, NodeSharedPtr const &p2);
    std::ostream &operator<<(std::ostream &os, const NodeSharedPtr &n);

    /**
     * @brief Defines a hash function for nodes.
     *
     * The hash of a node is straight-forward; a combination of the x, y,
     * and z co-ordinates in this order.
     */
    struct NodeHash : std::unary_function<NodeSharedPtr, std::size_t>
    {
        std::size_t operator()(NodeSharedPtr const& p) const
        {
            std::size_t seed = 0;
            boost::hash_combine(seed, p->m_x);
            boost::hash_combine(seed, p->m_y);
            boost::hash_combine(seed, p->m_z);
            return seed;
        }
    };
    typedef boost::unordered_set<NodeSharedPtr, NodeHash> NodeSet;
}
}

#endif

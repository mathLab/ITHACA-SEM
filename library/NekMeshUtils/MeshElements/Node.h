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
//  Description: Mesh node object.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKMESHUTILS_MESHELEMENTS_NODE
#define NEKMESHUTILS_MESHELEMENTS_NODE

#include <LibUtilities/BasicUtils/HashUtils.hpp>
#include <NekMeshUtils/NekMeshUtilsDeclspec.h>

#include <iomanip>

#include <NekMeshUtils/CADSystem/CADCurve.h>
#include <NekMeshUtils/CADSystem/CADSurf.h>
#include <NekMeshUtils/CADSystem/CADSystem.h>
#include <SpatialDomains/PointGeom.h>

namespace Nektar
{
namespace NekMeshUtils
{
class Node;
typedef std::shared_ptr<Node> NodeSharedPtr;

/**
 * @brief Represents a point in the domain.
 *
 * Such points may either be element vertices, or simply control
 * points on high-order edges/faces, although this information is not
 * contained within this class.
 */
class Node
{
public:
    /// Create a new node at a specified coordinate.
    NEKMESHUTILS_EXPORT Node(int pId, NekDouble pX, NekDouble pY, NekDouble pZ)
        : m_id(pId), m_x(pX), m_y(pY), m_z(pZ), m_geom()
    {
    }
    /// Copy an existing node.
    // Node(const Node& pSrc)
    //    : m_id(pSrc.m_id), m_x(pSrc.m_x), m_y(pSrc.m_y),
    //      m_z(pSrc.m_z), m_geom() {}
    /// create an empty node
    NEKMESHUTILS_EXPORT Node() : m_id(0), m_x(0.0), m_y(0.0), m_z(0.0), m_geom()
    {
    }
    NEKMESHUTILS_EXPORT ~Node()
    {
    }

    /// Reset the local id;
    NEKMESHUTILS_EXPORT void SetID(int pId)
    {
        m_id = pId;
    }

    /// Get the local id;
    NEKMESHUTILS_EXPORT int GetID(void)
    {
        return m_id;
    }

    /// Define node ordering based on ID.
    NEKMESHUTILS_EXPORT bool operator<(const Node &pSrc)
    {
        return (m_id < pSrc.m_id);
    }
    /// Define node equality based on coordinate.
    NEKMESHUTILS_EXPORT bool operator==(const Node &pSrc)
    {
        return LibUtilities::IsRealEqual(m_x, pSrc.m_x) &&
               LibUtilities::IsRealEqual(m_y, pSrc.m_y) &&
               LibUtilities::IsRealEqual(m_z, pSrc.m_z);
    }

    NEKMESHUTILS_EXPORT Node operator+(const Node &pSrc) const
    {
        return Node(m_id, m_x + pSrc.m_x, m_y + pSrc.m_y, m_z + pSrc.m_z);
    }

    NEKMESHUTILS_EXPORT Node operator-(const Node &pSrc) const
    {
        return Node(m_id, m_x - pSrc.m_x, m_y - pSrc.m_y, m_z - pSrc.m_z);
    }

    NEKMESHUTILS_EXPORT Node operator*(const Node &pSrc) const
    {
        return Node(m_id, m_x * pSrc.m_x, m_y * pSrc.m_y, m_z * pSrc.m_z);
    }

    NEKMESHUTILS_EXPORT Node operator*(const NekDouble &alpha) const
    {
        return Node(m_id, alpha * m_x, alpha * m_y, alpha * m_z);
    }

    NEKMESHUTILS_EXPORT Node operator/(const NekDouble &alpha) const
    {
        return Node(m_id, m_x / alpha, m_y / alpha, m_z / alpha);
    }

    NEKMESHUTILS_EXPORT void operator+=(const Node &pSrc)
    {
        m_x += pSrc.m_x;
        m_y += pSrc.m_y;
        m_z += pSrc.m_z;
    }

    NEKMESHUTILS_EXPORT void operator*=(const NekDouble &alpha)
    {
        m_x *= alpha;
        m_y *= alpha;
        m_z *= alpha;
    }

    NEKMESHUTILS_EXPORT void operator/=(const NekDouble &alpha)
    {
        m_x /= alpha;
        m_y /= alpha;
        m_z /= alpha;
    }

    NEKMESHUTILS_EXPORT NodeSharedPtr copy()
    {
        return std::shared_ptr<Node>(new Node(m_id, m_x, m_y, m_z));
    }

    NEKMESHUTILS_EXPORT NekDouble abs2() const
    {
        return m_x * m_x + m_y * m_y + m_z * m_z;
    }

    NEKMESHUTILS_EXPORT NekDouble dot(const Node &pSrc) const
    {
        return m_x * pSrc.m_x + m_y * pSrc.m_y + m_z * pSrc.m_z;
    }

    NEKMESHUTILS_EXPORT Node curl(const Node &pSrc) const
    {
        return Node(m_id, m_y * pSrc.m_z - m_z * pSrc.m_y,
                    m_z * pSrc.m_x - m_x * pSrc.m_z,
                    m_x * pSrc.m_y - m_y * pSrc.m_x);
    }

    NEKMESHUTILS_EXPORT Array<OneD, NekDouble> GetLoc()
    {
        Array<OneD, NekDouble> out(3);
        out[0] = m_x;
        out[1] = m_y;
        out[2] = m_z;
        return out;
    }

    /// Generate a %SpatialDomains::PointGeom for this node.
    NEKMESHUTILS_EXPORT SpatialDomains::PointGeomSharedPtr GetGeom(int coordDim)
    {
        SpatialDomains::PointGeomSharedPtr ret =
            MemoryManager<SpatialDomains::PointGeom>::AllocateSharedPtr(
                coordDim, m_id, m_x, m_y, m_z);

        return ret;
    }

    NEKMESHUTILS_EXPORT NekDouble Distance(NodeSharedPtr &p)
    {
        return sqrt((m_x - p->m_x) * (m_x - p->m_x) +
                    (m_y - p->m_y) * (m_y - p->m_y) +
                    (m_z - p->m_z) * (m_z - p->m_z));
    }

    NEKMESHUTILS_EXPORT NekDouble Angle(NodeSharedPtr &a, NodeSharedPtr &b)
    {
        Array<OneD, NekDouble> va(3), vb(3), cn(3);
        va[0] = a->m_x - m_x;
        va[1] = a->m_y - m_y;
        va[2] = a->m_z - m_z;
        vb[0] = b->m_x - m_x;
        vb[1] = b->m_y - m_y;
        vb[2] = b->m_z - m_z;

        NekDouble lva = sqrt(va[0] * va[0] + va[1] * va[1] + va[2] * va[2]);
        NekDouble lvb = sqrt(vb[0] * vb[0] + vb[1] * vb[1] + vb[2] * vb[2]);

        NekDouble aw = 1.0 / (lva * lvb);

        NekDouble cosw = (va[0] * vb[0] + va[1] * vb[1] + va[2] * vb[2]) * aw;

        cn[0] = vb[1] * va[2] - vb[2] * va[1];
        cn[1] = vb[2] * va[0] - vb[0] * va[2];
        cn[2] = vb[0] * va[1] - vb[1] * va[0];

        NekDouble lcn = sqrt(cn[0] * cn[0] + cn[1] * cn[1] + cn[2] * cn[2]);

        NekDouble sinw = aw * lcn;

        NekDouble an = atan2(sinw, cosw);

        if (an < 0)
            an += 6.2831853071796;

        return an;
    }

    // functions for cad information

    void SetCADCurve(CADCurveSharedPtr c, NekDouble t)
    {
        auto it = CADCurveList.find(c->GetId());
        if (it != CADCurveList.end())
        {
            // already in list so remove it
            CADCurveList.erase(it);
        }
        CADCurveList.insert(make_pair(c->GetId(), make_pair(c, t)));
    }

    void SetCADSurf(CADSurfSharedPtr s, Array<OneD, NekDouble> uv)
    {
        auto it = CADSurfList.find(s->GetId());
        if (it != CADSurfList.end())
        {
            // already in list so remove it
            CADSurfList.erase(it);
        }
        CADSurfList.insert(make_pair(s->GetId(), make_pair(s, uv)));
    }

    NekDouble GetCADCurveInfo(int i)
    {
        auto search = CADCurveList.find(i);
        ASSERTL0(search != CADCurveList.end(), "node not on this curve");

        return search->second.second;
    }

    Array<OneD, NekDouble> GetCADSurfInfo(int i)
    {
        auto search = CADSurfList.find(i);
        ASSERTL0(search != CADSurfList.end(), "surface not found");

        return search->second.second;
    }

    std::vector<CADCurveSharedPtr> GetCADCurves()
    {
        std::vector<CADCurveSharedPtr> lst;
        for (auto &c : CADCurveList)
        {
            lst.push_back(c.second.first.lock());
        }
        return lst;
    }

    std::vector<CADSurfSharedPtr> GetCADSurfs()
    {
        std::vector<CADSurfSharedPtr> lst;
        for (auto &s : CADSurfList)
        {
            lst.push_back(s.second.first.lock());
        }
        return lst;
    }

    int GetNumCadCurve()
    {
        return CADCurveList.size();
    }

    int GetNumCADSurf()
    {
        return CADSurfList.size();
    }

    void Move(Array<OneD, NekDouble> l, int s, Array<OneD, NekDouble> uv)
    {
        m_x                 = l[0];
        m_y                 = l[1];
        m_z                 = l[2];
        CADSurfSharedPtr su = CADSurfList[s].first.lock();
        SetCADSurf(su, uv);
    }

    void Move(NekDouble x, NekDouble y, NekDouble z, int s,
              Array<OneD, NekDouble> uv)
    {
        m_x                 = x;
        m_y                 = y;
        m_z                 = z;
        CADSurfSharedPtr su = CADSurfList[s].first.lock();
        SetCADSurf(su, uv);
    }

    void Move(Array<OneD, NekDouble> l, int c, NekDouble t)
    {
        m_x                  = l[0];
        m_y                  = l[1];
        m_z                  = l[2];
        CADCurveSharedPtr cu = CADCurveList[c].first.lock();
        SetCADCurve(cu, t);
    }

    void Move(NekDouble x, NekDouble y, NekDouble z, int c, NekDouble t)
    {
        m_x                  = x;
        m_y                  = y;
        m_z                  = z;
        CADCurveSharedPtr cu = CADCurveList[c].first.lock();
        SetCADCurve(cu, t);
    }

    void Rotate(std::string dir, NekDouble angle)
    {
        if (dir == "x")
        {
            NekDouble yrot = cos(angle) * m_y - sin(angle) * m_z;
            NekDouble zrot = sin(angle) * m_y + cos(angle) * m_z;

            m_y = yrot;
            m_z = zrot;
        }
        else if (dir == "y")
        {
            NekDouble zrot = cos(angle) * m_z - sin(angle) * m_x;
            NekDouble xrot = sin(angle) * m_z + cos(angle) * m_x;

            m_z = zrot;
            m_x = xrot;
        }
        else if (dir == "z")
        {
            NekDouble xrot = cos(angle) * m_x - sin(angle) * m_y;
            NekDouble yrot = sin(angle) * m_x + cos(angle) * m_y;

            m_x = xrot;
            m_y = yrot;
        }
        else
        {
            ASSERTL0(false, "Unrecognised rotational direction: " + dir);
        }
    }

    NekDouble Angle(Array<OneD, NekDouble> locA, Array<OneD, NekDouble> locB,
                    Array<OneD, NekDouble> N)
    {
        // calculates the angle between this node to a to this node to b
        // Uses the CAD surface to orientate the angle
        Array<OneD, NekDouble> A(3), B(3), CP(3);
        A[0] = locA[0] - m_x;
        A[1] = locA[1] - m_y;
        A[2] = locA[2] - m_z;
        B[0] = locB[0] - m_x;
        B[1] = locB[1] - m_y;
        B[2] = locB[2] - m_z;

        CP[0] = A[1] * B[2] - A[2] * B[1];
        CP[1] = -1.0 * (A[0] * B[2] - A[2] * B[0]);
        CP[2] = A[0] * B[1] - A[1] * B[0];

        NekDouble ang = sqrt(CP[0] * CP[0] + CP[1] * CP[1] + CP[2] * CP[2]);

        ang /= sqrt(A[0] * A[0] + A[1] * A[1] + A[2] * A[2]);
        ang /= sqrt(B[0] * B[0] + B[1] * B[1] + B[2] * B[2]);

        NekDouble dot = N[0] * CP[0] + N[1] * CP[1] + N[2] * CP[2];

        ang = asin(ang);

        if (dot < 0.0)
        {
            ang = 2.0 * M_PI - ang;
        }

        return ang;
    }

    /// ID of node.
    int m_id;
    /// X-coordinate.
    NekDouble m_x;
    /// Y-coordinate.
    NekDouble m_y;
    /// Z-coordinate.
    NekDouble m_z;

    /// list of cadcurves the node lies on
    std::map<int, std::pair<std::weak_ptr<CADCurve>, NekDouble>> CADCurveList;
    /// list of cadsurfs the node lies on
    std::map<int, std::pair<std::weak_ptr<CADSurf>, Array<OneD, NekDouble>>>
        CADSurfList;

private:
    SpatialDomains::PointGeomSharedPtr m_geom;
};
/// Shared pointer to a Node.

NEKMESHUTILS_EXPORT bool operator==(NodeSharedPtr const &p1,
                                    NodeSharedPtr const &p2);
NEKMESHUTILS_EXPORT bool operator<(NodeSharedPtr const &p1,
                                   NodeSharedPtr const &p2);
NEKMESHUTILS_EXPORT bool operator!=(NodeSharedPtr const &p1,
                                    NodeSharedPtr const &p2);
NEKMESHUTILS_EXPORT std::ostream &operator<<(std::ostream &os,
                                             const NodeSharedPtr &n);

/// Define node equality based on coordinate with optional custom tolerance factor.
NEKMESHUTILS_EXPORT bool IsNodeEqual(const Node &n1, const Node &n2,
    const unsigned int fact = NekConstants::kNekFloatCompFact);

/**
 * @brief Defines a hash function for nodes.
 *
 * The hash of a node is straight-forward; a combination of the x, y,
 * and z co-ordinates in this order.
 */
struct NodeHash : std::unary_function<NodeSharedPtr, std::size_t>
{
    std::size_t operator()(NodeSharedPtr const &p) const
    {
        return hash_combine(p->m_x, p->m_y, p->m_z);
    }
};
typedef std::unordered_set<NodeSharedPtr, NodeHash> NodeSet;
}
}

#endif

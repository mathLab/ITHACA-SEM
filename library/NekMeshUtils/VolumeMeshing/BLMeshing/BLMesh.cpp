////////////////////////////////////////////////////////////////////////////////
//
//  File: TetMesh.cpp
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
//  Description: tet meshing methods
//
////////////////////////////////////////////////////////////////////////////////

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>

#include <LibUtilities/BasicUtils/Progressbar.hpp>

#include <NekMeshUtils/CADSystem/CADSurf.h>
#include <NekMeshUtils/VolumeMeshing/BLMeshing/BLMesh.h>

#include <LibUtilities/Foundations/ManagerAccess.h>
#include <LibUtilities/Foundations/NodalUtil.h>

#include <algorithm>

namespace bg  = boost::geometry;
namespace bgi = boost::geometry::index;

typedef bg::model::point<double, 3, bg::cs::cartesian> point;
typedef bg::model::box<point> box;
typedef std::pair<box, unsigned int> boxI;

using namespace std;
namespace Nektar
{
namespace NekMeshUtils
{

inline box GetBox(ElementSharedPtr el, NekDouble ov)
{
    NekDouble xmin = numeric_limits<double>::max(),
              xmax = -1.0 * numeric_limits<double>::max(),
              ymin = numeric_limits<double>::max(),
              ymax = -1.0 * numeric_limits<double>::max(),
              zmin = numeric_limits<double>::max(),
              zmax = -1.0 * numeric_limits<double>::max();

    vector<NodeSharedPtr> ns = el->GetVertexList();
    for (int i = 0; i < ns.size(); i++)
    {
        xmin = min(xmin, ns[i]->m_x);
        xmax = max(xmax, ns[i]->m_x);
        ymin = min(ymin, ns[i]->m_y);
        ymax = max(ymax, ns[i]->m_y);
        zmin = min(zmin, ns[i]->m_z);
        zmax = max(zmax, ns[i]->m_z);
    }

    return box(point(xmin - ov, ymin - ov, zmin - ov),
               point(xmax + ov, ymax + ov, zmax + ov));
}

inline box GetBox(vector<ElementSharedPtr> els, NekDouble ov)
{
    NekDouble xmin = numeric_limits<double>::max(),
              xmax = -1.0 * numeric_limits<double>::max(),
              ymin = numeric_limits<double>::max(),
              ymax = -1.0 * numeric_limits<double>::max(),
              zmin = numeric_limits<double>::max(),
              zmax = -1.0 * numeric_limits<double>::max();

    for (int j = 0; j < els.size(); j++)
    {
        vector<NodeSharedPtr> ns = els[j]->GetVertexList();
        for (int i = 0; i < ns.size(); i++)
        {
            xmin = min(xmin, ns[i]->m_x);
            xmax = max(xmax, ns[i]->m_x);
            ymin = min(ymin, ns[i]->m_y);
            ymax = max(ymax, ns[i]->m_y);
            zmin = min(zmin, ns[i]->m_z);
            zmax = max(zmax, ns[i]->m_z);
        }
    }

    return box(point(xmin - ov, ymin - ov, zmin - ov),
               point(xmax + ov, ymax + ov, zmax + ov));
}

inline box GetBox(NodeSharedPtr n, NekDouble ov)
{
    return box(point(n->m_x - ov, n->m_y - ov, n->m_z - ov),
               point(n->m_x + ov, n->m_y + ov, n->m_z + ov));
}

void BLMesh::Mesh()
{
    Setup();

    BuildElements();

    GrowLayers();

    Shrink();

    map<NodeSharedPtr, blInfoSharedPtr>::iterator bit;
    for (bit = m_blData.begin(); bit != m_blData.end(); bit++)
    {
        vector<blInfoSharedPtr> infos = m_nToNInfo[bit->first];
        for (int i = 0; i < infos.size(); i++)
        {
            if (bit->second->bl > infos[i]->bl + 1)
            {
                cout << "non smooth error " << bit->second->bl << " "
                     << infos[i]->bl << endl;
            }
        }
    }

    for (int i = 0; i < m_mesh->m_element[3].size(); i++)
    {
        ElementSharedPtr el = m_mesh->m_element[3][i];
        if (!IsPrismValid(el))
        {
            cout << "validity error " << el->GetId() << endl;
        }
    }

    /*m_mesh->m_element[2] = m_psuedoSurface;
    m_mesh->m_element[3].clear();
    m_mesh->m_expDim = 2;*/
}

map<NodeSharedPtr, NodeSharedPtr> BLMesh::GetSymNodes()
{
    map<NodeSharedPtr, NodeSharedPtr> ret;

    map<NodeSharedPtr, blInfoSharedPtr>::iterator bit;
    for (bit = m_blData.begin(); bit != m_blData.end(); bit++)
    {
        if (!bit->second->onSym)
        {
            continue;
        }
        CADSurfSharedPtr s = m_mesh->m_cad->GetSurf(bit->second->symsurf);
        Array<OneD, NekDouble> loc = bit->second->pNode->GetLoc();
        Array<OneD, NekDouble> uv  = s->locuv(loc);
        bit->second->pNode->SetCADSurf(s, uv);
        ret[bit->first] = bit->second->pNode;
    }
    return ret;
}

inline bool Infont(NodeSharedPtr n, ElementSharedPtr el)
{
    vector<NodeSharedPtr> ns1 = el->GetVertexList();
    Array<OneD, NekDouble> N1 = el->Normal(true);

    Array<OneD, NekDouble> V(3);
    V[0] = n->m_x - ns1[0]->m_x;
    V[1] = n->m_y - ns1[0]->m_y;
    V[2] = n->m_z - ns1[0]->m_z;

    NekDouble Vmag = sqrt(V[0] * V[0] + V[1] * V[1] + V[2] * V[2]);

    NekDouble ang = (N1[0] * V[0] + N1[1] * V[1] + N1[2] * V[2]) / Vmag;

    return ang > 0.17;
}

void BLMesh::GrowLayers()
{
    map<NodeSharedPtr, blInfoSharedPtr>::iterator bit;

    // setup up a tree which is formed of boxes of each surface plus some
    // extra room (ideal bl thick)

    // in each iteration a tree is made for all the triangles in each surface
    // when considering to stop a bounary layer growing, it first
    // looks at the top tree to find surfaces which are canditates
    // it then searches the subtrees for each triangle canditate
    // it then does a distance calcation.
    // if a boundary layer should be close to that from another surface, it
    // should stop

    map<int, vector<ElementSharedPtr>> psElements;
    for (int i = 0; i < m_mesh->m_element[2].size(); i++)
    {
        ElementSharedPtr el = m_mesh->m_element[2][i];
        vector<unsigned int>::iterator f =
            find(m_blsurfs.begin(), m_blsurfs.end(), el->m_parentCAD->GetId());

        vector<unsigned int>::iterator s = find(
            m_symSurfs.begin(), m_symSurfs.end(), el->m_parentCAD->GetId());

        if (f == m_blsurfs.end() && s == m_symSurfs.end())
        {
            psElements[el->m_parentCAD->GetId()].push_back(el);
        }
    }
    for (int i = 0; i < m_psuedoSurface.size(); i++)
    {
        psElements[m_psuedoSurface[i]->m_parentCAD->GetId()].push_back(
            m_psuedoSurface[i]);
    }

    bgi::rtree<boxI, bgi::quadratic<16>> TopTree;
    map<int, bgi::rtree<boxI, bgi::quadratic<16>>> SubTrees;

    // ofstream file;
    // file.open("pts.3D");
    // file << "x y z value" << endl;

    for (int l = 1; l < m_layer; l++)
    {
        NekDouble delta = (m_layerT[l] - m_layerT[l - 1]);
        TopTree.clear();
        SubTrees.clear();
        map<int, vector<ElementSharedPtr>>::iterator it;
        for (it = psElements.begin(); it != psElements.end(); it++)
        {
            TopTree.insert(make_pair(GetBox(it->second, m_bl), it->first));
            vector<boxI> toInsert;
            for (int i = 0; i < it->second.size(); i++)
            {
                toInsert.push_back(make_pair(GetBox(it->second[i], m_bl), i));
            }
            SubTrees[it->first].insert(toInsert.begin(), toInsert.end());
        }

        for (bit = m_blData.begin(); bit != m_blData.end(); bit++)
        {
            if (bit->second->stopped)
            {
                continue;
            }

            vector<boxI> results;
            TopTree.query(bgi::intersects(point(bit->second->pNode->m_x,
                                                bit->second->pNode->m_y,
                                                bit->second->pNode->m_z)),
                          back_inserter(results));
            set<int> surfs;
            for (int i = 0; i < results.size(); i++)
            {
                set<int>::iterator f =
                    bit->second->surfs.find(results[i].second);
                if (f == bit->second->surfs.end())
                {
                    // hit
                    surfs.insert(results[i].second);
                }
            }

            set<int>::iterator iit;
            bool hit = false;
            for (iit = surfs.begin(); iit != surfs.end(); iit++)
            {
                results.clear();
                SubTrees[*iit].query(
                    bgi::intersects(GetBox(bit->second->pNode, m_bl)),
                    back_inserter(results));
                for (int i = 0; i < results.size(); i++)
                {
                    if (Infont(bit->second->pNode,
                               psElements[*iit][results[i].second]))
                    {
                        NekDouble prox =
                            Proximity(bit->second->pNode,
                                      psElements[*iit][results[i].second]);
                        if (prox < delta * 2.5)
                        {
                            hit = true;
                            // cout << "hit" << endl;
                            bit->second->stopped = true;
                            /*file << bit->first->m_x << " " << bit->first->m_y
                            << " " << bit->first->m_z << " " << l << endl;
                            file << bit->second->pNode->m_x << " " <<
                            bit->second->pNode->m_y << " " <<
                            bit->second->pNode->m_z << " " << l << endl;
                            m_mesh->m_element[2].clear();
                            m_mesh->m_expDim--;
                            m_mesh->m_element[2].push_back(psElements[*iit][results[i].second]);*/
                            break;
                            // return;
                        }
                    }
                }
                if (hit)
                    break;
            }
        }

        // if after proximity scanning all is okay, advance the layer
        for (bit = m_blData.begin(); bit != m_blData.end(); bit++)
        {
            if (bit->second->stopped)
            {
                continue;
            }

            // test the smoothness
            bool shouldStop            = false;
            vector<blInfoSharedPtr> ne = m_nToNInfo[bit->first];
            for (int i = 0; i < ne.size(); i++)
            {
                if (ne[i]->bl < bit->second->bl)
                {
                    shouldStop = true;
                    break;
                }
            }
            if (shouldStop)
            {
                bit->second->stopped = true;
                continue;
            }

            bit->second->AlignNode(m_layerT[l]);
            bit->second->bl = l;
        }
    }
    // file.close();
}

inline bool sign(NekDouble a, NekDouble b)
{
    return (a * b > 0.0);
}

inline NekDouble Dot(Array<OneD, NekDouble> a, Array<OneD, NekDouble> b)
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

NekDouble BLMesh::Proximity(NodeSharedPtr n, ElementSharedPtr el)
{
    vector<NodeSharedPtr> ns = el->GetVertexList();
    Array<OneD, NekDouble> B = ns[0]->GetLoc();
    Array<OneD, NekDouble> E0(3);
    E0[0] = ns[1]->m_x - ns[0]->m_x;
    E0[1] = ns[1]->m_y - ns[0]->m_y;
    E0[2] = ns[1]->m_z - ns[0]->m_z;
    Array<OneD, NekDouble> E1(3);
    E1[0]                    = ns[2]->m_x - ns[0]->m_x;
    E1[1]                    = ns[2]->m_y - ns[0]->m_y;
    E1[2]                    = ns[2]->m_z - ns[0]->m_z;
    Array<OneD, NekDouble> P = n->GetLoc();

    NekDouble a = Dot(E0, E0);
    NekDouble b = Dot(E0, E1);
    NekDouble c = Dot(E1, E1);

    Array<OneD, NekDouble> BP(3);
    Vmath::Vsub(3, B, 1, P, 1, BP, 1);

    NekDouble d = Dot(E0, BP);
    NekDouble e = Dot(E1, BP);

    NekDouble det = a * c - b * b;
    NekDouble s = b * e - c * d, t = b * d - a * e;

    if (s + t <= det)
    {
        if (s < 0)
        {
            if (t < 0)
            {
                t = 0;
                s = 0;
            }
            else
            {
                s = 0;
                t = (e >= 0 ? 0 : (-e >= c ? 1 : -e / c));
            }
        }
        else if (t < 0)
        {
            t = 0;
            s = (d >= 0 ? 0 : (-d >= a ? 1 : -d / a));
        }
        else
        {
            NekDouble invDet = 1.0 / det;
            s *= invDet;
            t *= invDet;
        }
    }
    else
    {
        if (s < 0)
        {
            s = 0;
            t = 1;
        }
        else if (t < 0)
        {
            s = 1;
            t = 0;
        }
        else
        {
            NekDouble numer = c + e - b - d;
            if (numer <= 0)
            {
                s = 0;
            }
            else
            {
                NekDouble denom = a - 2 * b + c;
                s               = (numer >= denom ? 1 : numer / denom);
            }
            t = 1 - s;
        }
    }

    NodeSharedPtr point = std::shared_ptr<Node>(
        new Node(0, B[0] + s * E0[0] + t * E1[0], B[1] + s * E0[1] + t * E1[1],
                 B[2] + s * E0[2] + t * E1[2]));

    return n->Distance(point);
}

bool BLMesh::TestIntersectionEl(ElementSharedPtr e1, ElementSharedPtr e2)
{
    vector<NodeSharedPtr> ns1 = e1->GetVertexList();
    vector<NodeSharedPtr> ns2 = e2->GetVertexList();
    if (ns1[0] == ns2[0] || ns1[0] == ns2[1] || ns1[0] == ns2[2] ||
        ns1[1] == ns2[0] || ns1[1] == ns2[1] || ns1[1] == ns2[2] ||
        ns1[2] == ns2[0] || ns1[2] == ns2[1] || ns1[2] == ns2[2])
    {
        return false;
    }

    Array<OneD, NekDouble> N1(3), N2(3);
    NekDouble d1, d2;

    N1[0] = (ns1[1]->m_y - ns1[0]->m_y) * (ns1[2]->m_z - ns1[0]->m_z) -
            (ns1[2]->m_y - ns1[0]->m_y) * (ns1[1]->m_z - ns1[0]->m_z);
    N1[1] = -1.0 * ((ns1[1]->m_x - ns1[0]->m_x) * (ns1[2]->m_z - ns1[0]->m_z) -
                    (ns1[2]->m_x - ns1[0]->m_x) * (ns1[1]->m_z - ns1[0]->m_z));
    N1[2] = (ns1[1]->m_x - ns1[0]->m_x) * (ns1[2]->m_y - ns1[0]->m_y) -
            (ns1[2]->m_x - ns1[0]->m_x) * (ns1[1]->m_y - ns1[0]->m_y);

    N2[0] = (ns2[1]->m_y - ns2[0]->m_y) * (ns2[2]->m_z - ns2[0]->m_z) -
            (ns2[2]->m_y - ns2[0]->m_y) * (ns2[1]->m_z - ns2[0]->m_z);
    N2[1] = -1.0 * ((ns2[1]->m_x - ns2[0]->m_x) * (ns2[2]->m_z - ns2[0]->m_z) -
                    (ns2[2]->m_x - ns2[0]->m_x) * (ns2[1]->m_z - ns2[0]->m_z));
    N2[2] = (ns2[1]->m_x - ns2[0]->m_x) * (ns2[2]->m_y - ns2[0]->m_y) -
            (ns2[2]->m_x - ns2[0]->m_x) * (ns2[1]->m_y - ns2[0]->m_y);

    d1 = -1.0 *
         (N1[0] * ns1[0]->m_x + N1[1] * ns1[0]->m_y + N1[2] * ns1[0]->m_z);
    d2 = -1.0 *
         (N2[0] * ns2[0]->m_x + N2[1] * ns2[0]->m_y + N2[2] * ns2[0]->m_z);

    Array<OneD, NekDouble> dv1(3), dv2(3);

    dv1[0] =
        N2[0] * ns1[0]->m_x + N2[1] * ns1[0]->m_y + N2[2] * ns1[0]->m_z + d2;
    dv1[1] =
        N2[0] * ns1[1]->m_x + N2[1] * ns1[1]->m_y + N2[2] * ns1[1]->m_z + d2;
    dv1[2] =
        N2[0] * ns1[2]->m_x + N2[1] * ns1[2]->m_y + N2[2] * ns1[2]->m_z + d2;

    dv2[0] =
        N1[0] * ns2[0]->m_x + N1[1] * ns2[0]->m_y + N1[2] * ns2[0]->m_z + d1;
    dv2[1] =
        N1[0] * ns2[1]->m_x + N1[1] * ns2[1]->m_y + N1[2] * ns2[1]->m_z + d1;
    dv2[2] =
        N1[0] * ns2[2]->m_x + N1[1] * ns2[2]->m_y + N1[2] * ns2[2]->m_z + d1;

    if (sign(dv1[0], dv1[1]) && sign(dv1[1], dv1[2]))
    {
        return false;
    }
    if (sign(dv2[0], dv2[1]) && sign(dv2[1], dv2[2]))
    {
        return false;
    }

    Array<OneD, NekDouble> D(3);
    D[0] = N1[1] * N2[2] - N1[2] * N2[1];
    D[1] = -1.0 * (N1[0] * N2[2] - N1[2] * N2[0]);
    D[2] = N1[0] * N2[1] - N1[0] * N2[1];

    int base1 = 0, base2 = 0;
    if (!sign(dv2[0], dv2[1]) && sign(dv2[1], dv2[2]))
    {
        base2 = 0;
    }
    else if (!sign(dv2[1], dv2[2]) && sign(dv2[2], dv2[0]))
    {
        base2 = 1;
    }
    else if (!sign(dv2[2], dv2[0]) && sign(dv2[0], dv2[1]))
    {
        base2 = 2;
    }
    else
    {
        cout << "base not set" << endl;
    }

    if (!sign(dv1[0], dv1[1]) && sign(dv1[1], dv1[2]))
    {
        base1 = 0;
    }
    else if (!sign(dv1[1], dv1[2]) && sign(dv1[2], dv1[0]))
    {
        base1 = 1;
    }
    else if (!sign(dv1[2], dv1[0]) && sign(dv1[0], dv1[1]))
    {
        base1 = 2;
    }
    else
    {
        cout << "base not set" << endl;
    }

    Array<OneD, NekDouble> p1(3), p2(3);

    p1[0] = D[0] * ns1[0]->m_x + D[1] * ns1[0]->m_y + D[2] * ns1[0]->m_z;
    p1[1] = D[0] * ns1[1]->m_x + D[1] * ns1[1]->m_y + D[2] * ns1[1]->m_z;
    p1[2] = D[0] * ns1[2]->m_x + D[1] * ns1[2]->m_y + D[2] * ns1[2]->m_z;

    p2[0] = D[0] * ns2[0]->m_x + D[1] * ns2[0]->m_y + D[2] * ns2[0]->m_z;
    p2[1] = D[0] * ns2[1]->m_x + D[1] * ns2[1]->m_y + D[2] * ns2[1]->m_z;
    p2[2] = D[0] * ns2[2]->m_x + D[1] * ns2[2]->m_y + D[2] * ns2[2]->m_z;

    NekDouble t11, t12, t21, t22;
    int o1 = 0, o2 = 0;
    if (base1 == 0)
    {
        o1 = 1;
        o2 = 2;
    }
    else if (base1 == 1)
    {
        o1 = 2;
        o2 = 0;
    }
    else if (base1 == 2)
    {
        o1 = 0;
        o2 = 1;
    }

    t11 = p1[o1] + (p1[base1] - p1[o1]) * dv1[o1] / (dv1[o1] - dv1[base1]);
    t12 = p1[o2] + (p1[base1] - p1[o2]) * dv1[o2] / (dv1[o2] - dv1[base1]);

    if (base2 == 0)
    {
        o1 = 1;
        o2 = 2;
    }
    else if (base2 == 1)
    {
        o1 = 2;
        o2 = 0;
    }
    else if (base2 == 2)
    {
        o1 = 0;
        o2 = 1;
    }

    t21 = p2[o1] + (p2[base2] - p2[o1]) * dv2[o1] / (dv2[o1] - dv2[base2]);
    t22 = p2[o2] + (p2[base2] - p2[o2]) * dv2[o2] / (dv2[o2] - dv2[base2]);

    if (t11 > t12)
    {
        swap(t11, t12);
    }
    if (t21 > t22)
    {
        swap(t21, t22);
    }

    if (t21 < t11)
    {
        swap(t11, t21);
        swap(t12, t22);
    }

    if (!sign(t21 - t11, t22 - t11) || !sign(t21 - t12, t22 - t12))
    {
        return true;
    }

    return false;
}

void BLMesh::Shrink()
{
    map<NodeSharedPtr, blInfoSharedPtr>::iterator bit;
    bool smsh = true;
    while (smsh)
    {
        smsh = false;

        vector<ElementSharedPtr> inv;
        for (int i = 0; i < m_mesh->m_element[3].size(); i++)
        {
            ElementSharedPtr el = m_mesh->m_element[3][i];
            if (!IsPrismValid(el))
            {
                inv.push_back(el);
            }
        }

        smsh = (inv.size() > 0);

        for (int i = 0; i < inv.size(); i++)
        {
            ElementSharedPtr t = m_priToTri[inv[i]];
            vector<blInfoSharedPtr> bls;
            vector<NodeSharedPtr> ns = t->GetVertexList();
            for (int j = 0; j < ns.size(); j++)
            {
                bls.push_back(m_blData[ns[j]]);
            }
            bool repeat = true;
            while (repeat)
            {
                repeat = false;
                int mx = 0;
                for (int j = 0; j < 3; j++)
                {
                    mx = max(mx, bls[j]->bl);
                }
                ASSERTL0(mx > 0, "shrinking to nothing");
                for (int j = 0; j < 3; j++)
                {
                    if (bls[j]->bl < mx)
                    {
                        continue;
                    }
                    bls[j]->bl--;
                    bls[j]->AlignNode(m_layerT[bls[j]->bl]);
                }
                if (!IsPrismValid(inv[i]))
                {
                    repeat = true;
                }
            }
        }

        bool repeat = true;
        while (repeat)
        {
            repeat = false;
            for (bit = m_blData.begin(); bit != m_blData.end(); bit++)
            {
                vector<blInfoSharedPtr> infos = m_nToNInfo[bit->first];
                for (int i = 0; i < infos.size(); i++)
                {
                    if (bit->second->bl > infos[i]->bl + 1)
                    {
                        bit->second->bl--;
                        bit->second->AlignNode(m_layerT[bit->second->bl]);
                        repeat = true;
                    }
                }
            }
        }
    }
}

bool BLMesh::IsPrismValid(ElementSharedPtr el)
{
    NekDouble mn             = numeric_limits<double>::max();
    NekDouble mx             = -1.0 * numeric_limits<double>::max();
    vector<NodeSharedPtr> ns = el->GetVertexList();
    NekVector<NekDouble> X(6), Y(6), Z(6);
    for (int j = 0; j < ns.size(); j++)
    {
        X(j) = ns[j]->m_x;
        Y(j) = ns[j]->m_y;
        Z(j) = ns[j]->m_z;
    }
    NekVector<NekDouble> x1(6), y1(6), z1(6), x2(6), y2(6), z2(6), x3(6), y3(6),
        z3(6);

    x1 = m_deriv[0] * X;
    y1 = m_deriv[0] * Y;
    z1 = m_deriv[0] * Z;
    x2 = m_deriv[1] * X;
    y2 = m_deriv[1] * Y;
    z2 = m_deriv[1] * Z;
    x3 = m_deriv[2] * X;
    y3 = m_deriv[2] * Y;
    z3 = m_deriv[2] * Z;

    for (int j = 0; j < 6; j++)
    {
        DNekMat dxdz(3, 3, 1.0, eFULL);
        dxdz(0, 0) = x1(j);
        dxdz(0, 1) = x2(j);
        dxdz(0, 2) = x3(j);
        dxdz(1, 0) = y1(j);
        dxdz(1, 1) = y2(j);
        dxdz(1, 2) = y3(j);
        dxdz(2, 0) = z1(j);
        dxdz(2, 1) = z2(j);
        dxdz(2, 2) = z3(j);

        NekDouble jacDet =
            dxdz(0, 0) * (dxdz(1, 1) * dxdz(2, 2) - dxdz(2, 1) * dxdz(1, 2)) -
            dxdz(0, 1) * (dxdz(1, 0) * dxdz(2, 2) - dxdz(2, 0) * dxdz(1, 2)) +
            dxdz(0, 2) * (dxdz(1, 0) * dxdz(2, 1) - dxdz(2, 0) * dxdz(1, 1));
        mn = min(mn, jacDet);
        mx = max(mx, jacDet);
    }

    /*SpatialDomains::GeometrySharedPtr geom = el->GetGeom(3);
    SpatialDomains::GeomFactorsSharedPtr gfac = geom->GetGeomFactors();

    cout << mn << " " << mx << " " << (mn > 0) << " " << gfac->IsValid() <<
    endl;*/

    return mn > 0;
}

void BLMesh::BuildElements()
{
    // make prisms
    map<CADOrientation::Orientation, vector<int>> baseTri;
    map<CADOrientation::Orientation, vector<int>> topTri;

    vector<int> tmp;
    // back-base
    tmp.push_back(0);
    tmp.push_back(4);
    tmp.push_back(1);
    baseTri[CADOrientation::eBackwards] = tmp;
    tmp.clear();
    // for-base
    tmp.push_back(0);
    tmp.push_back(1);
    tmp.push_back(4);
    baseTri[CADOrientation::eForwards] = tmp;
    // back-top
    tmp.clear();
    tmp.push_back(3);
    tmp.push_back(5);
    tmp.push_back(2);
    topTri[CADOrientation::eBackwards] = tmp;
    // for-top
    tmp.clear();
    tmp.push_back(3);
    tmp.push_back(2);
    tmp.push_back(5);
    topTri[CADOrientation::eForwards] = tmp;

    ElmtConfig pconf(LibUtilities::ePrism, 1, false, false);
    ElmtConfig tconf(LibUtilities::eTriangle, 1, false, false);

    for (int i = 0; i < m_mesh->m_element[2].size(); i++)
    {
        ElementSharedPtr el = m_mesh->m_element[2][i];
        vector<unsigned int>::iterator f =
            find(m_blsurfs.begin(), m_blsurfs.end(), el->m_parentCAD->GetId());

        if (f == m_blsurfs.end())
        {
            // if this triangle is not in bl surfs continue
            continue;
        }

        vector<NodeSharedPtr> tn(3); // nodes for pseduo surface
        vector<NodeSharedPtr> pn(6); // all prism nodes
        vector<NodeSharedPtr> n       = el->GetVertexList();
        CADOrientation::Orientation o = el->m_parentCAD->Orientation();

        for (int j = 0; j < 3; j++)
        {
            pn[baseTri[o][j]] = n[j];
            pn[topTri[o][j]]  = m_blData[n[j]]->pNode;
            tn[j]             = m_blData[n[j]]->pNode;
        }

        vector<int> tags;
        tags.push_back(m_id);
        ElementSharedPtr E = GetElementFactory().CreateInstance(
            LibUtilities::ePrism, pconf, pn, tags);
        E->SetId(i);

        m_mesh->m_element[3].push_back(E);

        // tag of this element doesnt matter so can just be 1
        ElementSharedPtr T = GetElementFactory().CreateInstance(
            LibUtilities::eTriangle, tconf, tn, tags);
        m_psuedoSurface.push_back(T);

        T->m_parentCAD = el->m_parentCAD;

        m_priToTri[E] = el;
    }
}

NekDouble BLMesh::Visability(vector<ElementSharedPtr> tris,
                             Array<OneD, NekDouble> N)
{
    NekDouble mn = numeric_limits<double>::max();

    for (int i = 0; i < tris.size(); i++)
    {
        Array<OneD, NekDouble> tmp = tris[i]->Normal(true);
        NekDouble dt = tmp[0] * N[0] + tmp[1] * N[1] + tmp[2] * N[2];
        mn           = min(mn, dt);
    }
    return mn;
}

Array<OneD, NekDouble> BLMesh::GetNormal(vector<ElementSharedPtr> tris)
{
    // compile list of normals
    vector<Array<OneD, NekDouble>> N;
    for (int i = 0; i < tris.size(); i++)
    {
        N.push_back(tris[i]->Normal(true));
    }

    vector<NekDouble> w(N.size());
    Array<OneD, NekDouble> Np(3, 0.0);

    for (int i = 0; i < N.size(); i++)
    {
        w[i] = 1.0 / N.size();
    }

    for (int i = 0; i < N.size(); i++)
    {
        Np[0] += w[i] * N[i][0];
        Np[1] += w[i] * N[i][1];
        Np[2] += w[i] * N[i][2];
    }
    NekDouble mag = sqrt(Np[0] * Np[0] + Np[1] * Np[1] + Np[2] * Np[2]);
    Np[0] /= mag;
    Np[1] /= mag;
    Np[2] /= mag;

    Array<OneD, NekDouble> Ninital = Np;

    NekDouble dot = 0.0;
    int ct        = 0;
    vector<NekDouble> a(N.size());
    while (fabs(dot - 1) > 1e-6)
    {
        ct++;
        Array<OneD, NekDouble> Nplast(3);
        Nplast[0] = Np[0];
        Nplast[1] = Np[1];
        Nplast[2] = Np[2];

        NekDouble aSum = 0.0;
        for (int i = 0; i < N.size(); i++)
        {
            NekDouble dot2 =
                Np[0] * N[i][0] + Np[1] * N[i][1] + Np[2] * N[i][2];
            if (fabs(dot2 - 1) < 1e-9)
            {
                a[i] = dot2 / fabs(dot2) * 1e-9;
            }
            else
            {
                a[i] = acos(dot2);
            }

            aSum += a[i];
        }

        NekDouble wSum = 0.0;
        for (int i = 0; i < N.size(); i++)
        {
            w[i] = w[i] * a[i] / aSum;
            wSum += w[i];
        }

        for (int i = 0; i < N.size(); i++)
        {
            w[i] = w[i] / wSum;
        }

        Array<OneD, NekDouble> NpN(3, 0.0);
        for (int i = 0; i < N.size(); i++)
        {
            NpN[0] += w[i] * N[i][0];
            NpN[1] += w[i] * N[i][1];
            NpN[2] += w[i] * N[i][2];
        }
        mag = sqrt(NpN[0] * NpN[0] + NpN[1] * NpN[1] + NpN[2] * NpN[2]);
        NpN[0] /= mag;
        NpN[1] /= mag;
        NpN[2] /= mag;

        Np[0] = 0.8 * NpN[0] + (1.0 - 0.8) * Np[0];
        Np[1] = 0.8 * NpN[1] + (1.0 - 0.8) * Np[1];
        Np[2] = 0.8 * NpN[2] + (1.0 - 0.8) * Np[2];
        mag   = sqrt(Np[0] * Np[0] + Np[1] * Np[1] + Np[2] * Np[2]);
        Np[0] /= mag;
        Np[1] /= mag;
        Np[2] /= mag;

        dot = Np[0] * Nplast[0] + Np[1] * Nplast[1] + Np[2] * Nplast[2];

        if (ct > 100000)
        {
            cout << "run out of iterations" << endl;
            Np = Ninital;
            break;
        }
    }

    return Np;
}

void BLMesh::Setup()
{
    NekDouble a = 2.0 * (1.0 - m_prog) / (1.0 - pow(m_prog, m_layer + 1));
    m_layerT    = Array<OneD, NekDouble>(m_layer);
    m_layerT[0] = a * m_prog * m_bl;
    for (int i = 1; i < m_layer; i++)
    {
        m_layerT[i] = m_layerT[i - 1] + a * pow(m_prog, i) * m_bl;
    }

    if (m_mesh->m_verbose)
    {
        cout << "First layer height " << m_layerT[0] << endl;
    }

    // this sets up all the boundary layer normals data holder
    set<int> symSurfs;
    NodeSet::iterator it;
    int ct     = 0;
    int failed = 0;

    // ofstream file1;
    // file1.open("pts.3D");
    // file1 << "X Y Z value" << endl;
    for (it = m_mesh->m_vertexSet.begin(); it != m_mesh->m_vertexSet.end();
         it++, ct++)
    {
        vector<CADSurfSharedPtr> ss = (*it)->GetCADSurfs();
        vector<unsigned int> surfs;
        for (int i = 0; i < ss.size(); i++)
        {
            surfs.push_back(ss[i]->GetId());
        }
        sort(surfs.begin(), surfs.end());
        vector<unsigned int> inter, diff;

        set_intersection(m_blsurfs.begin(), m_blsurfs.end(), surfs.begin(),
                         surfs.end(), back_inserter(inter));
        set_symmetric_difference(inter.begin(), inter.end(), surfs.begin(),
                                 surfs.end(), back_inserter(diff));

        // is somewhere on a bl surface
        if (inter.size() > 0)
        {
            // initialise a new bl boudnary node
            blInfoSharedPtr bln = std::shared_ptr<blInfo>(new blInfo);
            bln->oNode          = (*it);
            bln->stopped        = false;

            // file1 << (*it)->m_x << " " << (*it)->m_y << " " << (*it)->m_z <<
            // " " << ss.size() << endl;

            if (diff.size() > 0)
            {
                // if the diff size is greater than 1 there is a curve that
                // needs remeshing
                ASSERTL0(diff.size() <= 1, "not setup for curve bl refinement");
                symSurfs.insert(diff[0]);
                bln->symsurf = diff[0];
                bln->onSym   = true;
            }
            else
            {
                bln->onSym = false;
            }

            m_blData[(*it)] = bln;
        }
    }
    // file1.close();

    // need a map from vertex idx to surface elements
    // but do not care about triangles which are not in the bl
    for (int i = 0; i < m_mesh->m_element[2].size(); i++)
    {
        vector<unsigned int>::iterator f =
            find(m_blsurfs.begin(), m_blsurfs.end(),
                 m_mesh->m_element[2][i]->m_parentCAD->GetId());

        if (f == m_blsurfs.end())
        {
            // if this triangle is not in bl surfs continue
            continue;
        }

        vector<NodeSharedPtr> ns = m_mesh->m_element[2][i]->GetVertexList();
        for (int j = 0; j < ns.size(); j++)
        {
            m_blData[ns[j]]->els.push_back(m_mesh->m_element[2][i]);
            m_blData[ns[j]]->surfs.insert(
                m_mesh->m_element[2][i]->m_parentCAD->GetId());
        }
    }

    map<NodeSharedPtr, blInfoSharedPtr>::iterator bit;
    for (bit = m_blData.begin(); bit != m_blData.end(); bit++)
    {
        // calculate mesh normal
        bit->second->N = GetNormal(bit->second->els);

        if (Visability(bit->second->els, bit->second->N) < 0.0)
        {
            cerr << "failed " << bit->first->m_x << " " << bit->first->m_y
                 << " " << bit->first->m_z << " "
                 << Visability(bit->second->els, bit->second->N) << endl;
            failed++;
        }

        Array<OneD, NekDouble> loc = bit->first->GetLoc();
        for (int k = 0; k < 3; k++)
        {
            loc[k] += bit->second->N[k] * m_layerT[0];
        }

        bit->second->pNode = std::shared_ptr<Node>(
            new Node(m_mesh->m_numNodes++, loc[0], loc[1], loc[2]));
        bit->second->bl = 0;
    }

    m_symSurfs = vector<unsigned int>(symSurfs.begin(), symSurfs.end());

    // now need to enforce that all symmetry plane nodes have their normal
    // forced onto the symmetry surface
    for (bit = m_blData.begin(); bit != m_blData.end(); bit++)
    {
        if (!bit->second->onSym)
        {
            continue;
        }

        Array<OneD, NekDouble> loc = bit->second->pNode->GetLoc();
        Array<OneD, NekDouble> uv =
            m_mesh->m_cad->GetSurf(bit->second->symsurf)->locuv(loc);

        Array<OneD, NekDouble> nl =
            m_mesh->m_cad->GetSurf(bit->second->symsurf)->P(uv);

        Array<OneD, NekDouble> N(3);
        N[0] = nl[0] - bit->first->m_x;
        N[1] = nl[1] - bit->first->m_y;
        N[2] = nl[2] - bit->first->m_z;

        NekDouble mag = sqrt(N[0] * N[0] + N[1] * N[1] + N[2] * N[2]);
        N[0] /= mag;
        N[1] /= mag;
        N[2] /= mag;

        bit->second->N = N;
        bit->second->AlignNode(m_layerT[0]);
    }

    // now smooth all the normals by distance weighted average
    // keep normals on curves constant
    for (bit = m_blData.begin(); bit != m_blData.end(); bit++)
    {
        set<int> added;
        added.insert(bit->first->m_id);
        for (int i = 0; i < bit->second->els.size(); i++)
        {
            vector<NodeSharedPtr> ns = bit->second->els[i]->GetVertexList();
            for (int j = 0; j < ns.size(); j++)
            {
                set<int>::iterator t = added.find(ns[j]->m_id);
                if (t == added.end())
                {
                    m_nToNInfo[bit->first].push_back(m_blData[ns[j]]);
                }
            }
        }
    }

    for (int l = 0; l < 10; l++)
    {
        for (bit = m_blData.begin(); bit != m_blData.end(); bit++)
        {
            if (bit->first->GetNumCADSurf() > 1)
            {
                continue;
            }

            Array<OneD, NekDouble> sumV(3, 0.0);
            vector<blInfoSharedPtr> data = m_nToNInfo[bit->first];
            NekDouble Dtotal             = 0.0;
            for (int i = 0; i < data.size(); i++)
            {
                NekDouble d = bit->first->Distance(data[i]->oNode);
                Dtotal += d;
                sumV[0] += data[i]->N[0] / d;
                sumV[1] += data[i]->N[1] / d;
                sumV[2] += data[i]->N[2] / d;
            }
            sumV[0] *= Dtotal;
            sumV[1] *= Dtotal;
            sumV[2] *= Dtotal;
            NekDouble mag =
                sqrt(sumV[0] * sumV[0] + sumV[1] * sumV[1] + sumV[2] * sumV[2]);
            sumV[0] /= mag;
            sumV[1] /= mag;
            sumV[2] /= mag;

            Array<OneD, NekDouble> N(3);

            N[0] = (1.0 - 0.8) * bit->second->N[0] + 0.8 * sumV[0];
            N[1] = (1.0 - 0.8) * bit->second->N[1] + 0.8 * sumV[1];
            N[2] = (1.0 - 0.8) * bit->second->N[2] + 0.8 * sumV[2];

            mag = sqrt(N[0] * N[0] + N[1] * N[1] + N[2] * N[2]);
            N[0] /= mag;
            N[1] /= mag;
            N[2] /= mag;

            bit->second->N = N;
            bit->second->AlignNode(m_layerT[0]);
        }
    }

    /*ofstream file;
    file.open("bl.lines");
    for(bit = m_blData.begin(); bit != m_blData.end(); bit++)
    {
        NekDouble l = 0.05;
        file << bit->first->m_x << ", " << bit->first->m_y << ", " <<
    bit->first->m_z << endl;
        file << bit->first->m_x + bit->second->N[0]*l << ", "
             << bit->first->m_y + bit->second->N[1]*l << ", "
             << bit->first->m_z + bit->second->N[2]*l << endl;
        file << endl;
    }
    file.close();*/

    ASSERTL0(failed == 0, "some normals failed to generate");

    LibUtilities::PointsKey pkey1(2, LibUtilities::eNodalPrismElec);

    Array<OneD, NekDouble> u1, v1, w1;
    LibUtilities::PointsManager()[pkey1]->GetPoints(u1, v1, w1);

    LibUtilities::NodalUtilPrism nodalPrism(1, u1, v1, w1);

    NekMatrix<NekDouble> Vandermonde  = *nodalPrism.GetVandermonde();
    NekMatrix<NekDouble> VandermondeI = Vandermonde;
    VandermondeI.Invert();

    m_deriv[0] = *nodalPrism.GetVandermondeForDeriv(0) * VandermondeI;
    m_deriv[1] = *nodalPrism.GetVandermondeForDeriv(1) * VandermondeI;
    m_deriv[2] = *nodalPrism.GetVandermondeForDeriv(2) * VandermondeI;
}
} // namespace NekMeshUtils
} // namespace Nektar

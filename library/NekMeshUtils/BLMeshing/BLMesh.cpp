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
//  Description: tet meshing methods
//
////////////////////////////////////////////////////////////////////////////////

#include <NekMeshUtils/BLMeshing/BLMesh.h>
#include <NekMeshUtils/CADSystem/CADSurf.h>

using namespace std;
namespace Nektar
{
namespace NekMeshUtils
{

void BLMesh::Mesh()
{
    // At this stage the surface mesh is complete and the elements know their
    // neighbours through element links in the edges, this includes quads.

    // here elements are made for the boundary layer they will need to know
    // links (maybe facelinks), so that the tetmeshing module can extract the
    // surface upon which it needs to mesh (top of the bl and the rest of the
    // surface).

    vector<ElementSharedPtr> quad;
    vector<ElementSharedPtr> ptri; // triangles to grow prisms onto

    for (int i = 0; i < m_mesh->m_element[2].size(); i++)
    {
        bool onblsurf = false;
        for (int j = 0; j < m_blsurfs.size(); j++)
        {
            if (m_mesh->m_element[2][i]->CADSurfId == m_blsurfs[j])
            {
                onblsurf = true;
                break;
            }
        }
        if (m_mesh->m_element[2][i]->GetConf().m_e ==
            LibUtilities::eQuadrilateral)
        {
            quad.push_back(m_mesh->m_element[2][i]);
        }
        else if (onblsurf)
        {
            ptri.push_back(m_mesh->m_element[2][i]);
        }
    }

    map<int, NodeSharedPtr> blpair;
    for (int i = 0; i < quad.size(); i++)
    {
        vector<EdgeSharedPtr> e = quad[i]->GetEdgeList();
        for (int j = 0; j < e.size(); j++)
        {
            // if both or none are on curve skip
            if ((e[j]->m_n1->GetNumCadCurve() > 0 &&
                 e[j]->m_n2->GetNumCadCurve() > 0) ||
                (!(e[j]->m_n1->GetNumCadCurve() > 0) &&
                 !(e[j]->m_n2->GetNumCadCurve() > 0)))
            {
                continue;
            }

            if (e[j]->m_n1->GetNumCadCurve() > 0)
            {
                blpair[e[j]->m_n1->m_id] = e[j]->m_n2;
            }
            else if (e[j]->m_n2->GetNumCadCurve() > 0)
            {
                blpair[e[j]->m_n2->m_id] = e[j]->m_n1;
            }
            else
            {
                ASSERTL0(false, "that failed");
            }
        }
    }

    map<int, int> nm;
    nm[0] = 0;
    nm[1] = 3;
    nm[2] = 4;
    nm[3] = 5;
    nm[4] = 1;
    nm[5] = 2;

    for (int i = 0; i < ptri.size(); i++)
    {
        vector<NodeSharedPtr> pn(6); // all prism nodes
        vector<NodeSharedPtr> n = ptri[i]->GetVertexList();

        vector<pair<int, CADSurfSharedPtr> > tmpss = n[0]->GetCADSurfs();
        CADSurfSharedPtr tmps;

        for (int j = 0; j < tmpss.size(); j++)
        {
            if (tmpss[j].first == ptri[i]->CADSurfId)
            {
                tmps = tmpss[j].second;
                break;
            }
        }

        if (tmps->IsReversedNormal())
        {
            nm[0] = 0;
            nm[1] = 3;
            nm[2] = 1;
            nm[3] = 2;
            nm[4] = 4;
            nm[5] = 5;
        }

        for (int j = 0; j < n.size(); j++)
        {
            pn[nm[j * 2]] = n[j];

            map<int, NodeSharedPtr>::iterator it;
            it = blpair.find(n[j]->m_id);
            if (it != blpair.end())
            {
                pn[nm[j * 2 + 1]] = blpair[n[j]->m_id];
            }
            else
            {
                Array<OneD, NekDouble> AN(3);

                vector<pair<int, CADSurfSharedPtr> > surfs =
                    n[j]->GetCADSurfs();

                for (int s = 0; s < surfs.size(); s++)
                {
                    Array<OneD, NekDouble> N = surfs[s].second->N(
                        n[j]->GetCADSurfInfo(surfs[s].first));
                    for (int k = 0; k < 3; k++)
                    {
                        AN[k] += N[k];
                    }
                }

                NekDouble mag =
                    sqrt(AN[0] * AN[0] + AN[1] * AN[1] + AN[2] * AN[2]);

                for (int k = 0; k < 3; k++)
                {
                    AN[k] /= mag;
                }

                Array<OneD, NekDouble> loc = n[j]->GetLoc();
                Array<OneD, NekDouble> np(3);
                for (int k = 0; k < 3; k++)
                    np[k]        = loc[k] + AN[k] * m_bl;
                NodeSharedPtr nn = boost::shared_ptr<Node>(
                    new Node(m_mesh->m_numNodes++, np[0], np[1], np[2]));
                pn[nm[j * 2 + 1]]  = nn;
                blpair[n[j]->m_id] = nn;
            }
        }

        ElmtConfig conf(LibUtilities::ePrism, 1, false, false);
        vector<int> tags;
        tags.push_back(1); // all prisms are comp 1
        ElementSharedPtr E = GetElementFactory().CreateInstance(
            LibUtilities::ePrism, conf, pn, tags);

        m_mesh->m_element[3].push_back(E);

        // need to give the surface element some information about
        // which prism is above it
        // so that tetmesh can infer the pseudo surface
        vector<NodeSharedPtr> faceNodes;
        vector<EdgeSharedPtr> edgeList = ptri[i]->GetEdgeList();
        FaceSharedPtr F = boost::shared_ptr<Face>(new Face(
            n, faceNodes, edgeList, ptri[i]->GetConf().m_faceCurveType));
        vector<FaceSharedPtr> f = E->GetFaceList();
        for (int j = 0; j < f.size(); j++)
        {
            if (f[j]->m_vertexList.size() != 3) // quad
                continue;

            // only two triangle faces so if its not this one, this is the
            // pseudo surfaces
            if (!(F == f[j]))
            {
                m_surftopriface[ptri[i]->GetId()] = f[j];
            }
        }
    }
}
}
}

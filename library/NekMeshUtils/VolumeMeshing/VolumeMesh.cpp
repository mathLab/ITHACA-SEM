////////////////////////////////////////////////////////////////////////////////
//
//  File: VolumeMesh.cpp
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
//  Description: Process volume meshing.
//
////////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/ParseUtils.h>

#include "VolumeMesh.h"
#include <NekMeshUtils/CADSystem/CADCurve.h>
#include <NekMeshUtils/CADSystem/CADSurf.h>
#include <NekMeshUtils/SurfaceMeshing/CurveMesh.h>
#include <NekMeshUtils/SurfaceMeshing/FaceMesh.h>
#include <NekMeshUtils/VolumeMeshing/BLMeshing/BLMesh.h>
#include <NekMeshUtils/VolumeMeshing/TetMeshing/TetMesh.h>

using namespace std;
namespace Nektar
{
namespace NekMeshUtils
{

ModuleKey VolumeMesh::className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eProcessModule, "volumemesh"), VolumeMesh::create,
    "Generates a volume mesh");

VolumeMesh::VolumeMesh(MeshSharedPtr m) : ProcessModule(m)
{
    m_config["blsurfs"] =
        ConfigOption(false, "0", "Generate prisms on these surfs");
    m_config["blthick"]  = ConfigOption(false, "0", "Prism layer thickness");
    m_config["bllayers"] = ConfigOption(false, "0", "Prism layers");
    m_config["blprog"]   = ConfigOption(false, "0", "Prism progression");
}

VolumeMesh::~VolumeMesh()
{
}

void VolumeMesh::Process()
{
    if (m_mesh->m_verbose)
        cout << endl << "Volume meshing" << endl;

    bool makeBL;
    vector<unsigned int> blSurfs;

    if (m_config["blsurfs"].beenSet)
    {
        makeBL = true;
        ParseUtils::GenerateSeqVector(m_config["blsurfs"].as<string>(),
                                      blSurfs);
    }
    else
    {
        makeBL = false;
    }

    NekDouble prefix = 100;
    if (m_mesh->m_cad->GetNumSurf() > 100)
    {
        prefix *= 10;
    }

    TetMeshSharedPtr tet;
    if (makeBL)
    {
        BLMeshSharedPtr blmesh = MemoryManager<BLMesh>::AllocateSharedPtr(
            m_mesh, blSurfs, m_config["blthick"].as<NekDouble>(),
            m_config["bllayers"].as<int>(), m_config["blprog"].as<NekDouble>(),
            prefix + 1);

        blmesh->Mesh();

        // remesh the correct surfaces
        vector<unsigned int> symsurfs = blmesh->GetSymSurfs();
        vector<ElementSharedPtr> els  = m_mesh->m_element[2];
        m_mesh->m_element[2].clear();
        for (int i = 0; i < els.size(); i++)
        {
            vector<unsigned int>::iterator f = find(
                symsurfs.begin(), symsurfs.end(), els[i]->m_parentCAD->GetId());

            if (f == symsurfs.end())
            {
                m_mesh->m_element[2].push_back(els[i]);
            }
            else
            {
                // remove element from links
                vector<EdgeSharedPtr> es = els[i]->GetEdgeList();
                for (int j = 0; j < es.size(); j++)
                {
		    vector<pair<weak_ptr<Element>, int> > lk = es[j]->m_elLink;
                    es[j]->m_elLink.clear();
                    for (int k = 0; k < lk.size(); k++)
                    {
                        if (lk[k].first.lock() == els[i])
                        {
                            continue;
                        }
                        es[j]->m_elLink.push_back(lk[k]);
                    }
                }
            }
        }

        for (int i = 0; i < symsurfs.size(); i++)
        {
            set<int> cIds;
            vector<EdgeLoopSharedPtr> e =
                m_mesh->m_cad->GetSurf(symsurfs[i])->GetEdges();
            for (int k = 0; k < e.size(); k++)
            {
                for (int j = 0; j < e[k]->edges.size(); j++)
                {
                    cIds.insert(e[k]->edges[j]->GetId());
                }
            }

            // find the curve nodes which are on this symsurf
            map<int, vector<NodeSharedPtr> > curveNodeMap;
            NodeSet::iterator it;
            for (it = m_mesh->m_vertexSet.begin();
                 it != m_mesh->m_vertexSet.end(); it++)
            {
                vector<CADCurveSharedPtr> cc = (*it)->GetCADCurves();
                for (int j = 0; j < cc.size(); j++)
                {
                    set<int>::iterator f = cIds.find(cc[j]->GetId());
                    if (f != cIds.end())
                    {
                        curveNodeMap[cc[j]->GetId()].push_back((*it));
                    }
                }
            }

            // need to bubble sort the vectors
            map<int, vector<NodeSharedPtr> >::iterator cit;
            for (cit = curveNodeMap.begin(); cit != curveNodeMap.end(); cit++)
            {
                vector<NekDouble> ts;
                for (int i = 0; i < cit->second.size(); i++)
                {
                    ts.push_back(cit->second[i]->GetCADCurveInfo(cit->first));
                }
                bool repeat = true;
                while (repeat)
                {
                    repeat = false;
                    for (int i = 0; i < ts.size() - 1; i++)
                    {
                        if (ts[i] > ts[i + 1])
                        {
                            swap(ts[i], ts[i + 1]);
                            swap(cit->second[i], cit->second[i + 1]);
                            repeat = true;
                            break;
                        }
                    }
                }
            }

            // create quads
            map<NodeSharedPtr, NodeSharedPtr> nmap = blmesh->GetSymNodes();
            for (cit = curveNodeMap.begin(); cit != curveNodeMap.end(); cit++)
            {
                for (int j = 0; j < cit->second.size() - 1; j++)
                {
                    map<NodeSharedPtr, NodeSharedPtr>::iterator f1 =
                        nmap.find(cit->second[j]);
                    map<NodeSharedPtr, NodeSharedPtr>::iterator f2 =
                        nmap.find(cit->second[j + 1]);

                    if (f1 == nmap.end() || f2 == nmap.end())
                    {
                        continue;
                    }

                    NodeSharedPtr n1 = f1->second;
                    NodeSharedPtr n2 = f2->second;

                    vector<NodeSharedPtr> ns;
                    ns.push_back(cit->second[j]);
                    ns.push_back(n1);
                    ns.push_back(n2);
                    ns.push_back(cit->second[j + 1]);

                    ElmtConfig conf(LibUtilities::eQuadrilateral, 1, false,
                                    false);

                    vector<int> tags;
                    tags.push_back(prefix * 2 + symsurfs[i]);
                    ElementSharedPtr E = GetElementFactory().CreateInstance(
                        LibUtilities::eQuadrilateral, conf, ns, tags);
                    E->m_parentCAD = m_mesh->m_cad->GetSurf(symsurfs[i]);
                    m_mesh->m_element[2].push_back(E);

                    // need to dummy process the new elements
                    for (int k = 0; k < E->GetEdgeCount(); ++k)
                    {
                        pair<EdgeSet::iterator, bool> testIns;
                        EdgeSharedPtr ed = E->GetEdge(k);
                        testIns          = m_mesh->m_edgeSet.insert(ed);

                        if (testIns.second)
                        {
                            EdgeSharedPtr ed2 = *testIns.first;
                            ed2->m_elLink.push_back(
                                pair<ElementSharedPtr, int>(E, k));
                        }
                        else
                        {
                            EdgeSharedPtr e2 = *(testIns.first);
                            E->SetEdge(k, e2);
                            e2->m_elLink.push_back(
                                pair<ElementSharedPtr, int>(E, k));
                        }
                    }
                }
            }

            // swap nodes
            for (cit = curveNodeMap.begin(); cit != curveNodeMap.end(); cit++)
            {
                for (int j = 0; j < cit->second.size(); j++)
                {
                    map<NodeSharedPtr, NodeSharedPtr>::iterator f1 =
                        nmap.find(cit->second[j]);
                    if (f1 == nmap.end())
                    {
                        continue;
                    }
                    cit->second[j] = f1->second;
                }
            }
            map<int, CurveMeshSharedPtr> cm;
            for (cit = curveNodeMap.begin(); cit != curveNodeMap.end(); cit++)
            {
                cm[cit->first] = MemoryManager<CurveMesh>::AllocateSharedPtr(
                    cit->first, m_mesh, cit->second);
            }

            FaceMeshSharedPtr f = MemoryManager<FaceMesh>::AllocateSharedPtr(
                symsurfs[i], m_mesh, cm, symsurfs[i]);
            f->Mesh();
        }

        vector<unsigned int> blsurfs = blmesh->GetBLSurfs();

        // build the surface for tetgen to use.
        vector<ElementSharedPtr> tetsurface = blmesh->GetPseudoSurface();
        for (int i = 0; i < m_mesh->m_element[2].size(); i++)
        {
            if (m_mesh->m_element[2][i]->GetConf().m_e ==
                LibUtilities::eQuadrilateral)
            {
                continue;
            }

            vector<unsigned int>::iterator f =
                find(blsurfs.begin(), blsurfs.end(),
                     m_mesh->m_element[2][i]->m_parentCAD->GetId());

            if (f == blsurfs.end())
            {
                tetsurface.push_back(m_mesh->m_element[2][i]);
            }
        }

        tet = MemoryManager<TetMesh>::AllocateSharedPtr(m_mesh, prefix, tetsurface);
    }
    else
    {
        tet = MemoryManager<TetMesh>::AllocateSharedPtr(m_mesh, prefix);
    }

    tet->Mesh();

    ClearElementLinks();
    ProcessVertices();
    ProcessEdges();
    ProcessFaces();
    ProcessElements();
    ProcessComposites();

    if (m_mesh->m_verbose)
    {
        cout << endl;
    }
}
}
}

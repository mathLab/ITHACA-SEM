////////////////////////////////////////////////////////////////////////////////
//
//  File: SurfaceMeshing.cpp
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
//  Description: surfacemeshing object methods.
//
////////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/ParseUtils.hpp>

#include "VolumeMesh.h"
#include <NekMeshUtils/VolumeMeshing/BLMeshing/BLMesh.h>
#include <NekMeshUtils/VolumeMeshing/TetMeshing/TetMesh.h>

using namespace std;
namespace Nektar
{
namespace NekMeshUtils
{

ModuleKey VolumeMesh::className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eProcessModule, "volumemesh"),
    VolumeMesh::create,
    "Generates a volume mesh");

VolumeMesh::VolumeMesh(MeshSharedPtr m) : ProcessModule(m)
{
    m_config["blsurfs"] =
        ConfigOption(false, "0", "Generate prisms on these surfs");
    m_config["blthick"] =
        ConfigOption(false, "0", "Prism layer thickness");
    m_config["bllayers"] =
        ConfigOption(false, "0", "Prism layers");
    m_config["blprog"] =
        ConfigOption(false, "0", "Prism progression");
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

    if(m_config["blsurfs"].beenSet)
    {
        makeBL = true;
        m_mesh->m_numcomp = 2;
        ParseUtils::GenerateSeqVector(m_config["blsurfs"].as<string>().c_str(),
                                      blSurfs);
    }
    else
    {
        makeBL = false;
        m_mesh->m_numcomp = 1;
    }

    TetMeshSharedPtr tet;
    if (makeBL)
    {
        BLMeshSharedPtr blmesh = MemoryManager<BLMesh>::AllocateSharedPtr(
                                        m_mesh, blSurfs,
                                        m_config["blthick"].as<NekDouble>(),
                                        m_config["bllayers"].as<int>(),
                                        m_config["blprog"].as<NekDouble>());

        blmesh->Mesh();

        m_mesh->m_element[2].clear();
        /*m_mesh->m_vertexSet.clear();
        m_mesh->m_edgeSet.clear();
        int eid=0, fid=0, elid=0;
        for(int i = 0; i < m_mesh->m_element[3].size(); i++)
        {
            vector<NodeSharedPtr> ns = m_mesh->m_element[3][i]->GetVertexList();
            for(int j = 0; j < ns.size(); j++)
            {
                m_mesh->m_vertexSet.insert(ns[j]);
            }
            vector<EdgeSharedPtr> es = m_mesh->m_element[3][i]->GetEdgeList();
            for(int j = 0; j < es.size(); j++)
            {
                es[j]->m_id = eid++;
                m_mesh->m_edgeSet.insert(es[j]);
            }
            vector<FaceSharedPtr> fs = m_mesh->m_element[3][i]->GetFaceList();
            for(int j = 0; j < fs.size(); j++)
            {
                fs[j]->m_id = fid++;
                m_mesh->m_faceSet.insert(fs[j]);
            }
            m_mesh->m_element[3][i]->SetId(elid++);
        }
        set<NodeSharedPtr> tmp(m_mesh->m_vertexSet.begin(),
                               m_mesh->m_vertexSet.end());
        set<NodeSharedPtr>::iterator it;
        for(it = tmp.begin(); it != tmp.end(); it++)
        {
            cout << (*it)->m_id << endl;
        }*/
        ClearElementLinks();
        ProcessVertices();
        ProcessEdges();
        ProcessFaces();
        ProcessElements();
        ProcessComposites();
        return;
    }
    else
    {
        tet = MemoryManager<TetMesh>::AllocateSharedPtr(m_mesh);
        tet->Mesh();
    }

    ClearElementLinks();
    ProcessVertices();
    ProcessEdges();
    ProcessFaces();
    ProcessElements();
    ProcessComposites();

    if (m_mesh->m_verbose)
        cout << endl;
}
}
}

////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessInsertSurface.cpp
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
//  Description: Calculate Jacobians of elements.
//
////////////////////////////////////////////////////////////////////////////////

#include <NekMeshUtils/MeshElements/Element.h>
#include "ProcessInsertSurface.h"

#include <ANN/ANN.h>

using namespace std;
using namespace Nektar::NekMeshUtils;

namespace Nektar
{
namespace Utilities
{

ModuleKey ProcessInsertSurface::className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eProcessModule, "insertsurface"),
    ProcessInsertSurface::create,
    "Insert high-order surface mesh into current working mesh.");

ProcessInsertSurface::ProcessInsertSurface(MeshSharedPtr m) : ProcessModule(m)
{
    m_config["mesh"] =
        ConfigOption(false, "", "Mesh to be inserted.");
    m_config["nonconforming"] =
        ConfigOption(false,"", "Relax tests for nonconforming boundries");
}

ProcessInsertSurface::~ProcessInsertSurface()
{
}

void ProcessInsertSurface::Process()
{
    if (m_mesh->m_verbose)
    {
        cout << "ProcessInsertSurface: Inserting mesh... " << endl;
    }

    string file = m_config["mesh"].as<string>();
    bool nonconform = m_config["nonconforming"].beenSet;

    if (m_mesh->m_verbose)
    {
        cout << "inserting surface from " << file << endl;
    }
    MeshSharedPtr inMsh = boost::shared_ptr<Mesh>(new Mesh());
    inMsh->m_verbose = m_mesh->m_verbose;
    ModuleSharedPtr mod = GetModuleFactory().CreateInstance(
        ModuleKey(eInputModule, "xml"), inMsh);
    mod->RegisterConfig("infile", file);
    mod->Process();

    //build ann tree of surface verticies from inMsh
    //match surface vertices in ccm mesh to inMsh and copy information

    //tolerance of matching vertices
    NekDouble tol = 1e-5;

    NodeSet surfaceNodes;
    for(int i = 0; i < inMsh->m_element[2].size(); i++)
    {
        vector<NodeSharedPtr> ns = inMsh->m_element[2][i]->GetVertexList();
        for(int j = 0; j < ns.size(); j++)
        {
            surfaceNodes.insert(ns[j]);
        }
    }

    vector<NodeSharedPtr> inMshnodeList(surfaceNodes.begin(), surfaceNodes.end());

    ANNpointArray dataPts;
    ANNpoint queryPt;
    ANNidxArray nnIdx;
    ANNdistArray dists;
    ANNkd_tree* kdTree;

    queryPt = annAllocPt(3);
    dataPts = annAllocPts(inMshnodeList.size(),3);

    for(int i = 0; i < inMshnodeList.size(); i++)
    {
        dataPts[i][0] = inMshnodeList[i]->m_x;
        dataPts[i][1] = inMshnodeList[i]->m_y;
        dataPts[i][2] = inMshnodeList[i]->m_z;
    }

    kdTree = new ANNkd_tree(dataPts, inMshnodeList.size(), 3);

    int sample = 1;
    nnIdx = new ANNidx[sample];
    dists = new ANNdist[sample];

    surfaceNodes.clear();
    for(int i = 0; i < m_mesh->m_element[2].size(); i++)
    {
        vector<NodeSharedPtr> ns = m_mesh->m_element[2][i]->GetVertexList();
        for(int j = 0; j < ns.size(); j++)
        {
            surfaceNodes.insert(ns[j]);
        }
    }

    if(!nonconform)
    {
        ASSERTL0(surfaceNodes.size() == inMshnodeList.size(),
                 "surface mesh node count mismatch, will not work");
    }

    EdgeSet surfEdges;
    for(int i = 0; i < m_mesh->m_element[2].size(); i++)
    {
        FaceSharedPtr f = m_mesh->m_element[2][i]->GetFaceLink();
        vector<EdgeSharedPtr> es = f->m_edgeList;
        for(int j = 0; j < es.size(); j++)
        {
            surfEdges.insert(es[j]);
        }
    }

    EdgeSet::iterator it;
    for(it = surfEdges.begin(); it != surfEdges.end(); it++)
    {
        queryPt[0] = (*it)->m_n1->m_x;
        queryPt[1] = (*it)->m_n1->m_y;
        queryPt[2] = (*it)->m_n1->m_z;
        kdTree->annkSearch(queryPt,sample,nnIdx,dists);
        if(nonconform)
        {
            if(sqrt(dists[0]) > tol)
                continue;
        }
        else
        {
            ASSERTL0(sqrt(dists[0]) < tol, "cannot locate point accurately enough");
        }

        NodeSharedPtr inN1 = inMshnodeList[nnIdx[0]];

        queryPt[0] = (*it)->m_n2->m_x;
        queryPt[1] = (*it)->m_n2->m_y;
        queryPt[2] = (*it)->m_n2->m_z;
        kdTree->annkSearch(queryPt,sample,nnIdx,dists);
        if(nonconform)
        {
            if(sqrt(dists[0]) > tol)
                continue;
        }
        else
        {
            ASSERTL0(sqrt(dists[0]) < tol, "cannot locate point accurately enough");
        }
        NodeSharedPtr inN2 = inMshnodeList[nnIdx[0]];

        EdgeSharedPtr tst = boost::shared_ptr<Edge>(new Edge(inN1,inN2));

        EdgeSet::iterator f = inMsh->m_edgeSet.find(tst);

        ASSERTL0(f != inMsh->m_edgeSet.end(),"could not find edge in input");

        (*it)->m_edgeNodes = (*f)->m_edgeNodes;
        (*it)->m_curveType = (*f)->m_curveType;

        if((*f)->m_n1->Distance((*it)->m_n1) > tol)
        {
            reverse((*it)->m_edgeNodes.begin(),(*it)->m_edgeNodes.end());
        }
    }
}
}
}

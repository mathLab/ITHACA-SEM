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

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

#include <NekMeshUtils/MeshElements/Element.h>
#include "ProcessInsertSurface.h"

namespace bg  = boost::geometry;
namespace bgi = boost::geometry::index;

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
    typedef bg::model::point<NekDouble, 3, bg::cs::cartesian> Point;
    typedef pair<Point, unsigned int> PointI;

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
    MeshSharedPtr inMsh = std::shared_ptr<Mesh>(new Mesh());
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

    vector<PointI> dataPts;
    for(int i = 0; i < inMshnodeList.size(); i++)
    {
         dataPts.push_back(make_pair(Point( inMshnodeList[i]->m_x,
                                            inMshnodeList[i]->m_y,
                                            inMshnodeList[i]->m_z), i));
    }

    //Build tree
    bgi::rtree<PointI, bgi::rstar<16> > rtree;
    rtree.insert(dataPts.begin(), dataPts.end());

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

    for(auto &it : surfEdges)
    {
        Point queryPt1(it->m_n1->m_x, it->m_n1->m_y, it->m_n1->m_z);
        vector<PointI> result;
        rtree.query(bgi::nearest(queryPt1, 1), std::back_inserter(result));

        NekDouble dist1 = bg::distance(result[0].first, queryPt1);
        if(nonconform)
        {
            if(dist1 > tol)
            {
                continue;
            }
        }
        else
        {
            ASSERTL0(dist1 < tol, "cannot locate point accurately enough");
        }

        NodeSharedPtr inN1 = inMshnodeList[result[0].second];

        Point queryPt2(it->m_n2->m_x, it->m_n2->m_y, it->m_n2->m_z);
        result.clear();
        rtree.query(bgi::nearest(queryPt2, 1), std::back_inserter(result));

        NekDouble dist2 = bg::distance(result[0].first, queryPt2);
        if(nonconform)
        {
            if(dist2 > tol)
            {
                continue;
            }
        }
        else
        {
            ASSERTL0(dist2 < tol, "cannot locate point accurately enough");
        }
        NodeSharedPtr inN2 = inMshnodeList[result[0].second];

        EdgeSharedPtr tst = std::shared_ptr<Edge>(new Edge(inN1,inN2));

        auto f = inMsh->m_edgeSet.find(tst);

        ASSERTL0(f != inMsh->m_edgeSet.end(),"could not find edge in input");

        it->m_edgeNodes = (*f)->m_edgeNodes;
        it->m_curveType = (*f)->m_curveType;

        if((*f)->m_n1->Distance(it->m_n1) > tol)
        {
            reverse(it->m_edgeNodes.begin(), it->m_edgeNodes.end());
        }
    }
}
}
}

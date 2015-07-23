////////////////////////////////////////////////////////////////////////////////
//
//  File: InputCAD.cpp
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
//  Description: GMSH converter.
//
////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <iostream>
using namespace std;

#include <LibUtilities/CADSystem/CADSystem.h>
#include <MeshUtils/Octree.h>
#include <MeshUtils/SurfaceMeshing.h>
#include <MeshUtils/TetMesh.h>
#include <MeshUtils/MeshElem.hpp>

#include <LibUtilities/BasicUtils/SharedArray.hpp>

#include "../MeshElements.h"
#include "InputCAD.h"

namespace Nektar
{
namespace Utilities
{


    ModuleKey InputCAD::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eInputModule, "CAD"), InputCAD::create,
        "Reads CAD geometry and will generate the mesh file.");

    /**
     * @brief Set up InputCAD object.
     */
    InputCAD::InputCAD(MeshSharedPtr m) : InputModule(m)
    {
        m_config["min"] = ConfigOption(false,"-1",
                "minimum delta to be in mesh");
        m_config["max"] = ConfigOption(false,"-1",
                "maximum delta to be in mesh");
        m_config["eps"] = ConfigOption(false,"-1",
                "sensitivity to curvature");
        m_config["order"] = ConfigOption(false,"-1",
                "order of the mesh to be produced");
    }

    InputCAD::~InputCAD()
    {

    }


    void InputCAD::Process()
    {
        string CADName = m_config["infile"].as<string>();
        NekDouble m_minDelta = m_config["min"].as<NekDouble>();
        NekDouble m_maxDelta = m_config["max"].as<NekDouble>();
        NekDouble m_eps = m_config["eps"].as<NekDouble>();
        int m_order = m_config["order"].as<int>();

        ASSERTL0(!(m_minDelta == -1 || m_maxDelta == -1 ||
                   m_eps == -1 || m_order == -1),
                 "User parameters required");


        LibUtilities::CADSystemSharedPtr m_cad =
            MemoryManager<LibUtilities::CADSystem>::AllocateSharedPtr(CADName);

        if(m_mesh->m_verbose)
            cout << "Building mesh for: " << m_cad->GetName() << endl;

        ASSERTL0(m_cad->LoadCAD(),
                 "Failed to load CAD");

        if(m_mesh->m_verbose)
        {
            cout << "With parameters:" << endl;
            cout << "\tmin delta: " << m_minDelta << endl
                 << "\tmax delta: " << m_maxDelta << endl
                 << "\tesp: " << m_eps << endl
                 << "\torder: " << m_order << endl;
            m_cad->Report();
        }


        MeshUtils::OctreeSharedPtr m_octree =
            MemoryManager<MeshUtils::Octree>::AllocateSharedPtr(m_cad,
                                    m_mesh->m_verbose);

        m_octree->Build(m_minDelta, m_maxDelta, m_eps);

        MeshUtils::SurfaceMeshingSharedPtr m_surfacemeshing =
            MemoryManager<MeshUtils::SurfaceMeshing>::
                AllocateSharedPtr(m_mesh->m_verbose,m_cad,m_octree,m_order);

        m_surfacemeshing->Mesh();

        m_surfacemeshing->HOSurf();

        map<int, MeshUtils::MeshTriSharedPtr> Tris;
        map<int, MeshUtils::MeshEdgeSharedPtr> Edges;
        map<int, MeshUtils::MeshNodeSharedPtr> Nodes;
        m_surfacemeshing->Get(Nodes, Edges, Tris);

        MeshUtils::TetMeshSharedPtr m_tet =
            MemoryManager<MeshUtils::TetMesh>::AllocateSharedPtr(
                0, m_mesh->m_verbose, m_octree);

        m_mesh->m_expDim = 2;
        m_mesh->m_spaceDim = 3;
        m_mesh->m_order = m_order;

        m_mesh->m_fields.push_back("u");
        m_mesh->m_fields.push_back("v");
        m_mesh->m_fields.push_back("p");

        map<int, MeshUtils::MeshNodeSharedPtr>::iterator nit;
        map<int, MeshUtils::MeshTriSharedPtr>::iterator trit;
        map<int, NodeSharedPtr> allnodes;

        for(nit = Nodes.begin(); nit != Nodes.end(); nit++)
        {
            Array<OneD, NekDouble> loc = nit->second->GetLoc();
            NodeSharedPtr nn =
                    boost::shared_ptr<Node>(
                                new Node(nit->first,loc[0],
                                         loc[1],loc[2]));
            allnodes[nit->first] = nn;
        }


        for(trit = Tris.begin(); trit != Tris.end(); trit++)
        {
            bool add = false;

            Array<OneD, int> n = trit->second->GetN();



            vector<NodeSharedPtr> localnode;
            for(int j = 0; j < 3; j++)
            {
                if(Nodes[n[j]]->IsOnSurf(8) || Nodes[n[j]]->IsOnSurf(7) || Nodes[n[j]]->IsOnSurf(9))
                {
                    add = true;
                }
                localnode.push_back(allnodes[n[j]]);
            }
            if(!add)
                continue;

            Array<OneD, int> eg = trit->second->GetE();
            for(int j = 0; j < 3; j++)
            {
                MeshUtils::MeshEdgeSharedPtr e;
                e = Edges[eg[j]];
                vector<int> hon = e->GetHONodes(n[j]);
                for(int k = 0; k < hon.size(); k++)
                {
                    localnode.push_back(allnodes[hon[k]]);
                }
            }

            ElmtConfig conf(LibUtilities::eTriangle,m_order,false,false,false);
            vector<int> tags;
            tags.push_back(trit->second->Getcid());
            ElementSharedPtr E = GetElementFactory().
                        CreateInstance(LibUtilities::eTriangle,
                                       conf,localnode,tags);
            m_mesh->m_element[2].push_back(E);

        }

        /*for(int i = 1; i <=m_cad->GetNumSurf(); i++)
        {
            //m_surfacemeshing->Extract(i,numtris,numppt,points);

            for(int j = 0; j < numtris; j++)
            {
                vector<NodeSharedPtr> nodes;
                for(int k = 0; k < numppt; k++) //numppt
                {
                    NodeSharedPtr n =
                    boost::shared_ptr<Node>(
                                new Node(nodeCounter++,points[j][k][0],
                                         points[j][k][1],points[j][k][2]));
                    nodes.push_back(n);
                }


                ElmtConfig conf(LibUtilities::eTriangle,m_order,true,false,false);

                vector<int> tags;
                tags.push_back(i);
                ElementSharedPtr E = GetElementFactory().
                            CreateInstance(LibUtilities::eTriangle,
                                           conf,nodes,tags);
                m_mesh->m_element[2].push_back(E);

            }
        }*/

        ProcessVertices  ();
        ProcessEdges     ();
        ProcessFaces     ();
        ProcessElements  ();
        ProcessComposites();

    }


}
}

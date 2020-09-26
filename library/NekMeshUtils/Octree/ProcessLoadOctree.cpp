////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessLoadOctree.cpp
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
//  Description: Load octree.
//
////////////////////////////////////////////////////////////////////////////////

#include "ProcessLoadOctree.h"
#include <NekMeshUtils/Octree/Octree.h>

using namespace std;
namespace Nektar
{
namespace NekMeshUtils
{

ModuleKey ProcessLoadOctree::className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eProcessModule, "loadoctree"),
    ProcessLoadOctree::create,
    "Loads octree into m_mesh");

ProcessLoadOctree::ProcessLoadOctree(MeshSharedPtr m) : ProcessModule(m)
{
    m_config["mindel"] =
        ConfigOption(false, "0", "mindelta.");
    m_config["maxdel"] =
        ConfigOption(false, "0", "mindelta.");
    m_config["eps"] =
        ConfigOption(false, "0", "mindelta.");
    m_config["refinement"] =
        ConfigOption(false, "", "mindelta.");
    m_config["writeoctree"] =
        ConfigOption(true, "0", "dump octree as xml mesh");
}

ProcessLoadOctree::~ProcessLoadOctree()
{
}

void ProcessLoadOctree::Process()
{
    NekDouble minDelta, maxDelta, eps;

    minDelta = m_config["mindel"].as<NekDouble>();
    maxDelta = m_config["maxdel"].as<NekDouble>();
    eps = m_config["eps"].as<NekDouble>();

    if (m_mesh->m_verbose)
    {
        cout << endl << "Loading Octree with parameters:" << endl;
        cout << "\tmin delta: " << minDelta << endl
             << "\tmax delta: " << maxDelta << endl
             << "\tesp: " << eps << endl << endl;
    }

    ASSERTL0(minDelta > 0 && maxDelta > 0 && eps > 0, "invalid parameters");

    m_mesh->m_octree = MemoryManager<Octree>::AllocateSharedPtr(m_mesh);

    m_mesh->m_octree->SetParameters(minDelta, maxDelta, eps);

    if(m_config["refinement"].beenSet)
    {
        m_mesh->m_octree->Refinement(m_config["refinement"].as<string>());
    }

    m_mesh->m_octree->Process();

    if(m_config["writeoctree"].beenSet)
    {
        m_mesh->m_octree->WriteOctree(m_config["writeoctree"].as<string>());
    }
}
}
}

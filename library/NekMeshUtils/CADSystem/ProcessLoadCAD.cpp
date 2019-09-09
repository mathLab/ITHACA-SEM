////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessLoadCAD.cpp
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
//  Description: Load CAD module
//
////////////////////////////////////////////////////////////////////////////////

#include "ProcessLoadCAD.h"
#include <NekMeshUtils/CADSystem/CADSystem.h>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

using namespace std;
namespace Nektar
{
namespace NekMeshUtils
{

ModuleKey ProcessLoadCAD::className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eProcessModule, "loadcad"),
    ProcessLoadCAD::create,
    "Loads cad into m_mesh");

ProcessLoadCAD::ProcessLoadCAD(MeshSharedPtr m) : ProcessModule(m)
{
    m_config["filename"] =
        ConfigOption(false, "", "Generate prisms on these surfs");
    m_config["2D"] =
        ConfigOption(true, "", "allow 2d loading");
    m_config["CFIMesh"] =
        ConfigOption(true, "", "specifies that the CAD can be multibody");
    m_config["NACA"] =
        ConfigOption(false, "", "naca domain");
    m_config["verbose"] =
        ConfigOption(true, "", "verbose output from cadsystem");
}

ProcessLoadCAD::~ProcessLoadCAD()
{
}

void ProcessLoadCAD::Process()
{
    string name = m_config["filename"].as<string>();

    if (m_mesh->m_verbose)
    {
        cout << "Loading CAD for " << name << endl;
    }

    string ext = boost::filesystem::extension(name);

    if(boost::iequals(ext,".fbm"))
    {
        m_mesh->m_cad = GetEngineFactory().CreateInstance("cfi",name);
    }
    else
    {
        m_mesh->m_cad = GetEngineFactory().CreateInstance("oce",name);
    }

    if(m_config["2D"].beenSet)
    {
        m_mesh->m_cad->Set2D();
    }

    if(m_config["NACA"].beenSet)
    {
        m_mesh->m_cad->SetNACA(m_config["NACA"].as<string>());
    }

    if(m_config["CFIMesh"].beenSet)
    {
        m_mesh->m_cad->SetCFIMesh();
    }

    if(m_config["verbose"].beenSet)
    {
        m_mesh->m_cad->SetVerbose();
    }

    ASSERTL0(m_mesh->m_cad->LoadCAD(), "Failed to load CAD");
}
}
}

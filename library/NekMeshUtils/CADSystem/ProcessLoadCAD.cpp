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
    m_config["cfimesh"] =
        ConfigOption(true, "0", "load existing mesh from cfi");
}

ProcessLoadCAD::~ProcessLoadCAD()
{
}

void ProcessLoadCAD::Process()
{
    m_mesh->m_CADId = m_config["filename"].as<string>();

    if (m_mesh->m_verbose)
    {
        cout << "Loading CAD for " << m_mesh->m_CADId << endl;
    }

    string ext = boost::filesystem::extension(m_mesh->m_CADId);

    if(boost::iequals(ext,".fbm"))
    {
        m_mesh->m_cad = GetEngineFactory().CreateInstance("cfi",m_mesh->m_CADId);
    }
    else
    {
        m_mesh->m_cad = GetEngineFactory().CreateInstance("oce",m_mesh->m_CADId);
    }

    ASSERTL0(m_mesh->m_cad->LoadCAD(), "Failed to load CAD");

    if(boost::iequals(ext,".fbm") && m_config["cfimesh"].beenSet)
    {
        m_mesh->LoadMeshFromCAD();
        //vertices already done
        ProcessEdges();
        ProcessFaces();
        ProcessElements();
        ProcessComposites();
        m_mesh->BuildComps();
        ProcessEdges();
        ProcessFaces();
        ProcessElements();
        ProcessComposites();
        m_mesh->ReconstructSurfaceAndCADInfo();
        ProcessVertices();
        ProcessEdges();
        ProcessFaces();
        ProcessElements();
        ProcessComposites();
    }

    m_mesh->m_hasCAD = true;

    if (m_mesh->m_verbose)
    {
        m_mesh->m_cad->Report();
    }
}
}
}

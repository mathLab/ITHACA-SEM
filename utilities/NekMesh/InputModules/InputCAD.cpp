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
//  Description: create mesh from cad using mesh utils
//
////////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/ParseUtils.hpp>

#include <boost/filesystem.hpp>

#include <NekMeshUtils/MeshElements/Element.h>

#include <NekMeshUtils/Octree/Octree.h>
#include <NekMeshUtils/SurfaceMeshing/SurfaceMesh.h>
#include <NekMeshUtils/VolumeMeshing/VolumeMesh.h>

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/Communication/CommSerial.h>

#include "InputCAD.h"

using namespace std;
using namespace Nektar::NekMeshUtils;

namespace Nektar
{
namespace Utilities
{

ModuleKey InputCAD::className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eInputModule, "mcf"),
    InputCAD::create,
    "Reads CAD geometry and will generate the mesh file.");

/**
 * @brief Set up InputCAD object.
 */
InputCAD::InputCAD(MeshSharedPtr m) : InputModule(m)
{
}

InputCAD::~InputCAD()
{
}

void InputCAD::Process()
{
    vector<string> filename;
    filename.push_back(m_config["infile"].as<string>());
    string fn = filename[0].substr(0, filename[0].find("."));

    LibUtilities::SessionReaderSharedPtr pSession =
        LibUtilities::SessionReader::CreateInstance(0, NULL, filename);

    // these parameters must be defined for any mesh generation to work
    pSession->LoadParameter("MinDelta", m_minDelta);
    pSession->LoadParameter("MaxDelta", m_maxDelta);
    pSession->LoadParameter("EPS", m_eps);
    pSession->LoadParameter("Order", m_order);
    //m_mesh->m_CADId = pSession->GetSolverInfo("CADFile");

    m_mesh->m_expDim   = 3;
    m_mesh->m_spaceDim = 3;
    m_mesh->m_nummode = m_order + 1;

    vector<ModuleSharedPtr> mods;

    mods.push_back(GetModuleFactory().CreateInstance(
        ModuleKey(eProcessModule, "loadcad"), m_mesh));
    mods.back()->RegisterConfig("filename",pSession->GetSolverInfo("CADFile"));

    mods.push_back(GetModuleFactory().CreateInstance(
        ModuleKey(eProcessModule, "loadoctree"), m_mesh));
    mods.back()->RegisterConfig("mindel",boost::lexical_cast<std::string>(m_minDelta));
    mods.back()->RegisterConfig("maxdel",boost::lexical_cast<std::string>(m_maxDelta));
    mods.back()->RegisterConfig("eps",boost::lexical_cast<std::string>(m_eps));

    if(pSession->DefinesSolverInfo("SourcePoints"))
    {
        ASSERTL0(boost::filesystem::exists(pSession->GetSolverInfo("SourcePoints").c_str()),
                 "sourcepoints file does not exist");
        mods.back()->RegisterConfig("sourcefile",pSession->GetSolverInfo("SourcePoints"));
        NekDouble sp;
        pSession->LoadParameter("SPSize", sp);
        mods.back()->RegisterConfig("sourcesize",boost::lexical_cast<std::string>(sp));
    }

    if (pSession->DefinesSolverInfo("UserDefinedSpacing"))
    {
        string udsName = pSession->GetSolverInfo("UserDefinedSpacing");
        ASSERTL0(boost::filesystem::exists(udsName.c_str()),
                 "UserDefinedSpacing file does not exist");

        mods.back()->RegisterConfig("udsfile",udsName);
    }

    if (pSession->DefinesSolverInfo("WriteOctree"))
    {
        mods.back()->RegisterConfig("writeoctree",fn + "_oct.xml");
    }

    //create surface mesh

    mods.push_back(GetModuleFactory().CreateInstance(
        ModuleKey(eProcessModule, "surfacemesh"), m_mesh));

    mods.push_back(GetModuleFactory().CreateInstance(
        ModuleKey(eProcessModule, "volumemesh"), m_mesh));

    mods.push_back(GetModuleFactory().CreateInstance(
        ModuleKey(eProcessModule, "hosurface"), m_mesh));
    if(pSession->DefinesSolverInfo("SurfaceOpt"))
    {
        mods.back()->RegisterConfig("opti","");
    }

    for(int i = 0; i < mods.size(); i++)
    {
        mods[i]->Process();
    }
}
}
}

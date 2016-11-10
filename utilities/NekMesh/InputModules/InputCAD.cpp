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

#include <tinyxml.h>

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

void InputCAD::ParseFile(string nm)
{
    vector<string> filename;
    filename.push_back(nm);
    LibUtilities::SessionReaderSharedPtr pSession =
        LibUtilities::SessionReader::CreateInstance(0, NULL, filename);

    ASSERTL0(pSession->DefinesElement("NEKTAR/MESHING"),"no meshing tag");
    ASSERTL0(pSession->DefinesElement("NEKTAR/MESHING/INFORMATION"),"no information tag");
    ASSERTL0(pSession->DefinesElement("NEKTAR/MESHING/PARAMETERS"),"no parameters tag");

    TiXmlElement *mcf = pSession->GetElement("NEKTAR/MESHING");

    TiXmlElement *info = mcf->FirstChildElement("INFORMATION");
    TiXmlElement *I = info->FirstChildElement("I");
    map<string,string> information;
    while(I)
    {
        string tmp1,tmp2;
        I->QueryStringAttribute("PROPERTY",&tmp1);
        I->QueryStringAttribute("VALUE",&tmp2);
        information[tmp1] = tmp2;
        I = I->NextSiblingElement("I");
    }

    TiXmlElement *param = mcf->FirstChildElement("PARAMETERS");
    TiXmlElement *P = param->FirstChildElement("P");
    map<string,string> parameters;
    while(P)
    {
        string tmp1,tmp2;
        P->QueryStringAttribute("PARAM",&tmp1);
        P->QueryStringAttribute("VALUE",&tmp2);
        parameters[tmp1] = tmp2;
        P = P->NextSiblingElement("P");
    }

    set<string> boolparameters;

    if(pSession->DefinesElement("NEKTAR/MESHING/BOOLPARAMETERS"))
    {
        TiXmlElement *bparam = mcf->FirstChildElement("BOOLPARAMETERS");
        TiXmlElement *BP = bparam->FirstChildElement("P");

        while(BP)
        {
            string tmp;
            BP->QueryStringAttribute("VALUE",&tmp);
            boolparameters.insert(tmp);
            BP = BP->NextSiblingElement("P");
        }
    }

    map<string,string>::iterator it;

    it = information.find("CADFile");
    ASSERTL0(it != information.end(),"no cadfile defined");
    m_cadfile = it->second;
    it = information.find("MeshType");
    ASSERTL0(it != information.end(),"no meshtype defined");
    m_makeBL = it->second == "BL";

    it = information.find("UDSFile");
    if(it != information.end())
    {
        m_udsfile = it->second;
        m_uds = true;
    }

    it = parameters.find("MinDelta");
    ASSERTL0(it != parameters.end(),"no mindelta defined");
    m_minDelta = it->second;
    it = parameters.find("MaxDelta");
    ASSERTL0(it != parameters.end(),"no maxdelta defined");
    m_maxDelta = it->second;
    it = parameters.find("EPS");
    ASSERTL0(it != parameters.end(),"no eps defined");
    m_eps = it->second;
    it = parameters.find("Order");
    ASSERTL0(it != parameters.end(),"no order defined");
    m_order = it->second;
    if(m_makeBL)
    {
        it = parameters.find("BLSurfs");
        ASSERTL0(it != parameters.end(), "no blsurfs defined");
        m_blsurfs = it->second;
        it = parameters.find("BLThick");
        ASSERTL0(it != parameters.end(), "no blthick defined");
        m_blthick = it->second;
        it = parameters.find("BLLayers");
        ASSERTL0(it != parameters.end(), "no bllayer defined");
        m_bllayers = it->second;
        it = parameters.find("BLProg");
        ASSERTL0(it != parameters.end(), "no blprog defined");
        m_blprog = it->second;
    }

    set<string>::iterator sit;
    sit = boolparameters.find("SurfOpti");
    m_surfopti = sit != boolparameters.end();
    sit = boolparameters.find("WriteOctree");
    m_woct = sit != boolparameters.end();
}

void InputCAD::Process()
{
    ParseFile(m_config["infile"].as<string>());

    m_mesh->m_expDim   = 3;
    m_mesh->m_spaceDim = 3;
    m_mesh->m_nummode = boost::lexical_cast<int>(m_order) + 1;

    vector<ModuleSharedPtr> mods;

    ////**** CAD ****////
    mods.push_back(GetModuleFactory().CreateInstance(
        ModuleKey(eProcessModule, "loadcad"), m_mesh));
    mods.back()->RegisterConfig("filename", m_cadfile);

    ////**** Octree ****////
    mods.push_back(GetModuleFactory().CreateInstance(
        ModuleKey(eProcessModule, "loadoctree"), m_mesh));
    mods.back()->RegisterConfig("mindel", m_minDelta);
    mods.back()->RegisterConfig("maxdel", m_maxDelta);
    mods.back()->RegisterConfig("eps", m_eps);
    if(m_uds)
    {
        mods.back()->RegisterConfig("udsfile", m_udsfile);
    }
    if(m_woct)
    {
        mods.back()->RegisterConfig("writeoctree", "");
    }

    ////**** SurfaceMesh ****////
    mods.push_back(GetModuleFactory().CreateInstance(
        ModuleKey(eProcessModule, "surfacemesh"), m_mesh));

    ////**** VolumeMesh ****////
    mods.push_back(GetModuleFactory().CreateInstance(
        ModuleKey(eProcessModule, "volumemesh"), m_mesh));
    if(m_makeBL)
    {
        mods.back()->RegisterConfig("blsurfs",m_blsurfs);
        mods.back()->RegisterConfig("blthick",m_blthick);
        mods.back()->RegisterConfig("bllayers",m_bllayers);
        mods.back()->RegisterConfig("blprog",m_blprog);
    }

    ////**** HOSurface ****////
    mods.push_back(GetModuleFactory().CreateInstance(
        ModuleKey(eProcessModule, "hosurface"), m_mesh));
    if(m_surfopti)
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

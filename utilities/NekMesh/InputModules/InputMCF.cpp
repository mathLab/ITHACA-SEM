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

#include <boost/thread.hpp>

#include <tinyxml.h>

#include "InputMCF.h"

using namespace std;
using namespace Nektar::NekMeshUtils;

namespace Nektar
{
namespace Utilities
{

ModuleKey InputMCF::className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eInputModule, "mcf"), InputMCF::create,
    "Reads mesh configuration and will generate the mesh file.");

/**
 * @brief Set up InputCAD object.
 */
InputMCF::InputMCF(MeshSharedPtr m) : InputModule(m)
{
}

InputMCF::~InputMCF()
{
}

void InputMCF::ParseFile(string nm)
{
    vector<string> filename;
    filename.push_back(nm);
    LibUtilities::SessionReaderSharedPtr pSession =
        LibUtilities::SessionReader::CreateInstance(0, NULL, filename);

    ASSERTL0(pSession->DefinesElement("NEKTAR/MESHING"), "no meshing tag");
    ASSERTL0(pSession->DefinesElement("NEKTAR/MESHING/INFORMATION"),
             "no information tag");
    ASSERTL0(pSession->DefinesElement("NEKTAR/MESHING/PARAMETERS"),
             "no parameters tag");

    TiXmlElement *mcf = pSession->GetElement("NEKTAR/MESHING");

    TiXmlElement *info = mcf->FirstChildElement("INFORMATION");
    TiXmlElement *I    = info->FirstChildElement("I");
    map<string, string> information;
    while (I)
    {
        string tmp1, tmp2;
        I->QueryStringAttribute("PROPERTY", &tmp1);
        I->QueryStringAttribute("VALUE", &tmp2);
        information[tmp1] = tmp2;
        I                 = I->NextSiblingElement("I");
    }

    TiXmlElement *param = mcf->FirstChildElement("PARAMETERS");
    TiXmlElement *P     = param->FirstChildElement("P");
    map<string, string> parameters;
    while (P)
    {
        string tmp1, tmp2;
        P->QueryStringAttribute("PARAM", &tmp1);
        P->QueryStringAttribute("VALUE", &tmp2);
        parameters[tmp1] = tmp2;
        P                = P->NextSiblingElement("P");
    }

    set<string> boolparameters;

    if (pSession->DefinesElement("NEKTAR/MESHING/BOOLPARAMETERS"))
    {
        TiXmlElement *bparam = mcf->FirstChildElement("BOOLPARAMETERS");
        TiXmlElement *BP     = bparam->FirstChildElement("P");

        while (BP)
        {
            string tmp;
            BP->QueryStringAttribute("VALUE", &tmp);
            boolparameters.insert(tmp);
            BP = BP->NextSiblingElement("P");
        }
    }

    set<string> refinement;
    if(pSession->DefinesElement("NEKTAR/MESHING/REFINEMENT"))
    {
        TiXmlElement *refine = mcf->FirstChildElement("REFINEMENT");
        TiXmlElement *L     = refine->FirstChildElement("LINE");

        while (L)
        {
            stringstream ss;
            TiXmlElement *T = L->FirstChildElement("X1");
            ss << T->GetText() << ",";
            T = L->FirstChildElement("Y1");
            ss << T->GetText() << ",";
            T = L->FirstChildElement("Z1");
            ss << T->GetText() << ",";
            T = L->FirstChildElement("X2");
            ss << T->GetText() << ",";
            T = L->FirstChildElement("Y2");
            ss << T->GetText() << ",";
            T = L->FirstChildElement("Z2");
            ss << T->GetText() << ",";
            T = L->FirstChildElement("R");
            ss << T->GetText() << ",";
            T = L->FirstChildElement("D");
            ss << T->GetText();

            refinement.insert(ss.str());

            L = L->NextSiblingElement("LINE");
        }
    }

    map<string,string>::iterator it;

    it = information.find("CADFile");
    ASSERTL0(it != information.end(), "no cadfile defined");
    m_cadfile = it->second;

    it = information.find("MeshType");
    ASSERTL0(it != information.end(), "no meshtype defined");
    m_makeBL = it->second == "3DBndLayer";
    m_2D = it->second == "2D";
    if (it->second == "2DBndLayer")
    {
        m_makeBL = true;
        m_2D = true;
    }
    if(!m_makeBL && !m_2D)
    {
        ASSERTL0(it->second == "3D", "unsure on MeshType")
    }


    it = parameters.find("MinDelta");
    ASSERTL0(it != parameters.end(), "no mindelta defined");
    m_minDelta = it->second;

    it = parameters.find("MaxDelta");
    ASSERTL0(it != parameters.end(), "no maxdelta defined");
    m_maxDelta = it->second;

    it = parameters.find("EPS");
    ASSERTL0(it != parameters.end(), "no eps defined");
    m_eps = it->second;

    it = parameters.find("Order");
    ASSERTL0(it != parameters.end(), "no order defined");
    m_order = it->second;

    if (m_makeBL)
    {
        it = parameters.find("BndLayerSurfaces");
        ASSERTL0(it != parameters.end(), "no BndLayersurfs defined");
        m_blsurfs = it->second;

        it = parameters.find("BndLayerThickness");
        ASSERTL0(it != parameters.end(), "no BndLayerthick defined");
        m_blthick = it->second;

        it = parameters.find("BndLayerLayers");
        m_splitBL = it != parameters.end();
        if(m_splitBL)
        {
            m_bllayers = it->second;
            it = parameters.find("BndLayerProgression");
            m_blprog = it != parameters.end() ? it->second : "2.0";
        }
    }

    m_naca = false;
    if(m_2D && m_cadfile.find('.') == std::string::npos)
    {
        m_naca = true;

        stringstream ss;
        it = parameters.find("Xmin");
        ASSERTL0(it != parameters.end(), "no xmin defined");
        ss << it->second << ",";
        it = parameters.find("Ymin");
        ASSERTL0(it != parameters.end(), "no ymin defined");
        ss << it->second << ",";
        it = parameters.find("Xmax");
        ASSERTL0(it != parameters.end(), "no xmax defined");
        ss << it->second << ",";
        it = parameters.find("Ymax");
        ASSERTL0(it != parameters.end(), "no zmax defined");
        ss << it->second << ",";
        it = parameters.find("AOA");
        ASSERTL0(it != parameters.end(), "no aoa defined");
        ss << it->second;

        m_nacadomain = ss.str();
    }

    set<string>::iterator sit;
    sit        = boolparameters.find("SurfaceOptimiser");
    m_surfopti = sit != boolparameters.end();
    sit        = boolparameters.find("WriteOctree");
    m_woct     = sit != boolparameters.end();
    sit        = boolparameters.find("VariationalOptimiser");
    m_varopti  = sit != boolparameters.end();

    m_refine = refinement.size() > 0;
    if(m_refine)
    {
        stringstream ss;
        for(sit = refinement.begin(); sit != refinement.end(); sit++)
        {
            ss << *sit;
            ss << ":";
        }
        m_refinement = ss.str();
        m_refinement.erase(m_refinement.end()-1);
    }
}

void InputMCF::Process()
{
    ParseFile(m_config["infile"].as<string>());

    m_mesh->m_expDim   = 3;
    m_mesh->m_spaceDim = 3;
    m_mesh->m_nummode  = boost::lexical_cast<int>(m_order) + 1;

    vector<ModuleSharedPtr> mods;

    ////**** CAD ****////
    mods.push_back(GetModuleFactory().CreateInstance(
        ModuleKey(eProcessModule, "loadcad"), m_mesh));
    mods.back()->RegisterConfig("filename", m_cadfile);

    if(m_2D)
    {
        mods.back()->RegisterConfig("2D","");
    }
    if(m_naca)
    {
        mods.back()->RegisterConfig("NACA",m_nacadomain);
    }

    ////**** Octree ****////
    mods.push_back(GetModuleFactory().CreateInstance(
        ModuleKey(eProcessModule, "loadoctree"), m_mesh));
    mods.back()->RegisterConfig("mindel", m_minDelta);
    mods.back()->RegisterConfig("maxdel", m_maxDelta);
    mods.back()->RegisterConfig("eps", m_eps);
    if (m_refine)
    {
        mods.back()->RegisterConfig("refinement", m_refinement);
    }
    if (m_woct)
    {
        mods.back()->RegisterConfig("writeoctree", "");
    }

    if(m_2D)
    {
        m_mesh->m_expDim = 2;
        m_mesh->m_spaceDim = 2;
        mods.push_back(GetModuleFactory().CreateInstance(
            ModuleKey(eProcessModule, "2dgenerator"), m_mesh));
        if (m_makeBL)
        {
            mods.back()->RegisterConfig("blcurves", m_blsurfs);
            mods.back()->RegisterConfig("blthick", m_blthick);
        }
    }
    else
    {
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
    }

    ////**** HOSurface ****////
    mods.push_back(GetModuleFactory().CreateInstance(
        ModuleKey(eProcessModule, "hosurface"), m_mesh));
    if (m_surfopti)
    {
        mods.back()->RegisterConfig("opti", "");
    }

    ////*** VARIATIONAL OPTIMISATION ****////
    if(m_varopti)
    {
        unsigned int np = boost::thread::physical_concurrency();
        if(m_mesh->m_verbose)
        {
            cout << "Detecting 4 cores, will attempt to run in parrallel" << endl;
        }
        mods.push_back(GetModuleFactory().CreateInstance(
            ModuleKey(eProcessModule, "varopti"), m_mesh));
        mods.back()->RegisterConfig("hyperelastic","");
        mods.back()->RegisterConfig("maxiter","10");
        mods.back()->RegisterConfig("numthreads",boost::lexical_cast<string>(np));
    }

    ////**** SPLIT BL ****////
    if(m_splitBL)
    {
        mods.push_back(GetModuleFactory().CreateInstance(
            ModuleKey(eProcessModule, "bl"), m_mesh));
        mods.back()->RegisterConfig("layers",m_bllayers);
        mods.back()->RegisterConfig("surf",m_blsurfs);
        mods.back()->RegisterConfig("nq",boost::lexical_cast<string>(m_mesh->m_nummode));
        mods.back()->RegisterConfig("r",m_blprog);
    }

    for(int i = 0; i < mods.size(); i++)
    {
        mods[i]->SetDefaults();
        mods[i]->Process();
    }
}
}
}

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
#include <NekMesh/CADSystem/CADCurve.h>

#include <boost/thread.hpp>
#include <boost/algorithm/string.hpp>

#include <tinyxml.h>

#include "InputMCF.h"

using namespace std;
using namespace Nektar::NekMesh;

namespace Nektar
{
namespace NekMesh
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

    char *prgname = (char*)"NekMesh";
    LibUtilities::SessionReaderSharedPtr pSession =
        LibUtilities::SessionReader::CreateInstance(1, &prgname, filename);
    pSession->InitSession();

    auto comm = pSession->GetComm();
    if (comm->GetType().find("MPI") != std::string::npos)
    {
        m_mesh->m_comm = comm;
    }

    if (!pSession->DefinesElement("NEKTAR/MESHING"))
    {
        m_log(FATAL) << "Input file is missing a MESHING tag." << endl;
    }

    if (!pSession->DefinesElement("NEKTAR/MESHING/INFORMATION"))
    {
        m_log(FATAL) << "Input file is missing an INFORMATION tag." << endl;
    }

    if (!pSession->DefinesElement("NEKTAR/MESHING/PARAMETERS"))
    {
        m_log(FATAL) << "Input file is missing a PARAMETERS tag." << endl;
    }

    TiXmlElement *mcf = pSession->GetElement("NEKTAR/MESHING");

    // Save MESHING tag as provenance information.
    std::stringstream ss;
    ss << *mcf;
    m_mesh->m_metadata["XML_NekMeshMCF"] = ss.str();

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
    if (pSession->DefinesElement("NEKTAR/MESHING/REFINEMENT"))
    {
        TiXmlElement *refine = mcf->FirstChildElement("REFINEMENT");
        TiXmlElement *L      = refine->FirstChildElement("LINE");

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

    set<string> periodic;
    if (pSession->DefinesElement("NEKTAR/MESHING/PERIODIC"))
    {
        TiXmlElement *per  = mcf->FirstChildElement("PERIODIC");
        TiXmlElement *pair = per->FirstChildElement("P");

        while (pair)
        {
            string tmp;
            pair->QueryStringAttribute("PAIR", &tmp);
            periodic.insert(tmp);
            pair = pair->NextSiblingElement("P");
        }
    }

    if (pSession->DefinesElement("NEKTAR/MESHING/VOIDPOINTS"))
    {
        TiXmlElement *vpts = mcf->FirstChildElement("VOIDPOINTS");
        TiXmlElement *v = vpts->FirstChildElement("V");
        stringstream ss;
        while (v)
        {
            std::string tmp = v->GetText();
            boost::trim(tmp);
            ss << tmp << ";";
            v = v->NextSiblingElement("V");
        }
        m_voidPts = ss.str();
        m_voidPts.pop_back();
    }

    auto it = information.find("CADFile");

    if (it == information.end())
    {
        m_log(FATAL) << "No 'CADFile' property found in the INFORMATION "
                     << "section." << endl;
    }
    m_cadfile = it->second;

    it = information.find("MeshType");
    if (it == information.end())
    {
        m_log(FATAL) << "No 'MeshType' property found in the INFORMATION "
                     << "section." << endl;
    }

    m_makeBL   = it->second == "3DBndLayer";
    m_2D       = it->second == "2D";
    m_manifold = it->second == "Manifold";

    if (it->second == "2DBndLayer")
    {
        m_makeBL = true;
        m_2D     = true;
    }

    if (!m_makeBL && !m_2D && !m_manifold && it->second != "3D")
    {
        if (it == information.end())
        {
            m_log(FATAL) << "Invalid 'MeshType' property found in the "
                         << "INFORMATION section: should be one of 2D, 3D, "
                         << "2DBndLayer, 3DBndLayer or Manifold." << endl;
        }
    }

    it = parameters.find("MinDelta");
    if (it == parameters.end())
    {
        m_log(FATAL) << "No 'MinDelta' parameter found in the PARAMETERS "
                     << "section." << endl;
    }
    m_minDelta = it->second;

    it = parameters.find("MaxDelta");
    if (it == parameters.end())
    {
        m_log(FATAL) << "No 'MaxDelta' parameter found in the PARAMETERS "
                     << "section." << endl;
    }
    m_maxDelta = it->second;

    it = parameters.find("EPS");
    if (it == parameters.end())
    {
        m_log(FATAL) << "No 'EPS' parameter found in the PARAMETERS "
                     << "section." << endl;
    }
    m_eps = it->second;

    it = parameters.find("Order");
    if (it == parameters.end())
    {
        m_log(FATAL) << "No 'Order' parameter found in the PARAMETERS "
                     << "section." << endl;
    }
    m_order = it->second;

    if (m_makeBL)
    {
        it = parameters.find("BndLayerSurfaces");
        if (it == parameters.end())
        {
            m_log(FATAL) << "No 'BndLayerSurfaces' parameter found in the "
                         << "PARAMETERS section." << endl;
        }
        m_blsurfs = it->second;

        it = parameters.find("BndLayerThickness");
        if (it == parameters.end())
        {
            m_log(FATAL) << "No 'BndLayerThickness' parameter found in the "
                         << "PARAMETERS section." << endl;
        }
        m_blthick = it->second;

        it        = parameters.find("BndLayerLayers");
        m_splitBL = it != parameters.end();
        if (m_splitBL)
        {
            m_bllayers = it->second;
            it         = parameters.find("BndLayerProgression");
            m_blprog   = it != parameters.end() ? it->second : "2.0";
        }

        it = parameters.find("BndLayerAdjustment");
        if (it != parameters.end())
        {
            m_adjust     = true;
            m_adjustment = it->second;
        }
        else
        {
            m_adjust = false;
        }

        it = parameters.find("SpaceOutBndLayer");
        if (it != parameters.end())
        {
            m_spaceoutbl    = true;
            m_spaceoutblthr = it->second;

            it = parameters.find("NoSpaceOutSurf");
            if (it != parameters.end())
            {
                m_nospaceoutsurf = it->second;
            }
        }
        else
        {
            m_spaceoutbl = false;
        }
    }
    else
    {
        m_splitBL = false;
    }

    m_naca = false;
    if (m_2D && m_cadfile.find('.') == std::string::npos)
    {
        m_naca = true;

        stringstream ss;
        it = parameters.find("Xmin");
        if (it == parameters.end())
        {
            m_log(FATAL) << "No 'Xmin' parameter found in the "
                         << "PARAMETERS section." << endl;
        }
        ss << it->second << ",";
        it = parameters.find("Ymin");
        if (it == parameters.end())
        {
            m_log(FATAL) << "No 'Ymin' parameter found in the "
                         << "PARAMETERS section." << endl;
        }
        ss << it->second << ",";
        it = parameters.find("Xmax");
        if (it == parameters.end())
        {
            m_log(FATAL) << "No 'Xmax' parameter found in the "
                         << "PARAMETERS section." << endl;
        }
        ss << it->second << ",";
        it = parameters.find("Ymax");
        if (it == parameters.end())
        {
            m_log(FATAL) << "No 'Ymax' parameter found in the "
                         << "PARAMETERS section." << endl;
        }
        ss << it->second << ",";
        it = parameters.find("AOA");
        if (it == parameters.end())
        {
            m_log(FATAL) << "No 'AOA' parameter found in the "
                         << "PARAMETERS section." << endl;
        }
        ss << it->second;

        m_nacadomain = ss.str();
    }

    auto sit    = boolparameters.find("SurfaceOptimiser");
    m_surfopti  = sit != boolparameters.end();
    sit         = boolparameters.find("WriteOctree");
    m_woct      = sit != boolparameters.end();
    sit         = boolparameters.find("VariationalOptimiser");
    m_varopti   = sit != boolparameters.end();
    sit         = boolparameters.find("BndLayerAdjustEverywhere");
    m_adjustall = sit != boolparameters.end();
    sit         = boolparameters.find("SmoothBndLayer");
    m_smoothbl  = sit != boolparameters.end();

    m_refine = refinement.size() > 0;
    if (m_refine)
    {
        stringstream ss;
        for (sit = refinement.begin(); sit != refinement.end(); sit++)
        {
            ss << *sit;
            ss << ":";
        }
        m_refinement = ss.str();
        m_refinement.erase(m_refinement.end() - 1);
    }

    if (periodic.size() > 0)
    {
        stringstream ss;
        for (sit = periodic.begin(); sit != periodic.end(); ++sit)
        {
            ss << *sit;
            ss << ":";
        }
        m_periodic = ss.str();
        m_periodic.erase(m_periodic.end() - 1);
    }
}

void InputMCF::Process()
{
    m_log(VERBOSE) << "Reading mesh configuration file '"
                   << m_config["infile"].as<string>() << "'" << endl;

    ParseFile(m_config["infile"].as<string>());

    m_mesh->m_expDim   = 3;
    m_mesh->m_spaceDim = 3;
    m_mesh->m_nummode  = boost::lexical_cast<int>(m_order) + 1;

    ModuleSharedPtr module;

    ////**** CAD ****////
    module = GetModuleFactory().CreateInstance(
        ModuleKey(eProcessModule, "loadcad"), m_mesh);
    module->SetLogger(m_log);
    module->RegisterConfig("filename", m_cadfile);
    module->RegisterConfig("voidpoints", m_voidPts);
    if (m_2D)
    {
        module->RegisterConfig("2D", "");
    }
    if (m_naca)
    {
        module->RegisterConfig("NACA", m_nacadomain);
    }

    module->SetDefaults();
    module->Process();

    ////**** OCTREE ****////
    module = GetModuleFactory().CreateInstance(
        ModuleKey(eProcessModule, "loadoctree"), m_mesh);
    module->SetLogger(m_log);
    module->RegisterConfig("mindel", m_minDelta);
    module->RegisterConfig("maxdel", m_maxDelta);
    module->RegisterConfig("eps", m_eps);
    if (m_refine)
    {
        module->RegisterConfig("refinement", m_refinement);
    }
    if (m_woct)
    {
        module->RegisterConfig("writeoctree", "");
    }

    module->SetDefaults();
    module->Process();

    ////**** LINEAR MESHING ****////
    if (m_2D)
    {
        ////**** 2DGenerator ****////
        m_mesh->m_expDim   = 2;
        m_mesh->m_spaceDim = 2;
        module             = GetModuleFactory().CreateInstance(
            ModuleKey(eProcessModule, "2dgenerator"), m_mesh);

        module->SetLogger(m_log);

        if (m_makeBL)
        {
            module->RegisterConfig("blcurves", m_blsurfs);
            module->RegisterConfig("blthick", m_blthick);

            if (m_adjust)
            {
                module->RegisterConfig("bltadjust", m_adjustment);

                if (m_adjustall)
                {
                    module->RegisterConfig("adjustblteverywhere", "");
                }
            }

            if (m_smoothbl)
            {
                module->RegisterConfig("smoothbl", "");
            }

            if (m_spaceoutbl)
            {
                module->RegisterConfig("spaceoutbl", m_spaceoutblthr);
                module->RegisterConfig("nospaceoutsurf", m_nospaceoutsurf);
            }
        }
        if (m_periodic.size())
        {
            module->RegisterConfig("periodic", m_periodic);
        }

        try
        {
            module->SetDefaults();
            module->Process();
        }
        catch (runtime_error &e)
        {
            m_log(WARNING) << "2D linear mesh generator failed with message:"
                           << endl;
            m_log(WARNING) << e.what() << endl;
            m_log(FATAL)   << "No mesh file has been created." << endl;
        }
    }
    else
    {
        ////**** SurfaceMesh ****////
        module = GetModuleFactory().CreateInstance(
            ModuleKey(eProcessModule, "surfacemesh"), m_mesh);

        try
        {
            module->SetLogger(m_log);
            module->SetDefaults();
            module->Process();
        }
        catch (runtime_error &e)
        {
            m_log(WARNING) << "Surface meshing has failed with message:"
                           << endl;
            m_log(WARNING) << e.what() << endl;
            m_log(WARNING) << "Any surfaces which were successfully meshed will"
                           << " be written as a manifold mesh."
                           << endl;

            m_mesh->m_expDim = 2;

            ProcessVertices();
            ProcessEdges();
            ProcessFaces();
            ProcessElements();
            ProcessComposites();
            return;
        }

        if (m_manifold)
        {
            // Don't want to volume mesh.
            m_mesh->m_expDim = 2;
        }
        else
        {
            ////**** VolumeMesh ****////
            module = GetModuleFactory().CreateInstance(
                ModuleKey(eProcessModule, "volumemesh"), m_mesh);

            module->SetLogger(m_log);

            if (m_makeBL)
            {
                module->RegisterConfig("blsurfs", m_blsurfs);
                module->RegisterConfig("blthick", m_blthick);
                module->RegisterConfig("bllayers", m_bllayers);
                module->RegisterConfig("blprog", m_blprog);
            }

            try
            {
                module->SetDefaults();
                module->Process();
            }
            catch (runtime_error &e)
            {
                m_log(WARNING) << "Volume meshing has failed with message:"
                               << endl;
                m_log(WARNING) << e.what() << endl;
                m_log(WARNING) << "The linear surface mesh be written as a "
                               << "manifold mesh"
                               << endl;

                m_mesh->m_expDim = 2;
                m_mesh->m_element[3].clear();

                ProcessVertices();
                ProcessEdges();
                ProcessFaces();
                ProcessElements();
                ProcessComposites();

                return;
            }
        }
    }

    ////**** HOSurface ****////
    module = GetModuleFactory().CreateInstance(
        ModuleKey(eProcessModule, "hosurface"), m_mesh);

    module->SetLogger(m_log);

    if (m_surfopti)
    {
        module->RegisterConfig("opti", "");
    }

    try
    {
        module->SetDefaults();
        module->Process();
    }
    catch (runtime_error &e)
    {
        m_log(WARNING) << "High-order surface meshing has failed with message:"
                       << endl;
        m_log(WARNING) << e.what() << endl;
        m_log(WARNING) << "The mesh will be written as normal but the "
                       << "incomplete surface will remain faceted"
                       << endl;
        return;
    }

    ////*** VARIATIONAL OPTIMISATION ****////
    if (m_varopti)
    {
        unsigned int np = boost::thread::physical_concurrency();

        m_log(VERBOSE) << "Detecting " << np << " cores, will attempt to run "
                       << "in parallel" << endl;

        module = GetModuleFactory().CreateInstance(
            ModuleKey(eProcessModule, "varopti"), m_mesh);

        module->SetLogger(m_log);
        module->RegisterConfig("hyperelastic", "");
        module->RegisterConfig("maxiter", "10");
        module->RegisterConfig("numthreads", boost::lexical_cast<string>(np));

        try
        {
            module->SetDefaults();
            module->Process();
        }
        catch (runtime_error &e)
        {
            m_log(WARNING) << "Variational optimisation has failed with "
                           << "message:" << endl;
            m_log(WARNING) << e.what() << endl;
            m_log(WARNING) << "The mesh will be written as is, it may be "
                           << "invalid" << endl;
            return;
        }
    }

    ////**** SPLIT BL ****////
    if (m_splitBL)
    {
        module = GetModuleFactory().CreateInstance(
            ModuleKey(eProcessModule, "bl"), m_mesh);

        module->SetLogger(m_log);
        module->RegisterConfig("layers", m_bllayers);
        module->RegisterConfig("surf", m_blsurfs);
        module->RegisterConfig("nq",
                               boost::lexical_cast<string>(m_mesh->m_nummode));
        module->RegisterConfig("r", m_blprog);

        try
        {
            module->SetDefaults();
            module->Process();
        }
        catch (runtime_error &e)
        {
            m_log(WARNING) << "Boundary layer splitting has failed with "
                           << "message:" << endl;
            m_log(WARNING) << e.what() << endl;
            m_log(WARNING) << "The mesh will be written as is, it may be "
                           << "invalid" << endl;
            return;
        }
    }

    // apply surface labels
    for (auto &it : m_mesh->m_composite)
    {
        ElementSharedPtr el = it.second->m_items[0];
        if (el->m_parentCAD)
        {
            string name = el->m_parentCAD->GetName();
            if (name.size() > 0)
            {
                m_mesh->m_faceLabels.insert(
                    make_pair(el->GetTagList()[0], name));
            }
        }
    }
    ProcessComposites();

    ////**** Peralign ****////
    if (m_2D && m_periodic.size())
    {
        vector<string> lines;
        boost::split(lines, m_periodic, boost::is_any_of(":"));

        for (auto &il : lines)
        {
            module = GetModuleFactory().CreateInstance(
                ModuleKey(eProcessModule, "peralign"), m_mesh);

            vector<string> tmp(2);
            boost::split(tmp, il, boost::is_any_of(","));
            module->SetLogger(m_log);
            module->RegisterConfig("surf1", tmp[0]);
            module->RegisterConfig("surf2", tmp[1]);

            module->SetDefaults();
            module->Process();
        }
    }

    m_log(VERBOSE) << "Input mesh generation complete." << endl;
}
}
}

////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessHomogeneousPlane.cpp
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
//  Description: Extract a single plane of a 3DH1D field.
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>
using namespace std;

#include <boost/core/ignore_unused.hpp>

#include <LibUtilities/BasicUtils/SharedArray.hpp>

#include "ProcessHomogeneousPlane.h"

namespace Nektar
{
namespace FieldUtils
{

ModuleKey ProcessHomogeneousPlane::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "homplane"),
        ProcessHomogeneousPlane::create,
        "Extracts a plane from a 3DH1D expansion, requires planeid to be "
        "defined.");

ProcessHomogeneousPlane::ProcessHomogeneousPlane(FieldSharedPtr f)
    : ProcessModule(f)
{
    m_config["planeid"] = ConfigOption(false, "NotSet", "plane id to extract");
    m_config["wavespace"] =
        ConfigOption(true, "0", "Extract plane in Fourier space");
}

ProcessHomogeneousPlane::~ProcessHomogeneousPlane()
{
}

void ProcessHomogeneousPlane::Process(po::variables_map &vm)
{
    boost::ignore_unused(vm);

    ASSERTL0(m_f->m_numHomogeneousDir == 1,
             "ProcessHomogeneousPlane only works for Homogeneous1D.");
    m_f->m_numHomogeneousDir = 0;

    // Skip in case of empty partition
    if (m_f->m_exp[0]->GetNumElmts() == 0)
    {
        return;
    }

    ASSERTL0(m_config["planeid"].m_beenSet,
             "Missing parameter planeid for ProcessHomogeneousPlane");

    int planeid = m_config["planeid"].as<int>();
    int nfields = m_f->m_variables.size();

    int nstrips;
    m_f->m_session->LoadParameter("Strip_Z", nstrips, 1);

    // Look for correct plane (because of parallel case)
    int plane = -1;
    for (int i = 0; i < m_f->m_exp[0]->GetZIDs().size(); ++i)
    {
        if (m_f->m_exp[0]->GetZIDs()[i] == planeid)
        {
            plane = i;
        }
    }

    if (plane != -1)
    {
        for (int s = 0; s < nstrips; ++s)
        {
            for (int i = 0; i < nfields; ++i)
            {
                int n         = s * nfields + i;
                m_f->m_exp[n] = m_f->m_exp[n]->GetPlane(plane);

                if (m_config["wavespace"].as<bool>())
                {
                    m_f->m_exp[n]->BwdTrans(m_f->m_exp[n]->GetCoeffs(),
                                            m_f->m_exp[n]->UpdatePhys());
                }
                else
                {
                    m_f->m_exp[n]->FwdTrans_IterPerExp(m_f->m_exp[n]->GetPhys(),
                                            m_f->m_exp[n]->UpdateCoeffs());
                }
            }
        }

        // Create new SessionReader with RowComm. This is done because when
        //  using a module requiring m_f->m_declareExpansionAsContField and
        //  outputting to vtu/dat, a new ContField with equispaced points
        //  is created. Since creating ContFields require communication, we have
        //  to change m_session->m_comm to prevent mpi from hanging
        //  (RowComm will only be used when creating the new expansion,
        //   since in other places we use m_f->m_comm)
        std::vector<std::string> files;
        for (int i = 0; i < m_f->m_inputfiles["xml"].size(); ++i)
        {
            files.push_back(m_f->m_inputfiles["xml"][i]);
        }
        for (int j = 0; j < m_f->m_inputfiles["xml.gz"].size(); ++j)
        {
            files.push_back(m_f->m_inputfiles["xml.gz"][j]);
        }
        vector<string> cmdArgs;
        cmdArgs.push_back("FieldConvert");
        if (m_f->m_verbose)
        {
            cmdArgs.push_back("--verbose");
        }
        int argc          = cmdArgs.size();
        const char **argv = new const char *[argc];
        for (int i = 0; i < argc; ++i)
        {
            argv[i] = cmdArgs[i].c_str();
        }
        m_f->m_session = LibUtilities::SessionReader::CreateInstance(
            argc, (char **)argv, files, m_f->m_comm->GetRowComm());
        m_f->m_session->InitSession();
    }
    else
    {
        // Create empty expansion
        for (int s = 0; s < nstrips; ++s)
        {
            for (int i = 0; i < nfields; ++i)
            {
                int n = s * nfields + i;
                m_f->m_exp[n] =
                    MemoryManager<MultiRegions::ExpList>::AllocateSharedPtr();
            }
        }
    }
}
}
}

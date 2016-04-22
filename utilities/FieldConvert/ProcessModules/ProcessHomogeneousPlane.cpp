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
//  Description: Extract a single plane of a 3DH1D field.
//
////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <iostream>
using namespace std;

#include "ProcessHomogeneousPlane.h"

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/ParseUtils.hpp>

namespace Nektar
{
namespace Utilities
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
        ConfigOption(true, "NotSet", "Extract plane in Fourier space");
}

ProcessHomogeneousPlane::~ProcessHomogeneousPlane()
{
}

void ProcessHomogeneousPlane::Process(po::variables_map &vm)
{
    if (m_f->m_verbose)
    {
        cout << "ProcessHomogeneousPlane: Extracting plane..." << endl;
    }

    if ((m_f->m_fielddef[0]->m_numHomogeneousDir) != 1)
    {
        ASSERTL0(false,
                 "ProcessHomogeneousPlane only works for Homogeneous1D.");
    }

    ASSERTL0(m_config["planeid"].m_beenSet,
             "Missing parameter planeid for ProcessHomogeneousPlane");

    int planeid = m_config["planeid"].as<int>();
    int nfields = m_f->m_fielddef[0]->m_fields.size();

    int nstrips;
    m_f->m_session->LoadParameter("Strip_Z", nstrips, 1);

    for (int s = 0; s < nstrips; ++s)
    {
        for (int i = 0; i < nfields; ++i)
        {
            int n = s * nfields + i;
            m_f->m_exp[n] = m_f->m_exp[n]->GetPlane(planeid);

            if (m_config["wavespace"].m_beenSet)
            {
                m_f->m_exp[n]->BwdTrans(m_f->m_exp[n]->GetCoeffs(),
                                        m_f->m_exp[n]->UpdatePhys());
            }
            else
            {
                m_f->m_exp[n]->FwdTrans(m_f->m_exp[n]->GetPhys(),
                                        m_f->m_exp[n]->UpdateCoeffs());
            }
        }
    }
    std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef =
        m_f->m_exp[0]->GetFieldDefinitions();
    std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());

    for (int s = 0; s < nstrips; ++s)
    {
        for (int j = 0; j < nfields; ++j)
        {
            for (int i = 0; i < FieldDef.size() / nstrips; ++i)
            {
                int n = s * FieldDef.size() / nstrips + i;

                FieldDef[n]->m_fields.push_back(
                    m_f->m_fielddef[0]->m_fields[j]);
                m_f->m_exp[s * nfields + j]->AppendFieldData(FieldDef[n],
                                                             FieldData[n]);
            }
        }
    }

    m_f->m_fielddef = FieldDef;
    m_f->m_data     = FieldData;
}
}
}

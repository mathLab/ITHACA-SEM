////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessMeanMode.cpp
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
//  Description: Extract mean mode of 3DH1D field.
//
////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <iostream>
using namespace std;

#include "ProcessMeanMode.h"

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/ParseUtils.hpp>

namespace Nektar
{
namespace Utilities
{

ModuleKey ProcessMeanMode::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "meanmode"),
        ProcessMeanMode::create, "Extract mean mode from 3DH1D.");

ProcessMeanMode::ProcessMeanMode(FieldSharedPtr f) : ProcessModule(f)
{
}

ProcessMeanMode::~ProcessMeanMode()
{
}

void ProcessMeanMode::Process(po::variables_map &vm)
{
    if (m_f->m_verbose)
    {
        cout << "ProcessMeanMode: Extracting mean mode..." << endl;
    }

    if ((m_f->m_fielddef[0]->m_numHomogeneousDir) != 1)
    {
        ASSERTL0(false, "ProcessMeanMode only works for Homogeneous1D.");
    }
    
    if (m_f->m_fielddef[0]->m_homogeneousZIDs[0] != 0)
    {
        ASSERTL0(false, "ProcessMeanMode: mean mode not found.");
    }
    
    int nfields = m_f->m_fielddef[0]->m_fields.size();

    int nstrips;
    m_f->m_session->LoadParameter("Strip_Z",nstrips,1);

    for(int s = 0; s < nstrips; ++s)
    {
        for (int i = 0; i < nfields; ++i)
        {
            int n = s*nfields + i;
            m_f->m_exp[n] = m_f->m_exp[n]->GetPlane(0);
            m_f->m_exp[n]->BwdTrans(m_f->m_exp[n]->GetCoeffs(),
                                    m_f->m_exp[n]->UpdatePhys());
        }
    }
    std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef
        = m_f->m_exp[0]->GetFieldDefinitions();
    std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());

    for(int s = 0; s < nstrips; ++s)
    {
        for (int j = 0; j < nfields; ++j)
        {
            for (int i = 0; i < FieldDef.size()/nstrips; ++i)
            {
                int n = s * FieldDef.size()/nstrips + i;

                FieldDef[n]->m_fields.push_back(m_f->m_fielddef[0]->m_fields[j]);
                m_f->m_exp[s*nfields+j]->AppendFieldData(FieldDef[n], FieldData[n]);
            }
        }
    }

    m_f->m_fielddef = FieldDef;
    m_f->m_data     = FieldData;
}

}
}

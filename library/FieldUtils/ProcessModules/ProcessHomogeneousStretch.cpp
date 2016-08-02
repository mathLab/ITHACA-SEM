////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessHomogeneousStretch.cpp
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
//  Description: Stretch homogeneous direction by an integer factor.
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>
using namespace std;

#include "ProcessHomogeneousStretch.h"

#include <LibUtilities/BasicUtils/ParseUtils.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>

namespace Nektar
{
namespace FieldUtils
{

ModuleKey ProcessHomogeneousStretch::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "homstretch"),
        ProcessHomogeneousStretch::create,
        "Stretches homogeneous direction of a 3DH1D expansion, requires factor "
        "to be "
        "defined. The number of modes in the final expansion can be defined "
        "using"
        " --output-points-hom-z.");

ProcessHomogeneousStretch::ProcessHomogeneousStretch(FieldSharedPtr f)
    : ProcessModule(f)
{
    m_config["factor"] = ConfigOption(false, "NotSet", "stretch factor");
}

ProcessHomogeneousStretch::~ProcessHomogeneousStretch()
{
}

void ProcessHomogeneousStretch::Process(po::variables_map &vm)
{
    if (m_f->m_verbose)
    {
        if (m_f->m_comm->GetRank() == 0)
        {
            cout << "ProcessHomogeneousStretch: Stretching expansion..."
                 << endl;
        }
    }

    if ((m_f->m_fielddef[0]->m_numHomogeneousDir) != 1)
    {
        ASSERTL0(false,
                 "ProcessHomogeneousStretch only works for Homogeneous1D.");
    }

    ASSERTL0(m_config["factor"].m_beenSet,
             "Missing parameter factor for ProcessHomogeneousStretch");

    int factor  = m_config["factor"].as<int>();
    int nfields = m_f->m_fielddef[0]->m_fields.size();
    int nplanes = m_f->m_exp[0]->GetHomogeneousBasis()->GetZ().num_elements();

    ASSERTL0(factor > 1, "Parameter factor must be greater than 1.");

    int nstrips;
    m_f->m_session->LoadParameter("Strip_Z", nstrips, 1);

    for (int s = 0; s < nstrips; ++s)
    {
        for (int i = 0; i < nfields; ++i)
        {
            int n       = s * nfields + i;
            int ncoeffs = m_f->m_exp[n]->GetPlane(0)->GetNcoeffs();
            // Loop planes backwards so we can copy in place
            for (int p = nplanes - 1; p > 1; --p)
            {
                int newP = 2 * factor * (p / 2) + p % 2;
                if (newP < nplanes)
                {
                    Vmath::Vcopy(
                        ncoeffs, m_f->m_exp[n]->GetPlane(p)->GetCoeffs(), 1,
                        m_f->m_exp[n]->GetPlane(newP)->UpdateCoeffs(), 1);
                }
                Vmath::Zero(ncoeffs, m_f->m_exp[n]->GetPlane(p)->UpdateCoeffs(),
                            1);
            }

            m_f->m_exp[n]->BwdTrans(m_f->m_exp[n]->GetCoeffs(),
                                    m_f->m_exp[n]->UpdatePhys());
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
    for (int i = 0; i < FieldDef.size(); ++i)
    {
        FieldDef[i]->m_homogeneousLengths[0] =
            factor * m_f->m_fielddef[i]->m_homogeneousLengths[0];
    }
    m_f->m_fielddef = FieldDef;
    m_f->m_data     = FieldData;
}
}
}

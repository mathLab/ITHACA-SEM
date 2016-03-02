////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessNumModes.cpp
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
//  Description: Writes the number of modes in each direction.
//
////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <iostream>
using namespace std;

#include "ProcessNumModes.h"

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/ParseUtils.hpp>
#include <StdRegions/StdQuadExp.h>

namespace Nektar
{
namespace Utilities
{

ModuleKey ProcessNumModes::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "nummodes"), ProcessNumModes::create,
        "Computes number of modes in each direction for each element.");

ProcessNumModes::ProcessNumModes(FieldSharedPtr f) : ProcessModule(f)
{
}

ProcessNumModes::~ProcessNumModes()
{
}

void ProcessNumModes::Process(po::variables_map &vm)
{
    if (m_f->m_verbose)
    {
        cout << "ProcessNumModes: Calculating number of modes..." << endl;
    }

    int i, j, s;
    int expdim    = m_f->m_graph->GetMeshDimension();
    int nfields   = m_f->m_fielddef[0]->m_fields.size();
    int addfields = expdim;
    int npoints   = m_f->m_exp[0]->GetNpoints();
    Array<OneD, Array<OneD, NekDouble> > outfield(addfields);

    int nstrips;

    m_f->m_session->LoadParameter("Strip_Z", nstrips, 1);

    m_f->m_exp.resize(nfields * nstrips);

    for (i = 0; i < addfields; ++i)
    {
        outfield[i] = Array<OneD, NekDouble>(npoints);
    }

    vector<MultiRegions::ExpListSharedPtr> Exp(nstrips * addfields);

    int nExp, nq, offset;
    nExp = m_f->m_exp[0]->GetExpSize();

    for (int n = 0; n < nExp; n++)
    {
        offset = m_f->m_exp[0]->GetPhys_Offset(n);
        nq     = m_f->m_exp[0]->GetExp(n)->GetTotPoints();

        for (i = 0; i < expdim; i++)
        {
            int P = m_f->m_exp[0]->GetExp(n)->GetBasis(i)->GetNumModes();
            Array<OneD, NekDouble> result = outfield[i] + offset;
            Vmath::Fill(nq, 1.0 * P, result, 1);
        }
    }

    for (s = 0; s < nstrips; ++s)
    {
        for (i = 0; i < addfields; ++i)
        {
            int n = s * addfields + i;
            Exp[n] =
                m_f->AppendExpList(m_f->m_fielddef[0]->m_numHomogeneousDir);
            Exp[n]->UpdatePhys() = outfield[i];
            Exp[n]->FwdTrans_IterPerExp(outfield[i], Exp[n]->UpdateCoeffs());
        }
    }

    vector<MultiRegions::ExpListSharedPtr>::iterator it;
    for (s = 0; s < nstrips; ++s)
    {
        for (i = 0; i < addfields; ++i)
        {
            it = m_f->m_exp.begin() + s * (nfields + addfields) + nfields + i;
            m_f->m_exp.insert(it, Exp[s * addfields + i]);
        }
    }

    vector<string> outname;
    outname.push_back("P1");
    if (addfields >= 2)
    {
        outname.push_back("P2");
    }

    if (addfields == 3)
    {
        outname.push_back("P3");
    }

    std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef =
        m_f->m_exp[0]->GetFieldDefinitions();
    std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());

    // homogeneous strip variant
    for (s = 0; s < nstrips; ++s)
    {
        for (j = 0; j < nfields + addfields; ++j)
        {
            for (i = 0; i < FieldDef.size() / nstrips; ++i)
            {
                int n = s * FieldDef.size() / nstrips + i;

                if (j >= nfields)
                {
                    FieldDef[n]->m_fields.push_back(outname[j - nfields]);
                }
                else
                {
                    FieldDef[n]->m_fields.push_back(
                        m_f->m_fielddef[0]->m_fields[j]);
                }

                m_f->m_exp[s * (nfields + addfields) + j]->AppendFieldData(
                    FieldDef[n], FieldData[n]);
            }
        }
    }

    m_f->m_fielddef = FieldDef;
    m_f->m_data     = FieldData;
}
}
}

////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessGrad.cpp
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
//  Description: Computes gradient of fields.
//
////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <iostream>
using namespace std;

#include "ProcessGrad.h"

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/ParseUtils.hpp>

namespace Nektar
{
namespace Utilities
{

ModuleKey ProcessGrad::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "gradient"),
        ProcessGrad::create, "Computes gradient of fields.");

ProcessGrad::ProcessGrad(FieldSharedPtr f) : ProcessModule(f)
{
}

ProcessGrad::~ProcessGrad()
{
}

void ProcessGrad::Process(po::variables_map &vm)
{
    if (m_f->m_verbose)
    {
        cout << "ProcessGrad: Calculating gradients..." << endl;
    }

    int i, j;
    int expdim    = m_f->m_graph->GetMeshDimension();
    int spacedim  = m_f->m_fielddef[0]->m_numHomogeneousDir + expdim;
    int nfields   = m_f->m_fielddef[0]->m_fields.size();
    int addfields = nfields*spacedim;

    int npoints = m_f->m_exp[0]->GetNpoints();
    Array<OneD, Array<OneD, NekDouble> > grad(addfields);
    m_f->m_exp.resize(nfields+addfields);

    for (i = 0; i < addfields; ++i)
    {
        grad[i] = Array<OneD, NekDouble>(npoints);
    }

    // Calculate Gradient
    for (i = 0; i < nfields; ++i)
    {
        for (j = 0; j < spacedim; ++j)
        {
            m_f->m_exp[i]->PhysDeriv(MultiRegions::DirCartesianMap[j],
                                     m_f->m_exp[i]->GetPhys(),
                                     grad[i*spacedim+j]);
        }
    }

    for (i = 0; i < addfields; ++i)
    {
        m_f->m_exp[nfields + i] = m_f->AppendExpList(m_f->m_fielddef[0]->m_numHomogeneousDir);
        m_f->m_exp[nfields + i]->UpdatePhys() = grad[i];
        m_f->m_exp[nfields + i]->FwdTrans_IterPerExp(grad[i],
                            m_f->m_exp[nfields + i]->UpdateCoeffs());
    }

    vector<string > outname;
    for (i = 0; i<nfields; ++i)
    {
        if(spacedim == 1)
        {
            outname.push_back(m_f->m_fielddef[0]->m_fields[i]+"_x");
        }
        else if (spacedim == 2)
        {
            outname.push_back(m_f->m_fielddef[0]->m_fields[i]+"_x");
            outname.push_back(m_f->m_fielddef[0]->m_fields[i]+"_y");
        }
        else if (spacedim == 3)
        {
            outname.push_back(m_f->m_fielddef[0]->m_fields[i]+"_x");
            outname.push_back(m_f->m_fielddef[0]->m_fields[i]+"_y");
            outname.push_back(m_f->m_fielddef[0]->m_fields[i]+"_z");
        }
    }

    std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef
        = m_f->m_exp[0]->GetFieldDefinitions();
    std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());

    for (j = 0; j < nfields + addfields; ++j)
    {
        for (i = 0; i < FieldDef.size(); ++i)
        {
            if (j >= nfields)
            {
                FieldDef[i]->m_fields.push_back(outname[j-nfields]);
            }
            else
            {
                FieldDef[i]->m_fields.push_back(m_f->m_fielddef[0]->m_fields[j]);
            }
            m_f->m_exp[j]->AppendFieldData(FieldDef[i], FieldData[i]);
        }
    }

    m_f->m_fielddef = FieldDef;
    m_f->m_data     = FieldData;
}

}
}

////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessAddCompositeID.cpp
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
//  Description: Add composite ID as a variable to the field.
//
////////////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include <iostream>
#include <string>

using namespace std;

#include "ProcessAddCompositeID.h"

#include <LibUtilities/BasicUtils/ParseUtils.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>

namespace Nektar
{
namespace FieldUtils
{

ModuleKey ProcessAddCompositeID::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "addcompositeid"),
        ProcessAddCompositeID::create,
        "Add a field which contains the composite ID of each element");

ProcessAddCompositeID::ProcessAddCompositeID(FieldSharedPtr f)
    : ProcessModule(f)
{
}

ProcessAddCompositeID::~ProcessAddCompositeID()
{
}

void ProcessAddCompositeID::Process(po::variables_map &vm)
{
    if (m_f->m_verbose)
    {
        if (m_f->m_comm->GetRank() == 0)
        {
            cout << "ProcessAddCompositeID: Adding composite ID as a new field"
                 << endl;
        }
    }

    int nfields           = 0;
    int NumHomogeneousDir = 0;
    MultiRegions::ExpListSharedPtr exp;

    if (m_f->m_fielddef.size())
    {
        nfields           = m_f->m_fielddef[0]->m_fields.size();
        NumHomogeneousDir = m_f->m_fielddef[0]->m_numHomogeneousDir;

        m_f->m_exp.resize(nfields + 1);
        exp = m_f->AppendExpList(NumHomogeneousDir, "Composite ID");

        m_f->m_exp[nfields] = exp;
    }
    else
    {
        exp = m_f->m_exp[0];
    }

    // Get Composites
    const SpatialDomains::CompositeMap CompositeMap =
        m_f->m_graph->GetComposites();

    SpatialDomains::CompositeMapConstIter it;
    NekDouble compid=0;

    // loop over elements
    for (int n = 0; n < exp->GetNumElmts(); ++n)
    {
        LocalRegions::ExpansionSharedPtr elmt = exp->GetExp(n);

        // loop over composite list and search for geomtry pointer in list
        for (it = CompositeMap.begin(); it != CompositeMap.end(); ++it)
        {
            if (find(it->second->begin(), it->second->end(), elmt->GetGeom()) !=
                it->second->end())
            {
                compid = it->first;
                break;
            }
        }

        WARNINGL0(it != CompositeMap.end(),
                  "Failed to find composite ID for element: " +
                      boost::lexical_cast<string>(n));

        // Fill element with the value of the index
        int npts = elmt->GetTotPoints();
        Array<OneD, NekDouble> tmp;
        Vmath::Fill(npts, compid,
                    tmp = exp->UpdatePhys() + exp->GetPhys_Offset(n), 1);
    }

    // forward transform
    exp->FwdTrans_IterPerExp(exp->GetPhys(), exp->UpdateCoeffs());

    std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef =
        m_f->m_exp[0]->GetFieldDefinitions();
    std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());

    // copy in previous fields if they exist.
    for (int i = 0; i < nfields; ++i)
    {
        for (int j = 0; j < FieldDef.size(); ++j)
        {
            FieldDef[j]->m_fields.push_back(m_f->m_fielddef[0]->m_fields[i]);
            m_f->m_exp[i]->AppendFieldData(FieldDef[j], FieldData[j]);
        }
    }

    // append composite id field
    for (int j = 0; j < FieldDef.size(); ++j)
    {
        FieldDef[j]->m_fields.push_back("compositeID");
        m_f->m_exp[nfields]->AppendFieldData(FieldDef[j], FieldData[j]);
    }

    m_f->m_fielddef = FieldDef;
    m_f->m_data     = FieldData;
}
}
}

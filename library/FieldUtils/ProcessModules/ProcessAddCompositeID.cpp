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

#include <boost/core/ignore_unused.hpp>

#include <LibUtilities/BasicUtils/SharedArray.hpp>

#include "ProcessAddCompositeID.h"

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
    boost::ignore_unused(vm);

    int nfields           = m_f->m_variables.size();
    m_f->m_variables.push_back("compositeID");
    // Skip in case of empty partition
    if (m_f->m_exp[0]->GetNumElmts() == 0)
    {
        return;
    }

    int NumHomogeneousDir = m_f->m_numHomogeneousDir;
    MultiRegions::ExpListSharedPtr exp;

    if (nfields)
    {
        m_f->m_exp.resize(nfields + 1);
        exp = m_f->AppendExpList(NumHomogeneousDir);

        m_f->m_exp[nfields] = exp;
    }
    else
    {
        exp = m_f->m_exp[0];
    }

    // Get Composites
    const SpatialDomains::CompositeMap CompositeMap =
        m_f->m_graph->GetComposites();

    int compid = -1;

    // loop over elements
    for (int n = 0; n < exp->GetNumElmts(); ++n)
    {
        LocalRegions::ExpansionSharedPtr elmt = exp->GetExp(n);

        // loop over composite list and search for geometry pointer in list
        for (auto &it : CompositeMap)
        {
            if (find(it.second->m_geomVec.begin(),
                     it.second->m_geomVec.end(), elmt->GetGeom()) !=
                it.second->m_geomVec.end())
            {
                compid = it.first;
                break;
            }
        }

        WARNINGL0(compid != -1,
                  "Failed to find composite ID for element: " +
                      boost::lexical_cast<string>(n));

        // Fill element with the value of the index
        int npts = elmt->GetTotPoints();
        Array<OneD, NekDouble> tmp;
        Vmath::Fill(npts, (NekDouble)compid,
                    tmp = exp->UpdatePhys() + exp->GetPhys_Offset(n), 1);
    }

    // forward transform
    exp->FwdTrans_IterPerExp(exp->GetPhys(), exp->UpdateCoeffs());
}
}
}

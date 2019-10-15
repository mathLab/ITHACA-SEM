////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessPrintFldNorms.cpp
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
//  Description: Prints the L2 and LInf norms of the field variables.
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>
using namespace std;

#include <boost/core/ignore_unused.hpp>

#include <LibUtilities/BasicUtils/SharedArray.hpp>

#include "ProcessPrintFldNorms.h"

namespace Nektar
{
namespace FieldUtils
{

ModuleKey ProcessPrintFldNorms::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "printfldnorms"),
        ProcessPrintFldNorms::create,
        "Print L2 and LInf norms to stdout.");

ProcessPrintFldNorms::ProcessPrintFldNorms(FieldSharedPtr f) : ProcessModule(f)
{
}

ProcessPrintFldNorms::~ProcessPrintFldNorms()
{
}

void ProcessPrintFldNorms::Process(po::variables_map &vm)
{
    boost::ignore_unused(vm);

    // Skip in case of empty partition
    if (m_f->m_exp[0]->GetNumElmts() == 0)
    {
        return;
    }

    // Evaluate norms and print
    for (int j = 0; j < m_f->m_exp.size(); ++j)
    {
        NekDouble L2   = m_f->m_exp[j]->L2(m_f->m_exp[j]->GetPhys());
        NekDouble LInf = m_f->m_exp[j]->Linf(m_f->m_exp[j]->GetPhys());

        if (m_f->m_comm->TreatAsRankZero())
        {
            cout << "L 2 error (variable "   << m_f->m_variables[j]
                 << ") : " << L2 << endl;
            cout << "L inf error (variable " << m_f->m_variables[j]
                 << ") : " << LInf << endl;
        }
    }
}
}
}

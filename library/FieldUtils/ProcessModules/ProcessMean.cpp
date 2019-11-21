////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessMean.cpp
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
//  Description: Compute the mean of each field.
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>
using namespace std;

#include <boost/core/ignore_unused.hpp>

#include <LibUtilities/BasicUtils/SharedArray.hpp>

#include "ProcessMean.h"

namespace Nektar
{
namespace FieldUtils
{

ModuleKey ProcessMean::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "mean"), ProcessMean::create,
        "compute the mean of each field over the domain.");

ProcessMean::ProcessMean(FieldSharedPtr f) : ProcessModule(f)
{
}

ProcessMean::~ProcessMean()
{
}

void ProcessMean::Process(po::variables_map &vm)
{
    boost::ignore_unused(vm);

    int nfields  = m_f->m_variables.size();
    int spacedim = m_f->m_graph->GetMeshDimension() + m_f->m_numHomogeneousDir;
    int npoints  = m_f->m_exp[0]->GetNpoints();

    // Calculate volume (or area)
    Array<OneD, NekDouble> ones(npoints, 1.0);
    NekDouble scale = m_f->m_exp[0]->Integral(ones);

    // Output volume
    string name[3] = {"length", "area", "volume"};
    cout << "Domain " << name[spacedim - 1] << " : " << scale << endl;

    // Calculate integral and mean of each field
    for (int i = 0; i < nfields; ++i)
    {
        NekDouble integral = m_f->m_exp[0]->Integral(m_f->m_exp[i]->GetPhys());
        if (m_f->m_comm->GetRank() == 0)
        {
            cout << "Integral (variable " << m_f->m_variables[i]
                 << ") : " << integral << endl;
            cout << "Mean (variable " << m_f->m_variables[i]
                 << ") : " << integral / scale << endl;
        }
    }
}
}
}

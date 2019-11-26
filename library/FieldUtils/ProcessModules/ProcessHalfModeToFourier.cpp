////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessHalfModeToFourier.cpp
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
//  Description: Take a FourierHalfMode expansion and correct so it can be
//  understood as a two mode Fourier Expansion
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>
using namespace std;

#include <boost/core/ignore_unused.hpp>

#include <LibUtilities/BasicUtils/SharedArray.hpp>

#include "ProcessHalfModeToFourier.h"

namespace Nektar
{
namespace FieldUtils
{

ModuleKey ProcessHalfModeToFourier::className =
    GetModuleFactory().RegisterCreatorFunction(
            ModuleKey(eProcessModule, "halfmodetofourier"),
            ProcessHalfModeToFourier::create,
            "modify a FourierHalfMode into a Fourier expansion so it can be "
            "processed as a 3D field.");

ProcessHalfModeToFourier::ProcessHalfModeToFourier(FieldSharedPtr f)
    : ProcessModule(f)
{
    m_priority = eModifyFieldData;
}

ProcessHalfModeToFourier::~ProcessHalfModeToFourier()
{
}

void ProcessHalfModeToFourier::Process(po::variables_map &vm)
{
    boost::ignore_unused(vm);

    // modify field definition
    for (int i = 0; i < m_f->m_data.size(); ++i)
    {
        ASSERTL0((m_f->m_fielddef[i]->m_basis[2] ==
                  LibUtilities::eFourierHalfModeRe ||
                  m_f->m_fielddef[i]->m_basis[2] ==
                  LibUtilities::eFourierHalfModeIm),
                 "This module is only for fourier Half modes");

        // change HomogeneousID
        m_f->m_fielddef[i]->m_homogeneousZIDs[0] += 2;

        // change expansion
        m_f->m_fielddef[i]->m_numModes[2] = 4;

        // Set third expansion to Fourier
        m_f->m_fielddef[i]->m_basis[2] = LibUtilities::eFourier;
    }

}
}
}

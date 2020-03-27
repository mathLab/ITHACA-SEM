////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessScaleInFld.cpp
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
//  Description: Scale input fld
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>
using namespace std;

#include <boost/core/ignore_unused.hpp>

#include <LibUtilities/BasicUtils/SharedArray.hpp>

#include "ProcessScaleInFld.h"

namespace Nektar
{
namespace FieldUtils
{

ModuleKey ProcessScaleInFld::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "scaleinputfld"),
        ProcessScaleInFld::create,
        "rescale input field by a constant factor.");

ProcessScaleInFld::ProcessScaleInFld(FieldSharedPtr f) : ProcessModule(f)
{
    m_config["scale"] = ConfigOption(false, "NotSet", "scale factor");
    ASSERTL0(m_config["scale"].as<string>().compare("NotSet") != 0,
             "scaleinputfld: Need to specify a scale factor");
}

ProcessScaleInFld::~ProcessScaleInFld()
{
}

void ProcessScaleInFld::Process(po::variables_map &vm)
{
    boost::ignore_unused(vm);

    string scalestr = m_config["scale"].as<string>();
    NekDouble scale = boost::lexical_cast<NekDouble>(scalestr);

    for (int i = 0; i < m_f->m_data.size(); ++i)
    {
        int datalen = m_f->m_data[i].size();

        Vmath::Smul(datalen, scale, &(m_f->m_data[i][0]), 1,
                    &(m_f->m_data[i][0]), 1);
    }
}
}
}

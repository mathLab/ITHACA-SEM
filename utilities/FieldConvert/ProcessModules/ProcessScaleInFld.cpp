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
//  Description: Scale input fld
//
////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <iostream>
using namespace std;

#include "ProcessScaleInFld.h"

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/ParseUtils.hpp>

namespace Nektar
{
namespace Utilities
{

ModuleKey ProcessScaleInFld::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "scaleinputfld"),
        ProcessScaleInFld::create, "rescale input field by a constant factor.");

ProcessScaleInFld::ProcessScaleInFld(FieldSharedPtr f) : ProcessModule(f)
{
    if((f->m_inputfiles.count("fld") == 0) &&
       (f->m_inputfiles.count("rst") == 0) &&
       (f->m_inputfiles.count("chk") == 0))
    {
        cout << "A fld, chk or rst input file must be specified for the "
                "scaleinputfld module" << endl;
        exit(3);
    }

    m_config["scale"] = ConfigOption(false,"NotSet","scale factor");
    ASSERTL0(m_config["scale"].as<string>().compare("NotSet") != 0,
             "scaleinputfld: Need to specify a sacle factor");
}

ProcessScaleInFld::~ProcessScaleInFld()
{
}

void ProcessScaleInFld::Process(po::variables_map &vm)
{
    if (m_f->m_verbose)
    {
        cout << "ProcessScaleInFld: Rescaling input fld" << endl;
    }

    ASSERTL0(m_f->m_data.size() != 0,"No input data defined");

    string scalestr = m_config["scale"].as<string>();
    NekDouble scale = boost::lexical_cast<NekDouble>(scalestr);

    for(int i = 0; i < m_f->m_data.size(); ++i)
    {
        int datalen = m_f->m_data[i].size();

        Vmath::Smul(datalen, scale, &(m_f->m_data[i][0]), 1,
                                    &(m_f->m_data[i][0]), 1);
    }

    if(m_f->m_exp.size())// expansiosn are defined reload field
    {
        int nfields = m_f->m_fielddef[0]->m_fields.size();

        // import basic field again in case of rescaling
        for (int j = 0; j < nfields; ++j)
        {
            for (int i = 0; i < m_f->m_data.size(); ++i)
            {
                m_f->m_exp[j]->ExtractDataToCoeffs(
                                       m_f->m_fielddef[i],
                                       m_f->m_data[i],
                                       m_f->m_fielddef[i]->m_fields[j],
                                       m_f->m_exp[j]->UpdateCoeffs());
            }
        }
    }
}

}
}

////////////////////////////////////////////////////////////////////////////////
//
//  File: Concatenate field
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
//  Description: Concatenate parallel field
//
////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <iostream>
using namespace std;

#include "ProcessConcatenateFld.h"


#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/ParseUtils.hpp>

namespace Nektar
{
namespace Utilities
{
ModuleKey ProcessConcatenateFld::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "concatenate"),
        ProcessConcatenateFld::create,
        "Concatenate field file into single file");

ProcessConcatenateFld::ProcessConcatenateFld(FieldSharedPtr f)
        : ProcessModule(f)
{
    // check for correct input files
    if((f->m_inputfiles.count("xml")    == 0) &&
       (f->m_inputfiles.count("xml.gz") == 0))
    {
        cout << "An xml or xml.gz input file must be specified for the "
                "concatenate module" << endl;
        exit(3);
    }

    if((f->m_inputfiles.count("fld") == 0) &&
       (f->m_inputfiles.count("chk") == 0) &&
       (f->m_inputfiles.count("rst") == 0))
    {
        cout << "A fld or chk or rst input file must be specified for the "
                "concatenate module" << endl;

        exit(3);
    }

}

ProcessConcatenateFld::~ProcessConcatenateFld()
{
}

void ProcessConcatenateFld::Process(po::variables_map &vm)
{
    if (m_f->m_verbose)
    {
        cout << "ProcessConcatenateFld: Concatenating field file" << endl;
    }

    std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef
        = m_f->m_exp[0]->GetFieldDefinitions();
    std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());

    // Copy Data into FieldData and set variable
    for(int j = 0; j < m_f->m_exp.size(); ++j)
    {
        for(int i = 0; i < FieldDef.size(); ++i)
        {
            // Could do a search here to find correct variable
            FieldDef[i]->m_fields.push_back(m_f->m_fielddef[0]->m_fields[j]);
            m_f->m_exp[0]->AppendFieldData(FieldDef[i], FieldData[i],
                                           m_f->m_exp[j]->UpdateCoeffs());
        }
    }

    m_f->m_fielddef  = FieldDef;
    m_f->m_data = FieldData;
}
}
}



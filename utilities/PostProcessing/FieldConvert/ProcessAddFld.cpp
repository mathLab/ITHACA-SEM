////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessAddFld.cpp
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
//  Description: Add a field to the intput field
//
////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <iostream>
using namespace std;

#include "ProcessAddFld.h"

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/ParseUtils.hpp>

namespace Nektar
{
namespace Utilities
{
ModuleKey ProcessAddFld::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "addfld"),
        ProcessAddFld::create, "rescale input field by a constant factor.");

ProcessAddFld::ProcessAddFld(FieldSharedPtr f) : ProcessModule(f)
{
    if((f->m_inputfiles.count("fld") == 0) &&
       (f->m_inputfiles.count("rst") == 0) &&
       (f->m_inputfiles.count("chk") == 0))
    {
        cout << "A fld, chk or rst input file must be specified for the "
                "scaleinputfld module" << endl;
        exit(3);
    }

    m_config["scale"]   = ConfigOption(false, "1.0", "scale factor");

    m_config["fromfld"] = ConfigOption(false, "NotSet",
                                "Fld file form which to interpolate field");

    ASSERTL0(m_config["fromfld"].as<string>().compare("NotSet") != 0,
             "Need to specify fromfld=file.fld ");

}

ProcessAddFld::~ProcessAddFld()
{
}

void ProcessAddFld::Process(po::variables_map &vm)
{
    if (m_f->m_verbose)
    {
        cout << "ProcessAddFld: Adding new fld to input fld" << endl;
    }

    ASSERTL0(m_f->m_data.size() != 0,"No input data defined");

    string scalestr = m_config["scale"].as<string>();
    NekDouble scale = boost::lexical_cast<NekDouble>(scalestr);

    string fromfld = m_config["fromfld"].as<string>();
    m_fromField =  boost::shared_ptr<Field>(new Field());

    if(m_f->m_exp.size())
    {
        // Set up ElementGIDs in case of parallel processing
        Array<OneD,int> ElementGIDs(m_f->m_exp[0]->GetExpSize());
        for (int i = 0; i < m_f->m_exp[0]->GetExpSize(); ++i)
        {
            ElementGIDs[i] = m_f->m_exp[0]->GetExp(i)->GetGeom()->GetGlobalID();
        }
        m_f->m_fld->Import(fromfld,m_fromField->m_fielddef,
                           m_fromField->m_data,
                           LibUtilities::NullFieldMetaDataMap,
                           ElementGIDs);
    }
    else
    {
        m_f->m_fld->Import(fromfld,m_fromField->m_fielddef,
                           m_fromField->m_data,
                           LibUtilities::NullFieldMetaDataMap);
    }

    bool samelength = true;
    if(m_fromField->m_data.size() != m_f->m_data.size())
    {
        samelength = false;
    }

    // scale input field
    for(int i = 0; i < m_fromField->m_data.size(); ++i)
    {
        int datalen = m_fromField->m_data[i].size();

        Vmath::Smul(datalen, scale, &(m_fromField->m_data[i][0]), 1,
                                    &(m_fromField->m_data[i][0]), 1);

        if(samelength)
        {
            if(datalen != m_f->m_data[i].size())
            {
                samelength = false;
            }
        }
    }

    if(samelength == true)
    {

        for(int i = 0; i < m_f->m_data.size(); ++i)
        {
            int datalen = m_f->m_data[i].size();

            Vmath::Vadd(datalen, &(m_f->m_data[i][0]),         1,
                                 &(m_fromField->m_data[i][0]), 1,
                                 &(m_f->m_data[i][0]),         1);
        }
    }
    else
    {
        ASSERTL0(m_f->m_exp.size() != 0 ,
                 "Input fields have partitions of different length and so xml "
                 "file needs to be specified");

        int nfields = m_f->m_fielddef[0]->m_fields.size();
        int ncoeffs = m_f->m_exp[0]->GetNcoeffs();
        Array<OneD, NekDouble> SaveFld(ncoeffs);

        for (int j = 0; j < nfields; ++j)
        {
            Vmath::Vcopy(ncoeffs,m_f->m_exp[j]->GetCoeffs(),1, SaveFld,1);

            // load new field
            for (int i = 0; i < m_fromField->m_data.size(); ++i)
            {
                m_f->m_exp[j]->ExtractDataToCoeffs(
                                       m_fromField->m_fielddef[i],
                                       m_fromField->m_data[i],
                                       m_fromField->m_fielddef[i]->m_fields[j],
                                       m_f->m_exp[j]->UpdateCoeffs());
            }

            Vmath::Vadd(ncoeffs, m_f->m_exp[j]->GetCoeffs(),    1,
                                 SaveFld,                       1,
                                 m_f->m_exp[j]->UpdateCoeffs(), 1);
        }

        std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef
            = m_f->m_exp[0]->GetFieldDefinitions();
        std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());

        for(int i = 0; i < nfields; ++i)
        {
            for (int j = 0; j < FieldDef.size(); ++j)
            {
                FieldDef[j]->m_fields.push_back(m_f->m_fielddef[0]->m_fields[i]);
                m_f->m_exp[i]->AppendFieldData(FieldDef[j], FieldData[j]);
            }
        }

        m_f->m_fielddef = FieldDef;
        m_f->m_data     = FieldData;

    }






}

}
}

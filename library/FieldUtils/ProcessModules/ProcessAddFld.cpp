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

#include <iostream>
#include <string>
using namespace std;

#include <boost/core/ignore_unused.hpp>

#include <LibUtilities/BasicUtils/SharedArray.hpp>

#include "ProcessAddFld.h"

namespace Nektar
{
namespace FieldUtils
{

ModuleKey ProcessAddFld::className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eProcessModule, "addfld"),
    ProcessAddFld::create,
    "add two fields together with optional scaling. Must specify fromfld and "
    "scaling is optionally specified with input option scale.");

ProcessAddFld::ProcessAddFld(FieldSharedPtr f) : ProcessModule(f)
{
    m_config["scale"] = ConfigOption(false, "1.0", "scale factor");

    m_config["fromfld"] =
        ConfigOption(false, "NotSet", "Fld file form which to add field");

    if(f->m_inputfiles.count("xml"))
    {
        m_priority = eModifyExp;
    }
    else
    {
        m_priority = eModifyFieldData;
    }
}

ProcessAddFld::~ProcessAddFld()
{
}

void ProcessAddFld::Process(po::variables_map &vm)
{
    boost::ignore_unused(vm);

    string scalestr = m_config["scale"].as<string>();
    NekDouble scale = boost::lexical_cast<NekDouble>(scalestr);

    ASSERTL0(m_config["fromfld"].as<string>().compare("NotSet") != 0,
             "Need to specify fromfld=file.fld ");
    string fromfld           = m_config["fromfld"].as<string>();

    vector<LibUtilities::FieldDefinitionsSharedPtr> fromFieldDef;
    vector<vector<double> >                         fromFieldData;

    if (m_f->m_graph)
    {
        const SpatialDomains::ExpansionMap &expansions =
            m_f->m_graph->GetExpansions();

        // if Range has been speficied it is possible to have a
        // partition which is empty so check this and return if
        // no elements present.

        if (!expansions.size())
        {
            return;
        }

        Array<OneD, int> ElementGIDs(expansions.size());

        int i = 0;
        for (auto &expIt : expansions)
        {
            ElementGIDs[i++] = expIt.second->m_geomShPtr->GetGlobalID();
        }
        m_f->FieldIOForFile(fromfld)->Import(
            fromfld, fromFieldDef, fromFieldData,
            LibUtilities::NullFieldMetaDataMap, ElementGIDs);
    }
    else
    {
        m_f->FieldIOForFile(fromfld)->Import(
            fromfld, fromFieldDef, fromFieldData,
            LibUtilities::NullFieldMetaDataMap);
    }

    bool samelength = true;
    if (fromFieldData.size() != m_f->m_data.size())
    {
        samelength = false;
    }

    // scale input field
    for (int i = 0; i < fromFieldData.size(); ++i)
    {
        int datalen = fromFieldData[i].size();

        Vmath::Smul(datalen, scale, &(fromFieldData[i][0]), 1,
                    &(fromFieldData[i][0]), 1);

        if (samelength)
        {
            if (datalen != m_f->m_data[i].size())
            {
                samelength = false;
            }
        }
    }

    if (m_priority == eModifyFieldData)
    {
        ASSERTL0(samelength == true,
                "Input fields have partitions of different length and so xml "
                 "file needs to be specified");
        for (int i = 0; i < m_f->m_data.size(); ++i)
        {
            int datalen = m_f->m_data[i].size();

            Vmath::Vadd(datalen, &(m_f->m_data[i][0]), 1,
                        &(fromFieldData[i][0]), 1, &(m_f->m_data[i][0]), 1);
        }
        
    }
    else
    {
        // Skip in case of empty partition
        if (m_f->m_exp[0]->GetNumElmts() == 0)
        {
            return;
        }

        int nfields = m_f->m_variables.size();
        int ncoeffs = m_f->m_exp[0]->GetNcoeffs();
        Array<OneD, NekDouble> SaveFld(ncoeffs);

        for (int j = 0; j < nfields; ++j)
        {
            Vmath::Vcopy(ncoeffs, m_f->m_exp[j]->GetCoeffs(), 1, SaveFld, 1);

            // Check if new field has this variable
            auto it = find (fromFieldDef[0]->m_fields.begin(),
                            fromFieldDef[0]->m_fields.end(),
                            m_f->m_variables[j]);

            ASSERTL0(it != fromFieldDef[0]->m_fields.end(),
              "Could not find field " + m_f->m_variables[j] + " in from field");

            // load new field
            for (int i = 0; i < fromFieldData.size(); ++i)
            {
                m_f->m_exp[j]->ExtractDataToCoeffs(
                    fromFieldDef[i], fromFieldData[i],
                    m_f->m_variables[j],
                    m_f->m_exp[j]->UpdateCoeffs());
            }

            Vmath::Vadd(ncoeffs, m_f->m_exp[j]->GetCoeffs(), 1, SaveFld, 1,
                        m_f->m_exp[j]->UpdateCoeffs(), 1);
            m_f->m_exp[j]->BwdTrans(
                        m_f->m_exp[j]->GetCoeffs(),
                        m_f->m_exp[j]->UpdatePhys());
        }
    }
}
}
}

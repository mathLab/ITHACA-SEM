////////////////////////////////////////////////////////////////////////////////
//
//  File: InputFld.cpp
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
//  Description: Reads a Nektar++ FLD file.
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>
using namespace std;

#include <boost/core/ignore_unused.hpp>

#include "InputFld.h"
using namespace Nektar;

namespace Nektar
{
namespace FieldUtils
{

ModuleKey InputFld::m_className[4] = {
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eInputModule, "fld"),
        InputFld::create, "Reads Fld file."),
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eInputModule, "chk"),
        InputFld::create, "Reads checkpoint file."),
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eInputModule, "rst"),
        InputFld::create, "Reads restart file."),
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eInputModule, "bse"),
        InputFld::create, "Reads stability base-flow file.")
};

/**
 * @brief Set up InputFld object.
 *
 */
InputFld::InputFld(FieldSharedPtr f) : InputModule(f)
{
    m_allowedFiles.insert("fld");
    m_allowedFiles.insert("chk");
    m_allowedFiles.insert("rst");
    m_allowedFiles.insert("bse");
}

/**
 *
 */
InputFld::~InputFld()
{
}

/**
 *
 */
void InputFld::Process(po::variables_map &vm)
{
    boost::ignore_unused(vm);

    int i;
    string fileName = m_config["infile"].as<string>();

    LibUtilities::FieldIOSharedPtr fld =
                m_f->FieldIOForFile(fileName);

    int oldSize = m_f->m_fielddef.size();
    if(m_f->m_graph)
    {
        // Determine IDs of elements in the domain
        vector<int> IDs;
        auto domain = m_f->m_graph->GetDomain();
        for(int d = 0; d < domain.size(); ++d)
        {
            for (auto &compIter : domain[d])
            {
                for (auto &x : compIter.second->m_geomVec)
                {
                    IDs.push_back(x->GetGlobalID());
                }
            }
        }

        if (!IDs.size())
        {
            return;
        }

        // Move to an array to match FieldIO interface
        Array<OneD, int> ElementGIDs(IDs.size(), IDs.data());

        fld->Import(
            fileName, m_f->m_fielddef, m_f->m_data,
            m_f->m_fieldMetaDataMap, ElementGIDs);
    }
    else // load all data.
    {
        fld->Import(
            fileName, m_f->m_fielddef, m_f->m_data,
            m_f->m_fieldMetaDataMap);
    }

    // save field names
    for(i = 0; i < m_f->m_fielddef[oldSize]->m_fields.size(); ++i)
    {
        // check for multiple fld files
        auto it = find (m_f->m_variables.begin(),
                        m_f->m_variables.end(),
                        m_f->m_fielddef[oldSize]->m_fields[i]);

        if(it == m_f->m_variables.end())
        {
            m_f->m_variables.push_back(m_f->m_fielddef[oldSize]->m_fields[i]);
        }
    }
}
}
}

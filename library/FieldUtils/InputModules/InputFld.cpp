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
//  Description: Reads a Nektar++ FLD file.
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>
using namespace std;

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
    int i;
    string fldending;
    // Determine appropriate field input
    if (m_f->m_inputfiles.count("fld") != 0)
    {
        fldending = "fld";
    }
    else if (m_f->m_inputfiles.count("chk") != 0)
    {
        fldending = "chk";
    }
    else if (m_f->m_inputfiles.count("rst") != 0)
    {
        fldending = "rst";
    }
    else if (m_f->m_inputfiles.count("bse") != 0)
    {
        fldending = "bse";
    }
    else
    {
        ASSERTL0(false, "no input file found");
    }

    if(m_f->m_graph)
    {
        if (m_f->m_data.size() == 0)
        {
            // currently load all field (possibly could read data from
            //  expansion list but it is re-arranged in expansion)

            const SpatialDomains::ExpansionMap &expansions =
                m_f->m_graph->GetExpansions();

            // if Range has been specified it is possible to have a
            // partition which is empty so check this and return if
            // no elements present.

            if (!expansions.size())
            {
                return;
            }

            Array<OneD, int> ElementGIDs(expansions.size());
            SpatialDomains::ExpansionMap::const_iterator expIt;

            i = 0;
            for (expIt = expansions.begin(); expIt != expansions.end(); ++expIt)
            {
                ElementGIDs[i++] = expIt->second->m_geomShPtr->GetGlobalID();
            }

            m_f->m_fielddef.clear();
            m_f->m_data.clear();

            m_f->FieldIOForFile(m_f->m_inputfiles[fldending][0])->Import(
                m_f->m_inputfiles[fldending][0], m_f->m_fielddef, m_f->m_data,
                m_f->m_fieldMetaDataMap, ElementGIDs);
        }
    }
    else // load all data.
    {
        m_f->FieldIOForFile(m_f->m_inputfiles[fldending][0])->Import(
            m_f->m_inputfiles[fldending][0], m_f->m_fielddef, m_f->m_data,
            m_f->m_fieldMetaDataMap);
    }

    // save field names
    m_f->m_variables = m_f->m_fielddef[0]->m_fields;
}
}
}

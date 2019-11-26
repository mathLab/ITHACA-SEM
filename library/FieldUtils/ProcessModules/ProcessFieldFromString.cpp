///////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessFieldFromString.cpp
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
//  Description: Modify an existing or add a new field from a string based on existing variable
//
///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <string>
using namespace std;

#include <boost/core/ignore_unused.hpp>

#include <LibUtilities/BasicUtils/ParseUtils.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>

#include "ProcessFieldFromString.h"

namespace Nektar
{
namespace FieldUtils
{

ModuleKey ProcessFieldFromString::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "fieldfromstring"),
        ProcessFieldFromString::create,
        "Modify an existing or create a new field from the existing fields as "
        "specified by a string using a required argument of the form "
        "fieldstr=\"x + y + u\" ");

ProcessFieldFromString::ProcessFieldFromString(FieldSharedPtr f)
    : ProcessModule(f)
{
    m_config["fieldstr"] = ConfigOption(
        false, "NotSet", "Analytic expression");
    m_config["fieldname"] =
        ConfigOption(false,
                     "newfield",
                     "name for modified new field, default is \"newfield\" (optional)");
}

ProcessFieldFromString::~ProcessFieldFromString(void)
{
}

void ProcessFieldFromString::Process(po::variables_map &vm)
{
    boost::ignore_unused(vm);

    // Check if required parameter fieldstr was provided
    ASSERTL0(m_config["fieldstr"].m_beenSet, "fieldstr must be specified");

    // Get number of fields (before adding new entry)
    int nfields = m_f->m_variables.size();

    // Set up new field name
    string fieldName = m_config["fieldname"].as<string>();

    int fieldID;
    bool addField;
    // check if field exists
    auto it =
        std::find(m_f->m_variables.begin(), m_f->m_variables.end(), fieldName);
    if (it != m_f->m_variables.end())
    {
        addField = false;
        fieldID = std::distance(m_f->m_variables.begin(), it);
    }
    else
    {
        // Create new expansion
        addField = true;
        fieldID  = nfields;
        m_f->m_variables.push_back(fieldName);
    }

    // Skip in case of empty partition
    if (m_f->m_exp[0]->GetNumElmts() == 0)
    {
        return;
    }

    // Check if using strips
    int nstrips;
    m_f->m_session->LoadParameter("Strip_Z", nstrips, 1);
    ASSERTL0(nstrips == 1,
             "Routine is currently only setup for non-strip files");

    if (addField)
    {
        m_f->m_exp.resize(nfields + 1);
        m_f->m_exp[nfields] = m_f->AppendExpList(m_f->m_numHomogeneousDir);
    }

    // Variables for storing names and values for evaluating the function
    string varstr;
    vector<Array<OneD, const NekDouble> > interpfields;

    // Add the coordinate values
    varstr += "x y z";
    int npoints = m_f->m_exp[0]->GetTotPoints();
    Array<OneD, NekDouble> x(npoints, 0.0);
    Array<OneD, NekDouble> y(npoints, 0.0);
    Array<OneD, NekDouble> z(npoints, 0.0);
    m_f->m_exp[0]->GetCoords(x, y, z);
    interpfields.push_back(x);
    interpfields.push_back(y);
    interpfields.push_back(z);

    // Add the field values
    for (int i = 0; i < nfields; ++i)
    {
        varstr += " " + m_f->m_variables[i];
        interpfields.push_back(m_f->m_exp[i]->GetPhys());
    }

    // Create new function
    LibUtilities::Interpreter strEval;
    int exprId      = -1;
    string fieldstr = m_config["fieldstr"].as<string>();
    exprId          = strEval.DefineFunction(varstr.c_str(), fieldstr);

    // Evaluate function
    strEval.Evaluate(exprId, interpfields, m_f->m_exp[fieldID]->UpdatePhys());

    // Update coeffs
    m_f->m_exp[fieldID]->FwdTrans_IterPerExp(
        m_f->m_exp[fieldID]->GetPhys(), m_f->m_exp[fieldID]->UpdateCoeffs());
}
}
}

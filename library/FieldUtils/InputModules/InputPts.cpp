////////////////////////////////////////////////////////////////////////////////
//
//  File: InputPts.cpp
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
//  Description: Read xml file of a series of points and hold
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>
using namespace std;

#include <boost/core/ignore_unused.hpp>

#include <LibUtilities/BasicUtils/PtsField.h>
#include <LibUtilities/BasicUtils/PtsIO.h>
#include <LibUtilities/BasicUtils/CsvIO.h>

#include <tinyxml.h>

#include "InputPts.h"

namespace Nektar
{
namespace FieldUtils
{

ModuleKey InputPts::m_className[5] = {
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eInputModule, "pts"), InputPts::create, "Reads Pts file."),
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eInputModule, "pts.gz"), InputPts::create, "Reads Pts file."),
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eInputModule, "csv"), InputPts::create, "Reads csv file."),
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eInputModule, "csv.gz"), InputPts::create, "Reads csv file."),
};

/**
 * @brief Set up InputPts object.
 *
 */
InputPts::InputPts(FieldSharedPtr f) : InputModule(f)
{
    m_allowedFiles.insert("pts");
    m_allowedFiles.insert("csv");
}

/**
 *
 */
InputPts::~InputPts()
{
}

/**
 *
 */
void InputPts::Process(po::variables_map &vm)
{
    boost::ignore_unused(vm);

    string inFile = m_config["infile"].as<string>();

    // Determine appropriate field input
    if (m_f->m_inputfiles.count("csv") != 0)
    {
        LibUtilities::CsvIOSharedPtr csvIO =
            MemoryManager<LibUtilities::CsvIO>::AllocateSharedPtr(m_f->m_comm);
        csvIO->Import(inFile, m_f->m_fieldPts);
    }
    else if (m_f->m_inputfiles.count("pts") != 0)
    {
        LibUtilities::PtsIOSharedPtr ptsIO =
            MemoryManager<LibUtilities::PtsIO>::AllocateSharedPtr(m_f->m_comm);
        ptsIO->Import(inFile, m_f->m_fieldPts);
    }
    else
    {
        ASSERTL0(false, "unknown input file type");
    }

    // save field names
    for (int j = 0; j < m_f->m_fieldPts->GetNFields(); ++j)
    {
        m_f->m_variables.push_back(m_f->m_fieldPts->GetFieldName(j));
    }
}
}
}

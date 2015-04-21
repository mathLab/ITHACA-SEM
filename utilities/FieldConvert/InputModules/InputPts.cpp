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
//  Description: Read xml file of a series of points and hold
//
////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <iostream>
using namespace std;

#include <LibUtilities/BasicUtils/PtsIO.h>
#include <LibUtilities/BasicUtils/PtsField.h>

#include <tinyxml.h>

#include "InputPts.h"

namespace Nektar
{
namespace Utilities
{

ModuleKey InputPts::m_className[5] = {
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eInputModule, "pts"), InputPts::create, "Reads Pts file."),
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eInputModule, "pts.gz"), InputPts::create, "Reads Pts file."),
};


/**
 * @brief Set up InputPts object.
 *
 */
InputPts::InputPts(FieldSharedPtr f) : InputModule(f)
{
    m_allowedFiles.insert("pts");
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
    if (m_f->m_verbose)
    {
        cout << "Processing input pts file" << endl;
    }

    string inFile = (m_f->m_inputfiles["pts"][0]).c_str();

    m_f->m_ptsIO = MemoryManager<LibUtilities::PtsIO>::AllocateSharedPtr();
    m_f->m_ptsIO->Import(inFile,  m_f->m_fieldPts);
}

}
}

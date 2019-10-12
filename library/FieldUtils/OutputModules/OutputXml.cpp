////////////////////////////////////////////////////////////////////////////////
//
//  File: OutputXml.cpp
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
//  Description: VTK file format output.
//
////////////////////////////////////////////////////////////////////////////////

#include <set>
#include <string>
using namespace std;

#include <boost/core/ignore_unused.hpp>

#include "OutputXml.h"

namespace Nektar
{
namespace FieldUtils
{

ModuleKey OutputXml::m_className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eOutputModule, "xml"), OutputXml::create, "Writes an XML file.");

OutputXml::OutputXml(FieldSharedPtr f) : OutputModule(f)
{
}

OutputXml::~OutputXml()
{
}

void OutputXml::Process(po::variables_map &vm)
{
    boost::ignore_unused(vm);

    if (!m_f->m_exp.size()) // do nothing if no expansion defined
    {
        return;
    }

    // Extract the output filename and extension
    string filename = m_config["outfile"].as<string>();

    m_f->m_graph->WriteGeometry(filename);
    cout << "Written file: " << filename << endl;
}
}
}

////////////////////////////////////////////////////////////////////////////////
//
//  File: OutputXml.h
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
//  Description: Vtk output module
//
////////////////////////////////////////////////////////////////////////////////

#ifndef FIELDUTILS_OUTPUTXML
#define FIELDUTILS_OUTPUTXML

#include "../Module.h"
#include <tinyxml.h>

namespace Nektar
{
namespace FieldUtils
{
/// Converter from fld to vtk.
class OutputXml : public OutputModule
{
public:
    /// Creates an instance of this class
    static std::shared_ptr<Module> create(FieldSharedPtr f)
    {
        return MemoryManager<OutputXml>::AllocateSharedPtr(f);
    }
    static ModuleKey m_className;

    OutputXml(FieldSharedPtr f);
    virtual ~OutputXml();

    /// Write fld to output file.
    virtual void Process(po::variables_map &vm);

    virtual std::string GetModuleName()
    {
        return "OutputXml";
    }

    virtual std::string GetModuleDescription()
    {
        return "Writing file";
    }

    virtual ModulePriority GetModulePriority()
    {
        return eOutput;
    }

};
}
}

#endif

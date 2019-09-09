////////////////////////////////////////////////////////////////////////////////
//
//  File: InputPts.h
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
//  Description: Xml reader and expansion setup
//
////////////////////////////////////////////////////////////////////////////////

#ifndef FIELDUTILS_INPUTPTS
#define FIELDUTILS_INPUTPTS

#include "../Module.h"

namespace Nektar
{
namespace FieldUtils
{

/**
 *  Input module for Xml files.
 */
class InputPts : public InputModule
{
public:
    InputPts(FieldSharedPtr f);
    virtual ~InputPts();
    virtual void Process(po::variables_map &vm);

    /// Creates an instance of this class
    static ModuleSharedPtr create(FieldSharedPtr f)
    {
        return MemoryManager<InputPts>::AllocateSharedPtr(f);
    }
    /// %ModuleKey for class.
    static ModuleKey m_className[];

    virtual std::string GetModuleName()
    {
        return "InputPts";
    }

    virtual std::string GetModuleDescription()
    {
        return "Processing input pts file";
    }

    virtual ModulePriority GetModulePriority()
    {
        return eCreatePts;
    }

private:
};
}
}

#endif

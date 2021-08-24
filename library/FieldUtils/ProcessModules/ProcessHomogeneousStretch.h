////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessHomogeneousStretch.h
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
//  Description: Stretch homogeneous direction by an integer factor.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef FIELDUTILS_PROCESSHOMOGENEOUSSTRETCH
#define FIELDUTILS_PROCESSHOMOGENEOUSSTRETCH

#include "../Module.h"

namespace Nektar
{
namespace FieldUtils
{

/**
 * @brief This processing module stretches the homogeneous direction of a
 * 3DH1D expansion by an integer factor.
 */
class ProcessHomogeneousStretch : public ProcessModule
{
public:
    /// Creates an instance of this class
    static std::shared_ptr<Module> create(FieldSharedPtr f)
    {
        return MemoryManager<ProcessHomogeneousStretch>::AllocateSharedPtr(f);
    }
    static ModuleKey className;

    ProcessHomogeneousStretch(FieldSharedPtr f);
    virtual ~ProcessHomogeneousStretch();

    /// Write mesh to output file.
    virtual void Process(po::variables_map &vm);

    virtual std::string GetModuleName()
    {
        return "ProcessHomogeneousStretch";
    }

    virtual std::string GetModuleDescription()
    {
        return "Stretching expansion";
    }

    virtual ModulePriority GetModulePriority()
    {
        return eModifyExp;
    }

};
}
}

#endif

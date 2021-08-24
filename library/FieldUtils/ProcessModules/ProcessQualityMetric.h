////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessJacobianEnergy.h
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
//  Description: Computes energy of Jacobian.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef FIELDUTILS_PROCESSQUALITYMETRIC
#define FIELDUTILS_PROCESSQUALITYMETRIC

#include "../Module.h"

namespace Nektar
{
namespace FieldUtils
{

/// This processing module scales the input fld file
class ProcessQualityMetric : public ProcessModule
{
public:
    /// Creates an instance of this class
    static std::shared_ptr<Module> create(FieldSharedPtr f)
    {
        return MemoryManager<ProcessQualityMetric>::AllocateSharedPtr(f);
    }
    static ModuleKey className;

    ProcessQualityMetric(FieldSharedPtr f);
    virtual ~ProcessQualityMetric();

    /// Write mesh to output file.
    virtual void Process(po::variables_map &vm);

    virtual std::string GetModuleName()
    {
        return "ProcessQualityMetric";
    }

    virtual std::string GetModuleDescription()
    {
        return "Adding quality metric to field";
    }

    virtual ModulePriority GetModulePriority()
    {
        return eModifyExp;
    }

private:
    Array<OneD, NekDouble> GetQ(LocalRegions::ExpansionSharedPtr e, bool s);
};
}
}

#endif

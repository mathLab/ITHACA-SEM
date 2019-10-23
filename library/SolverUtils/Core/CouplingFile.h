///////////////////////////////////////////////////////////////////////////////
//
// File CouplingFile.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2017 Kilian Lackhove
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: File based Coupling class. Rather pointless, only for
// demonstration.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_COUPLINGFILE
#define NEKTAR_COUPLINGFILE

#include <SolverUtils/Core/Coupling.h>
#include <SolverUtils/Core/SessionFunction.h>

namespace Nektar
{
namespace SolverUtils
{

class CouplingFile;

class CouplingFile : public Coupling
{

public:
    static std::string className;

    /// Creates an instance of this class
    SOLVER_UTILS_EXPORT static CouplingSharedPtr create(
        MultiRegions::ExpListSharedPtr field)
    {
        CouplingSharedPtr p =
            MemoryManager<CouplingFile>::AllocateSharedPtr(field);
        p->Init();
        return p;
    }

    SOLVER_UTILS_EXPORT CouplingFile(MultiRegions::ExpListSharedPtr field);

    SOLVER_UTILS_EXPORT virtual ~CouplingFile();

protected:
    SOLVER_UTILS_EXPORT virtual void v_Init();

    SOLVER_UTILS_EXPORT virtual void v_Send(
        const int step,
        const NekDouble time,
        const Array<OneD, const Array<OneD, NekDouble> > &field,
        std::vector<std::string> &varNames);

    SOLVER_UTILS_EXPORT virtual void v_Receive(
        const int step,
        const NekDouble time,
        Array<OneD, Array<OneD, NekDouble> > &field,
        std::vector<std::string> &varNames);

private:
    int m_lastSend;
    int m_lastReceive;

    SessionFunctionSharedPtr m_inputFunction;
};
}
}

#endif
